# Reference paper: "Being Robust in High Dimension Can Be Practical"
# by Diakonikolas et al: https://arxiv.org/pdf/1703.00893.pdf
#
# matlab repo https://github.com/hoonose/robust-filter (no license)
#
# XXX note: this is a sandbox to translate the matlab code and see if it makes
# sense

# DEV NOTE:
# - assumption 0-mean gaussian (therefore no centering)
# - assumption that n > p at the start

function filter_gcov(X::AbstractMatrix{<:Real};
                     eps::Float64=0.05,
                     tau::Float64=0.1,
                     maxcond::Int=10_000)

    n, p   = size(X)
    thresh = eps * log(1.0/eps)^2
    C1, C2 = 0.4, 0.0

    Xk = X

    Sk = nothing # cache for Sk and its square root

    while true
        nk = size(Xk, 1)
        Sk = cov(Xk, Simple(corrected=false); mean=0)
        κ  = cond(Sk)

        if (κ > maxcond) || (nk < p)
            @warn "Ill conditioned case $κ, dims ($nk, $p)"
            return zeros(p, p)
        end

        cholesky!(Sk) # now UpperTriangular(Sk)' is L and previous Sk = LL'
        xΩxt = sum(x->x^2, UpperTriangular(Sk)' \ Xk', dims=1) # size 1 x nk

        # is any further than threshold distance?
        mask = (xΩxt .> C1 * p * log(nk / tau))

        any(mask) || break

        # update Xk, dropping the "outliers"
        Xk = Xk[mask, :]
    end # end while
end
