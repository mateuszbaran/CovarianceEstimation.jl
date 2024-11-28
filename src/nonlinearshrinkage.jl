"""
    AnalyticalNonlinearShrinkage

Analytical nonlinear shrinkage estimator. See docs for
`analytical_nonlinear_shrinkage` for details.
"""
struct AnalyticalNonlinearShrinkage{TEigen<:Union{Eigen,Nothing},Talpha<:Real} <: CovarianceEstimator
    decomp::TEigen
    corrected::Bool
    alpha::Talpha
    warn_if_small_eigenvalue::Bool
    function AnalyticalNonlinearShrinkage(;
        corrected=false,
        alpha::Real=0.0,
        warn_if_small_eigenvalue::Bool=false,
    )
        new{Nothing,typeof(alpha)}(nothing, corrected, alpha, warn_if_small_eigenvalue)
    end
    function AnalyticalNonlinearShrinkage(decomp::Eigen;
        corrected=false,
        alpha::Real=0.0,
        warn_if_small_eigenvalue::Bool=false,
    )
        new{Eigen,typeof(alpha)}(decomp, corrected, alpha, warn_if_small_eigenvalue)
    end
end


const SQRT5  = sqrt(5.0)
const INVPI  = 1.0 / π
const EPAN_1 = 3.0/(4.0 * SQRT5)
const EPAN_2 = EPAN_1 * INVPI
const EPAN_3 = 0.3 * INVPI

"""
    epanechnikov(x)

Return the Epanechnikov kernel evaluated at `x`.
"""
epanechnikov(x::T) where T<:Real = float(T)(EPAN_1 * max(0.0, 1.0 - x^2/5.0))

"""
    epnanechnikov_HT(x)

Return the Hilbert Transform of the Epanechnikov kernel evaluated at `x`
if `|x|≂̸√5`.
"""
function epanechnikov_HT1(x::T) where T <: Real
    float(T)(-EPAN_3 * x + EPAN_2 * (1.0 - x^2/5.0) * log(abs((SQRT5 - x)/(SQRT5 + x))))
end

"""
    epnanechnikov_HT2(x)
Return the Hilbert Transform of the Epanechnikov kernel evaluated at `x`
if `|x|=√5`.
"""
epanechnikov_HT2(x::T) where T <: Real = float(T)(-EPAN_3*x)

"""
    analytical_nonlinear_shrinkage(S, n, p; decomp;
        alpha::Real=0.0,
        warn_if_small_eigenvalue::Bool=false,
    )

Internal implementation of the analytical nonlinear shrinkage. The
implementation is inspired from the Matlab code given in section C of
Olivier Ledoit and Michael Wolf's paper "Analytical Nonlinear
Shrinkage of Large-Dimensional Covariance Matrices". (Nov 2018)
http://www.econ.uzh.ch/static/wp/econwp264.pdf

Shrinked eigenvalues smaller than `alpha` are replaced with `alpha`.
If `warn_if_small_eigenvalue` is `true`, additionally a warning is emitted
if any eigenvalue is replaced with `alpha`.
"""
function analytical_nonlinear_shrinkage(S::AbstractMatrix{<:Real},
                                        n::Int, p::Int, est_mean::Bool;
                                        decomp::Union{Nothing,Eigen}=nothing,
                                        alpha::Real=0.0,
                                        warn_if_small_eigenvalue::Bool=false)

    η = ifelse(p < n, n, n - Int(est_mean)) # effective sample size
    # sample eigenvalues sorted in ascending order and eigenvectors
    F    = isa(decomp, Nothing) ? eigen(S) : decomp # O(p^3)
    perm = sortperm(F.values)
    sample_perm = @view perm[max(1, (p - η) + 1):p]
    λ    = @view F.values[sample_perm]
    U    = F.vectors[:, perm]
    T    = float(eltype(F))

    # dominant cost forming of S or eigen(S) --> O(max{np^2, p^3})

    # compute analytical nonlinear shrinkage kernel formula
    L = repeat(λ, outer=(1, min(p, η)))
    # Equation (4.9)
    h = T(η^(-1//3))
    H = h * L'
    x = (L .- L') ./ H

    # additional useful definitions
    γ  = T(p/η)
    πλ = π * λ

    # Equation (4.7) in http://www.econ.uzh.ch/static/wp/econwp264.pdf
    f̃ = vec(mean(epanechnikov.(x) ./ H, dims=2))

    # Equation (4.8)
    Hf̃_tmp = epanechnikov_HT1.(x)
    # if by any chance there are x such that |x| ≈ √5 ...
    mask = (@. abs(x) ≈ SQRT5)
    any(mask) && (Hf̃_tmp[mask] = epanechnikov_HT2.(@view x[mask]))
    Hf̃ = vec(mean(Hf̃_tmp ./ H, dims=2))

    # dominant cost up to here: elementwise ops on x --> O(max{p^2, η^2})

    if p <= η
        # Equation (4.3)
        πγλ = γ * πλ
        denom = @. (πγλ * f̃)^2 + (one(T) - γ - πγλ * Hf̃)^2
        d̃ = λ ./ denom
    else
        # Equation (C.8)
        hs5  = T(h * SQRT5)
        Hf̃0  = T((0.3/h^2 + 0.75/hs5 * (1.0 - 0.2/h^2) * log((1+hs5)/(1-hs5))))
        Hf̃0 *= mean(one(T) ./ πλ)
        # Equation (C.5)
        d̃0 = T(INVPI / ((γ - one(T)) * Hf̃0))
        # Eq. (C.4)
        d̃1 = @. one(T) / (π * πλ * (f̃^2 + Hf̃^2))
        d̃  = [fill(d̃0, (p - η, 1)); d̃1]
    end
    if alpha > 0
        for i = eachindex(d̃)
            if d̃[i] < alpha
                if warn_if_small_eigenvalue
                    @warn "analytical_nonlinear_shrinkage: Eingenvalue at position smaller than " i alpha
                end
                d̃[i] = alpha
            end
        end
    end

    # Equation (4.4)
    return Symmetric(U*(d̃ .* U'))
end

"""
    cov(ans::AnalyticalNonlinearShrinkage, X, [weights]; dims=1, mean=nothing)

Nonlinear covariance estimator derived from the sample covariance estimator `S`
and its eigenvalue decomposition (which can be given through `decomp`).
See Ledoit and Wolf's paper http://www.econ.uzh.ch/static/wp/econwp264.pdf
The keyword `mean` can be `nothing` (centering via estimated mean),
zero (no centering) or a provided vector. In the first case, a rank-1
modification is applied and therefore the effective sample size is decreased
by one (see `analytical_nonlinear_shrinkage`). In the latter two case the mean
cannot have been estimated on the data (otherwise the effective sample size
will be 1 larger than it should be resulting in numerical instabilities).
If you are unsure, use either `nothing` or provide an explicit (non-estimated)
vector (possibly a zero vector) and avoid the use of `mean=0`.

* Time complexity (including formation of `S`)
    - (p<n): O(np^2 + n^2) with moderate constant
    - (p>n): O(p^3) with low constant (dominated by eigendecomposition of S)
"""
function cov(ans::AnalyticalNonlinearShrinkage, X::AbstractMatrix{<:Real}, weights::FrequencyWeights...;
             dims::Int=1, mean=nothing)

    @assert dims ∈ [1, 2] "Argument dims can only be 1 or 2 (given: $dims)"

    szx    = size(X)
    (n, p) = ifelse(dims==1, szx, reverse(szx))
    wn = floor(Int, totalweight(n, weights...))

    # explained in the paper there must be at least 12 samples
    (wn < 12) && throw(ArgumentError("The (weighted) number of samples `n` must be at " *
                                    "least 12 (given: $wn)."))

    S  = cov(SimpleCovariance(corrected=ans.corrected), X, weights...; dims=dims, mean=mean)
    return analytical_nonlinear_shrinkage(S, wn, p, mean === nothing;
                decomp=ans.decomp,
                alpha=ans.alpha,
                warn_if_small_eigenvalue=ans.warn_if_small_eigenvalue)
end
