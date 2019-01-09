# NOTE: these are older implementations/functions which were used in the past
# but have now been optimised heavily and ultimately replaced by other
# functions they're only kept here for the tests.

function sum_var_sij(Z::AbstractMatrix, CovZ::AbstractMatrix, n::Int,
                     corrected=false; with_diag=true)
    scale = ifelse(corrected, n-1, n)
    Z²    = Z.^2
    π̂mat  = ((Z²'*Z²)/n - CovZ.^2) / scale
    with_diag && return sum(π̂mat)
    return sum(π̂mat) - sum(diag(π̂mat))
end

function centercols(X::AbstractMatrix, μ=nothing; dims=1)
    Xc = (dims == 1) ? copy(X) : transpose(X)
    if μ === nothing
        Xc .-= mean(X, dims=1)
    elseif iszero(μ)
        # assume X is centered already
    elseif μ isa AbstractVector
        @assert length(μ) == size(Xc, 2) "Dimension mismatch, given mean has length $(length(μ)) while p=$(size(Xc, 2))"
        Xc .-= μ
    else
        throw(ArgumentError("Incorrect assignment for `mean`, expected either 0, nothing or a vector, got $μ."))
    end
    Xc
end
