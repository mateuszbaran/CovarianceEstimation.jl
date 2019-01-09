"""
    Simple(corrected::Bool, weights::Union{AbstractWeights, Nothing} = nothing)

Simple covariance estimator (scaled by `n-1` if corrected and `n` otherwise
where `n` is the number of samples). Optionally supports weights (see
http://juliastats.github.io/StatsBase.jl/stable/cov.html).
"""
struct Simple{TWeights<:Union{AbstractWeights, Nothing}} <: CovarianceEstimator
    corrected::Bool
    weights::TWeights
end

Simple(;corrected::Bool = false) = Simple(corrected, nothing)
Simple(w::AbstractWeights; corrected::Bool = false) = Simple(corrected, w)


"""
    cov(x::AbstractVector, y::AbstractVector, sc::Simple)

Compute the covariance of the vectors `x` and `y` using formula
``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate. If `sc.corrected` is true then the fraction ``\\frac{1}{n}``
is replaced with ``\\frac{1}{n-1}``.
"""
cov(x::AbstractVector, y::AbstractVector, sc::Simple) =
    cov(x, y; corrected=sc.corrected)


"""
    cov(X::AbstractMatrix, sc::Simple; dims=1, mean=nothing)

Compute the sample variance of `X`. The sum is scaled with the number of
observations `n` if `sc.corrected` is false and with `n-1` otherwise.
If `dims=1` the rows are assumed to be observations and the columns the
features. If `dims=2` the opposite is assumed.
"""
function cov(X::AbstractMatrix, sc::Simple; dims::Int=1, mean=nothing)
    @assert dims ∈ [1, 2] "Argument dims can only be 1 or 2 (given: $dims)"

    # weighted case via StatsBase (no argument mean supported)
    sc.weights === nothing || return StatsBase.cov(X, sc.weights, dims;
                                        corrected=sc.corrected)
    # unweighted case via Statistics
    if mean === nothing
        return Statistics.cov(X; dims=dims, corrected=sc.corrected)
    elseif iszero(mean)
        return Statistics.covzm(X, dims, corrected=sc.corrected)
    elseif mean isa AbstractArray
        pdim = 2 - mod(dims+1, 2)
        length(mean) == size(X, pdim) || throw(DimensionMismatch(
                                "Provided mean length must match the " *
                                "dimensions of `X`. Got $(length(mean)), " *
                                "expected $(size(X, dims))."))
        return Statistics.covm(X, mean, dims; corrected=sc.corrected)
   else
       throw(ArgumentError("`mean` kw expects `0`, `nothing` or a vector."))
   end
end


"""
    cov(X::AbstractVector, sc::Simple; mean=nothing)

Compute the sample variance of `X`. The sum is scaled with the number of
observations `n` if `sc.corrected` is false and with `n-1` otherwise.
"""
function cov(X::AbstractVector, sc::Simple; mean=nothing)
    # weighted case via StatsBase (no argument mean supported)
    sc.weights === nothing || return StatsBase.cov(X, sc.weights;
                                        corrected=sc.corrected)
    # unweighted case via Statistics
    if mean === nothing
        return Statistics.cov(X; corrected=sc.corrected)
    elseif iszero(mean)
        return Statistics.covzm(X; corrected=sc.corrected)
    elseif mean isa Number
        return Statistics.covm(X, mean; corrected=sc.corrected)
   else
       throw(ArgumentError("`mean` kw expects `0`, `nothing` or a vector."))
   end
end
