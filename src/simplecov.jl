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
    cov(x::AbstractVector, sc::Simple)

Compute the sample variance of the vector `x`. The sum is scaled with `n`
where `n = length(x)` if `sc.corrected` is false and with `n-1` otherwise.
"""
cov(x::AbstractVector, sc::Simple) = cov(x; corrected = sc.corrected)


"""
    cov(x::AbstractVector, y, sc::Simple)

Compute the covariance of the vectors `x` and `y` using formula
``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate. If `sc.corrected` is true then the fraction ``\\frac{1}{n}``
is replaced with ``\\frac{1}{n-1}``.
"""
cov(x::AbstractVector, y::AbstractVector, sc::Simple) = cov(x, y; corrected = sc.corrected)


"""
    cov(X::AbstractMatrix, sc::Simple; dims::Int=1)

Compute the covariance matrix associated with `X` along dimension `dims`.
The sum is scaled with `n` where  `n = length(x)` if `sc.corrected` is false
and with `n-1` otherwise.
"""
cov(X::AbstractMatrix, sc::Simple; dims::Int=1) = cov(X; dims=dims, corrected = sc.corrected)

"""
    cov(X, sc::Simple{<:AbstractWeights}; dims=1)

Compute the weighted covariance matrix.
"""
function cov(X::DenseMatrix, sc::Simple{<:AbstractWeights}; dims::Int=1)
    return cov(X, sc.weights, dims, corrected = sc.corrected)
end
