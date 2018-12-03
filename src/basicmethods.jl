
"""
    Simple

Uncorrected covariance estimator (scaled by `n` where `n` is the number of
elements in a vector).
"""
struct Simple <: CovarianceEstimator
end

"""
    cov(::Simple, x::AbstractVector)

Compute variance of the vector `x`. The sum is scaled with `n`
where `n = length(x)`.
"""
cov(::Simple, x::AbstractVector) = cov(x; corrected = false)

"""
    cov(::Simple, x::AbstractVector, y)

Compute covariance of the vectors `x` and `y` using formula
``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate.
"""
cov(::Simple, x::AbstractVector, y::AbstractVector) = cov(x, y; corrected = false)

"""
    cov(::Simple, X::AbstractMatrix)

Compute the covariance matrix of the matrix `X` along dimension `dims`.
The sum is scaled with `n` where `n = length(x)`.
"""
cov(::Simple, X::AbstractMatrix; dims::Int=1) = cov(X; dims=dims, corrected = false)


"""
    Corrected

Corrected covariance estimator (scaled by `n-1` where `n` is the number of
elements in a vector).
"""
struct Corrected <: CovarianceEstimator
end

"""
    cov(::Corrected, x::AbstractVector)

Compute variance of the vector `x`. The sum is scaled with `n-1`
where `n = length(x)`.
"""
cov(::Corrected, x::AbstractVector) = cov(x; corrected = true)

"""
    cov(::Corrected, x::AbstractVector, y)

Compute covariance of the vectors `x` and `y` using formula
``\\frac{1}{n-1}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate.
"""
cov(::Corrected, x::AbstractVector, y::AbstractVector) = cov(x, y; corrected = true)

"""
    cov(::Corrected, X::AbstractMatrix)

Compute the covariance matrix of the matrix `X` along dimension `dims`.
The sum is scaled with `n-1` where `n = length(x)`.
"""
cov(::Corrected, X::AbstractMatrix; dims::Int=1) = cov(X; dims=dims, corrected = true)
