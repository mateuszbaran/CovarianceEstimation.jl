
"""
    SimpleCovariance

Uncorrected covariance estimator (scaled by `n` where `n` is the number of
elements in a vector).
"""
struct SimpleCovariance <: CovarianceEstimator
end

"""
    cov(SimpleCovariance(), x::AbstractVector)

Compute variance of the vector `x`. The sum is scaled with `n`
where `n = length(x)`.
"""
cov(::SimpleCovariance, x::AbstractVector) = cov(x; corrected = false)

"""
    cov(SimpleCovariance(), x::AbstractVector, y)

Compute covariance of the vectors `x` and `y` using formula
``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate.
"""
cov(::SimpleCovariance, x::AbstractVector, y::AbstractVector) = cov(x, y; corrected = false)

"""
    cov(SimpleCovariance(), X::AbstractMatrix)

Compute the covariance matrix of the matrix `X` along dimension `dims`.
The sum is scaled with `n` where `n = length(x)`.
"""
cov(::SimpleCovariance, X::AbstractMatrix; dims::Int=1) = cov(X; dims=dims, corrected = false)


"""
    CorrectedCovariance

Corrected covariance estimator (scaled by `n-1` where `n` is the number of
elements in a vector).
"""
struct CorrectedCovariance <: CovarianceEstimator
end

"""
    cov(CorrectedCovariance(), x::AbstractVector)

Compute variance of the vector `x`. The sum is scaled with `n-1`
where `n = length(x)`.
"""
cov(::CorrectedCovariance, x::AbstractVector) = cov(x; corrected = true)

"""
    cov(CorrectedCovariance(), x::AbstractVector, y)

Compute covariance of the vectors `x` and `y` using formula
``\\frac{1}{n-1}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate.
"""
cov(::CorrectedCovariance, x::AbstractVector, y::AbstractVector) = cov(x, y; corrected = true)

"""
    cov(CorrectedCovariance(), X::AbstractMatrix)

Compute the covariance matrix of the matrix `X` along dimension `dims`.
The sum is scaled with `n-1` where `n = length(x)`.
"""
cov(::CorrectedCovariance, X::AbstractMatrix; dims::Int=1) = cov(X; dims=dims, corrected = true)
