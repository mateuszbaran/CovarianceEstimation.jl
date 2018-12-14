"""
    Simple(corrected::Bool)

Uncorrected covariance estimator (scaled by `n` where `n` is the number of
elements in a vector).
"""
struct Simple <: CovarianceEstimator
    corrected::Bool
    function Simple(;corrected::Bool = false)
        new(corrected)
    end
end

"""
    cov(sc::Simple, x::AbstractVector)

Compute variance of the vector `x`. The sum is scaled with `n`
where `n = length(x)`.
"""
cov(x::AbstractVector, sc::Simple) = cov(x; corrected = sc.corrected)

"""
    cov(sc::Simple, x::AbstractVector, y)

Compute covariance of the vectors `x` and `y` using formula
``\\frac{1}{n}\\sum_{i=1}^n (x_i-\\bar x) (y_i-\\bar y)^*`` where ``*`` denotes
complex conjugate.
"""
cov(x::AbstractVector, y::AbstractVector, sc::Simple) = cov(x, y; corrected = sc.corrected)

"""
    cov(sc::Simple, X::AbstractMatrix)

Compute the covariance matrix of the matrix `X` along dimension `dims`.
The sum is scaled with `n` where `n = length(x)`.
"""
cov(X::AbstractMatrix, sc::Simple; dims::Int=1) = cov(X; dims=dims, corrected = sc.corrected)
