
using LinearAlgebra

struct RBLWCovariance <: CovarianceEstimator
end

function chenshrinkagetarget(X::DenseMatrix{<:Real})
    p = size(X, 2)
    C = cov(X; dims=2)
    (tr(C)/p) * one(C)
end

"""
    cov(RBLWCovariance(), X, shrinkage)

Calculates shrunk covariance matrix for centered data `X` with shrinkage
parameter `shrinkage` and Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::RBLWCovariance, X::DenseMatrix{T}, shrinkage::Number) where T<:Real
    (0 ≤ shrinkage ≤ 1) || throw(ArgumentError("Shinkage must be in [0,1] (given shrinkage: $shrinkage)"))
    C = cov(X; dims=2)
    F = chenshrinkagetarget(X)
    (1-shrinkage)*C + shrinkage*F
end

"""
    cov(::RBLWCovariance, X::DenseMatrix)

Calculates shrunk covariance matrix for centered data `X` with optimal shrinkage
for Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::RBLWCovariance, X::DenseMatrix{T}) where T<:Real
    C = cov(X; dims=2)
    n, p = size(X)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = ((n-2)/n * trS2 + tr2S)/((n+2) * (trS2 - tr2S/p))
    ρhat = min(ρhat, 1)
    F = chenshrinkagetarget(X)
    (1-ρhat)*C + ρhat*F
end

struct OASCovariance <: CovarianceEstimator
end

"""
    cov(OASCovariance(), X, shrinkage)

Calculates shrunk covariance matrix for centered data `X` with shrinkage
parameter `shrinkage` and Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::OASCovariance, X::DenseMatrix{T}, shrinkage::Number) where T<:Real
    cov(RBLWCovariance(), X, shrinkage)
end

"""
    cov(::OASCovariance, X::DenseMatrix)

Calculates shrunk covariance matrix for centered data `X` with optimal
oracle approximating shrinkage estimator.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::OASCovariance, X::DenseMatrix{T}) where T<:Real
    C = cov(X; dims=2)
    n, p = size(X)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = (-1.0/p * trS2 + tr2S)/((n-1.0)/p * (trS2 - tr2S/p))
    ρhat = min(ρhat, 1)
    F = chenshrinkagetarget(X)
    (1-ρhat)*C + ρhat*F
end
