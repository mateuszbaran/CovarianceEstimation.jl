
using LinearAlgebra

struct RBLWCovariance <: CovarianceEstimator
end

function chenshrinkagetarget(X::AbstractMatrix{<:Real})
    p = size(X, 2)
    C = cov(X; dims=2)
    (tr(C)/p) * one(C)
end

"""
    cov(RBLWCovariance(), X::AbstractMatrix, shrinkage; dims::Int=1)

Calculates shrunk covariance matrix for centered data `X` with shrinkage
parameter `shrinkage` and Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::RBLWCovariance, X::AbstractMatrix{T}, shrinkage::Number; dims::Int=1) where T<:Real
    (0 ≤ shrinkage ≤ 1) || throw(ArgumentError("Shinkage must be in [0,1] (given shrinkage: $shrinkage)"))
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    C = cov(Xint; dims=2)
    F = chenshrinkagetarget(Xint)
    (1-shrinkage)*C + shrinkage*F
end

"""
    cov(::RBLWCovariance, X::AbstractMatrix; dims::Int=1)

Calculates shrunk covariance matrix for centered data `X` with optimal shrinkage
for Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::RBLWCovariance, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end

    C = cov(Xint; dims=2)
    n, p = size(Xint)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = ((n-2)/n * trS2 + tr2S)/((n+2) * (trS2 - tr2S/p))
    ρhat = min(ρhat, 1)
    F = chenshrinkagetarget(Xint)
    (1-ρhat)*C + ρhat*F
end

struct OASCovariance <: CovarianceEstimator
end

"""
    cov(OASCovariance(), X::AbstractMatrix, shrinkage; dims::Int=1)

Calculates shrunk covariance matrix for centered data `X` with shrinkage
parameter `shrinkage` and Rao-Blackwell theorem Ledoit-Wolf shrinkage target.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::OASCovariance, X::AbstractMatrix{T}, shrinkage::Number; dims::Int=1) where T<:Real
    cov(RBLWCovariance(), X, shrinkage; dims = dims)
end

"""
    cov(::OASCovariance, X::AbstractMatrix; dims::Int=1)

Calculates shrunk covariance matrix for centered data `X` with optimal
oracle approximating shrinkage estimator.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function cov(::OASCovariance, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end

    C = cov(Xint; dims=2)
    n, p = size(Xint)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = (-1.0/p * trS2 + tr2S)/((n-1.0)/p * (trS2 - tr2S/p))
    ρhat = min(ρhat, 1)
    F = chenshrinkagetarget(Xint)
    (1-ρhat)*C + ρhat*F
end
