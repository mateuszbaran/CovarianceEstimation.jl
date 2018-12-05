
using LinearAlgebra

"""
    RaoBlackwellLedoitWolf(shrinkage)

Rao-Blackwell theorem Ledoit-Wolf covariance estimator. The parameter
`shrinkage` is either equal to `:auto` and optimal shrinkage is calculated,
or it is a number between 0 and 1.
"""
struct RaoBlackwellLedoitWolf{S<:Union{Symbol, Real}} <: CovarianceEstimator
    shrinkage::S
end

RaoBlackwellLedoitWolf() = RaoBlackwellLedoitWolf{Symbol}(:auto)

function chenshrinkagetarget(X::AbstractMatrix{<:Real}; dims=2)
    p = size(X, dims)
    C = cov(X; dims=dims)
    (tr(C)/p) * one(C)
end

"""
    cov(rblw::RaoBlackwellLedoitWolf, X::AbstractMatrix; dims::Int=1)

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
function cov(rblw::RaoBlackwellLedoitWolf, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end

    C = cov(Xint; dims=2)
    shrinkage =
    if rblw.shrinkage isa Symbol
        n, p = size(Xint)
        trS2 = dot(C, transpose(C)) # trace of C*C
        tr2S = tr(C)^2
        ρhat = ((n-2)/n * trS2 + tr2S)/((n+2) * (trS2 - tr2S/p))
        min(ρhat, 1) # assigned to variable `shrinkage`
    else
        rblw.shrinkage # assigned to variable `shrinkage`
    end

    F = chenshrinkagetarget(Xint)
    (1-shrinkage)*C + shrinkage*F
end

"""
    OracleApproximatingShrinkage(shrinkage)

Oracle approximating shrinkage covariance estimator. The parameter
`shrinkage` is either equal to `:auto` and optimal shrinkage is calculated,
or it is a number between 0 and 1.
"""
struct OracleApproximatingShrinkage{S<:Union{Symbol, Real}} <: CovarianceEstimator
    shrinkage::S
end

OracleApproximatingShrinkage() = OracleApproximatingShrinkage{Symbol}(:auto)

"""
    cov(oas::OracleApproximatingShrinkage, X::AbstractMatrix; dims::Int=1)

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
function cov(oas::OracleApproximatingShrinkage, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end

    C = cov(Xint; dims=2)
    shrinkage =
    if oas.shrinkage isa Symbol
        p, n = size(Xint)
        trS2 = dot(C, transpose(C)) # trace of C*C
        tr2S = tr(C)^2
        ρhat = ((1-2.0/p) * trS2 + tr2S)/((n+1-2.0/p) * (trS2 - tr2S/p))
        min(ρhat, 1) # assigned to variable `shrinkage`
    else
        oas.shrinkage # assigned to variable `shrinkage`
    end
    F = chenshrinkagetarget(Xint)
    (1-shrinkage)*C + shrinkage*F
end
