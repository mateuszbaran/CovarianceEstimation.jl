
using LinearAlgebra

"""
    RaoBlackwellLedoitWolf(shrinkage)

Rao-Blackwell theorem Ledoit-Wolf covariance estimator. The parameter
`shrinkage` is either equal to `:optimal` and optimal shrinkage is calculated,
or it is a number between 0 and 1.
"""
struct RaoBlackwellLedoitWolf{S<:Union{Symbol, Real}} <: CovarianceEstimator
    shrinkage::S
    function RaoBlackwellLedoitWolf(s::T) where T<:Real
        @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
        new{T}(s)
    end
    function RaoBlackwellLedoitWolf(s::Symbol=:optimal)
        @assert s ∈ [:optimal] "Shrinkage setting not supported"
        new{Symbol}(s)
    end
end

function rblw_optimalshrinkage(n, p, C)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = ((n-2)/n * trS2 + tr2S)/((n+2) * (trS2 - tr2S/p))
    return min(ρhat, 1)
end

function rblw_shrinkagetarget(X; dims=2)
    p = size(X, dims)
    C = cov(X; dims=dims)
    return (tr(C)/p) * one(C)
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
    n, p = size(Xint)
    shrinkage = rblw.shrinkage
    (shrinkage == :optimal) && (shrinkage = rblw_optimalshrinkage(n, p, C))
    F = rblw_shrinkagetarget(Xint)
    return (1.0 - shrinkage) * C + shrinkage * F
end

"""
    OracleApproximatingShrinkage(shrinkage)

Oracle approximating shrinkage covariance estimator. The parameter
`shrinkage` is either equal to `:optimal` and optimal shrinkage is calculated,
or it is a number between 0 and 1.
"""
struct OracleApproximatingShrinkage{S<:Union{Symbol, Real}} <: CovarianceEstimator
    shrinkage::S
    function OracleApproximatingShrinkage(s::T) where T<:Real
        @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
        new{T}(s)
    end
    function OracleApproximatingShrinkage(s::Symbol=:optimal)
        @assert s ∈ [:optimal] "Shrinkage setting not supported"
        new{Symbol}(s)
    end
end

function oas_optimalshrinkage(n, p, C)
    trS2 = dot(C, transpose(C)) # trace of C*C
    tr2S = tr(C)^2
    ρhat = ((1.0-2.0/p) * trS2 + tr2S)/((n+1.0-2.0/p) * (trS2 - tr2S/p))
    return min(ρhat, 1) # assigned to variable `shrinkage`
end

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
    n, p = size(Xint)
    shrinkage = oas.shrinkage
    (shrinkage == :optimal) && (shrinkage = oas_optimalshrinkage(n, p, C))
    F = rblw_shrinkagetarget(Xint)
    return (1.0 - shrinkage) * C + shrinkage * F
end
