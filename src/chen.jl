"""
    RaoBlackwellLedoitWolf()

Rao-Blackwell theorem Ledoit-Wolf covariance estimator.
"""
struct RaoBlackwellLedoitWolf <: TargetType end

function rblw_optimalshrinkage(S::AbstractMatrix, n::Int, p::Int)
    # https://arxiv.org/pdf/0907.4698.pdf equations 17, 19
    trS² = dot(S, transpose(S))
    tr²S = tr(S)^2
    ρ    = ((n-2)/n * trS² + tr²S)/((n+2) * (trS² - tr²S/p))
    return min(ρ, one(ρ))
end

rblw_shrinkagetarget(S::AbstractMatrix, p::Int) = (tr(S)/p) * I

"""
    targetandshrinkage(rblw::RaoBlackwellLedoitWolf, S::AbstractMatrix X::AbstractMatrix)

Calculates shrunk covariance matrix for centered data `X` and simple
covariance `S` with optimal shrinkage for Rao-Blackwell theorem Ledoit-Wolf
shrinkage target.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function targetandshrinkage(::RaoBlackwellLedoitWolf, S::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real})
    # shrinkage
    n, p = size(X)
    F = rblw_shrinkagetarget(S, p)
    ρ = rblw_optimalshrinkage(S, n, p)
    return F, ρ
end

"""
    OracleApproximatingShrinkage(shrinkage)

Oracle approximating shrinkage covariance estimator.
"""
struct OracleApproximatingShrinkage <: TargetType end

function oas_optimalshrinkage(S::AbstractMatrix, n::Int, p::Int)
    # https://arxiv.org/pdf/0907.4698.pdf equation 23
    trS² = dot(S, transpose(S))
    tr²S = tr(S)^2
    ρ    = ((1.0-2.0/p) * trS² + tr²S)/((n+1.0-2.0/p) * (trS² - tr²S/p))
    return min(ρ, one(ρ))
end

"""
    targetandshrinkage(::OracleApproximatingShrinkage, S::AbstractMatrix, X::AbstractMatrix)

Calculates shrunk covariance matrix for centered data `X` and simple
covariance `S` with optimal oracle approximating shrinkage estimator.

Implements shrinkage target and optimal shrinkage according to
Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero,
“Shrinkage Algorithms for MMSE Covariance Estimation,”
IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
"""
function targetandshrinkage(::OracleApproximatingShrinkage, S::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real})
    n, p = size(X)
    F = rblw_shrinkagetarget(S, p)
    ρ = oas_optimalshrinkage(S, n, p)
    return F, ρ
end
