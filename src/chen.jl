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

function rblw_optimalshrinkage(Ŝ::AbstractMatrix, n::Int, p::Int)
    # https://arxiv.org/pdf/0907.4698.pdf equations 17, 19
    trŜ² = dot(Ŝ, transpose(Ŝ))
    tr²Ŝ = tr(Ŝ)^2
    ρ̂    = ((n-2)/n * trŜ² + tr²Ŝ)/((n+2) * (trŜ² - tr²Ŝ/p))
    return min(ρ̂, 1.0)
end

rblw_shrinkagetarget(Ŝ::AbstractMatrix, p::Int) = (tr(Ŝ)/p) * I

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
    Xc = copy(X)
    if dims == 2
        Xc = transpose(Xc)
    elseif dims != 1
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    centercols!(Xc)
    # sample covariance of size (p x p)
    n, p = size(Xc)
    Ŝ    = (Xc'*Xc)/n
    # shrinkage
    F = rblw_shrinkagetarget(Ŝ, p)
    ρ = rblw.shrinkage
    (ρ == :optimal) && (ρ = rblw_optimalshrinkage(Ŝ, n, p))
    return shrink(Ŝ, F, ρ)
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

function oas_optimalshrinkage(Ŝ::AbstractMatrix, n::Int, p::Int)
    # https://arxiv.org/pdf/0907.4698.pdf equation 23
    trŜ² = dot(Ŝ, transpose(Ŝ))
    tr²Ŝ = tr(Ŝ)^2
    ρ̂ = ((1.0-2.0/p) * trŜ² + tr²Ŝ)/((n+1.0-2.0/p) * (trŜ² - tr²Ŝ/p))
    return min(ρ̂, 1)
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
    Xc = copy(X)
    if dims == 2
        Xc = transpose(Xc)
    elseif dims != 1
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    centercols!(Xc)
    # sample covariance of size (p x p)
    n, p = size(Xc)
    Ŝ    = (Xc'*Xc)/n
    # shrinkage
    F = rblw_shrinkagetarget(Ŝ, p)
    ρ = oas.shrinkage
    (ρ == :optimal) && (ρ = oas_optimalshrinkage(Ŝ, n, p))
    return shrink(Ŝ, F, ρ)
end
