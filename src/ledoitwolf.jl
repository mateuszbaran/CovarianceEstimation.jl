using LinearAlgebra

"""
    LedoitWolfCovariance(shrinkage)

Ledoit-Wolf covariance estimator. The parameter `shrinkage` is either equal
to `:optimal` and optimal shrinkage is calculated, or it is a number between
0 and 1.
"""
struct LedoitWolf{S<:Union{Symbol, Real}} <: CovarianceEstimator
    shrinkage::S
    function LedoitWolf(s::T) where T<:Real
        @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
        new{T}(s)
    end
    function LedoitWolf(s::Symbol=:optimal)
        @assert s ∈ [:optimal] "Shrinkage setting not supported"
        new{Symbol}(s)
    end
end

function lw_optimalshrinkage(X, C, F, r̄)
    # steps leading to equation 5 of http://www.ledoit.net/honey.pdf in
    # appendix B. (notations follow the paper)
    N, T  = size(X)
    tmat  = [(X[:,t] * X[:,t]' - C) for t ∈ 1:T]
    π̂mat  = sum(tmat[t].^2 for t ∈ 1:T) / T
    π̂     = sum(π̂mat)
    tdiag = [Diagonal(X[:,t].^2 - diag(C)) for t ∈ 1:T]
    ϑ̂ᵢᵢ   = sum(tdiag[t] * tmat[t] for t ∈ 1:T) / T # row scaling
    ϑ̂ⱼⱼ   = sum(tmat[t] * tdiag[t] for t ∈ 1:T) / T # col scaling
    ρ̂₂    = zero(eltype(X))
    # TODO: inbounds/simd?
    for i ∈ 1:N, j ∈ 1:N
        (j == i) && continue
        αᵢⱼ = sdC[j]/sdC[i]
        ρ̂₂ += ϑ̂ᵢᵢ[i,j]*αᵢⱼ + ϑ̂ⱼⱼ[i,j]/αᵢⱼ
    end
    ρ̂ = sum(diag(π̂mat)) + (r̄/2)*ρ̂₂
    γ̂ = sum((F - C).^2)
    # if γ̂ is very small it may lead to NaNs or infinities
    (γ̂ ≤ eps()) && return ifelse(π̂ ≤ ρ̂, 0.0, 1.0)
    κ̂ = (π̂ - ρ̂)/γ̂
    return clamp(κ̂/T, 0.0, 1.0)
end

function lw_shrinkagetarget(C)
    N = size(C, 1)
    Cs = [sqrt(C[i,i]*C[j,j]) for i in 1:N, j in 1:N]
    r = C ./ Cs
    r̄ = (sum(r)-N)/(N*(N-1))
    Finterm = Cs .* r̄
    F = Finterm + Diagonal(diag(Cs) .- diag(Finterm))
    return F, r̄
end

"""
    cov(::LedoitWolf, X::AbstractMatrix; dims::Int=1)

Calculates shrunk covariance matrix for data in `X` with Ledoit-Wolf
optimal shrinkage.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,”
The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
"""
function cov(lw::LedoitWolf, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 2
        X = transpose(X)
    elseif dims != 1
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    n, p = size(X)
    μ = mean(X, dims=1)
    for i ∈ 1:n, j ∈ 1:p
        X[i, j] -= μ[j]
    end
    C = cov(Simple(), X; dims=1)
    sdC = sqrt.(diag(C))
    F, r̄ = lw_shrinkagetarget(C, sdC)
    shrinkage = lw.shrinkage
    (shrinkage == :optimal) && (shrinkage = lw_optimalshrinkage(X, C, sdC, F, r̄))
    return (1.0 - shrinkage) * C + shrinkage * F
end
