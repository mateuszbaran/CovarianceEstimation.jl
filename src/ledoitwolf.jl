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

function lw_optimalshrinkage(Xc::AbstractMatrix, Ŝ::AbstractMatrix,
                             V̂::AbstractVector, F::AbstractMatrix,
                             r̄::Real, n::Int, p::Int)
    # steps leading to equation 5 of http://www.ledoit.net/honey.pdf in
    # appendix B. (notations follow the paper)
    tmat  = @inbounds [(Xc[t,:] * Xc[t,:]' - Ŝ) for t ∈ 1:n]
    π̂mat  = sum(tmat[t].^2 for t ∈ 1:n) / n
    π̂     = sum(π̂mat)
    dŜ    = diag(Ŝ)
    tdiag = @inbounds [Diagonal(Xc[t,:].^2 - dŜ) for t ∈ 1:n]
    ϑ̂ᵢᵢ   = @inbounds sum(tdiag[t] * tmat[t] for t ∈ 1:n) / n # row scaling
    ϑ̂ⱼⱼ   = @inbounds sum(tmat[t] * tdiag[t] for t ∈ 1:n) / n # col scaling
    ρ̂₂    = zero(eltype(Xc))
    @inbounds for i ∈ 1:p, j ∈ 1:p
        (j == i) && continue
        αᵢⱼ = V̂[j]/V̂[i]
        ρ̂₂ += ϑ̂ᵢᵢ[i,j]*αᵢⱼ + ϑ̂ⱼⱼ[i,j]/αᵢⱼ
    end
    ρ̂ = sum(diag(π̂mat)) + (r̄/2)*ρ̂₂
    γ̂ = sum((F - Ŝ).^2)
    # if γ̂ is very small it may lead to NaNs or infinities
    (γ̂ ≤ eps()) && return ifelse(π̂ ≤ ρ̂, 0.0, 1.0)
    κ̂ = (π̂ - ρ̂)/γ̂
    return clamp(κ̂/n, 0.0, 1.0)
end

function lw_shrinkagetarget(Ŝ::AbstractMatrix, V̂::AbstractVector, p::Int)
    V̂_ = @inbounds [V̂[i]*V̂[j] for i ∈ 1:p, j ∈ 1:p]
    r  = Ŝ ./ V̂_
    r̄  = (sum(r) - p)/(p * (p - 1))
    F_ = V̂_ .* r̄
    F  = F_ + Diagonal(diag(V̂_) .- diag(F_))
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
    V̂    = sqrt.(diag(Ŝ))
    # shrinkage
    F, r̄ = lw_shrinkagetarget(Ŝ, V̂, p)
    ρ    = lw.shrinkage
    (ρ == :optimal) && (ρ = lw_optimalshrinkage(Xc, Ŝ, V̂, F, r̄, n, p))
    return shrink(Ŝ, F, ρ)
end
