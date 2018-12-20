const Shrinkage = Union{Symbol, Real}
abstract type LinearShrinkageTarget end

# Taxonomy from http://strimmerlab.org/publications/journals/shrinkcov2005.pdf
# page 13

struct DiagonalUnitVariance       <: LinearShrinkageTarget end
struct DiagonalCommonVariance     <: LinearShrinkageTarget end
struct DiagonalUnequalVariance    <: LinearShrinkageTarget end
struct CommonCovariance           <: LinearShrinkageTarget end
struct PerfectPositiveCorrelation <: LinearShrinkageTarget end
struct ConstantCorrelation        <: LinearShrinkageTarget end


struct LinearShrinkageEstimator{T<:LinearShrinkageTarget, S<:Shrinkage} <: CovarianceEstimator
    target::T
    shrinkage::S
    function LinearShrinkageEstimator(t::TT, s::SS) where {TT<:LinearShrinkageTarget, SS<:Real}
        @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
        new{TT, SS}(t, s)
    end
    function LinearShrinkageEstimator(t::TT, s::Symbol=:auto) where TT <: LinearShrinkageTarget
        @assert s ∈ [:auto, :lw, :ss, :rblw, :oas] "Shrinkage method not supported"
        new{TT, Symbol}(t, s)
    end
end

LinearShrinkageEstimator(;
    target::LinearShrinkageTarget=DiagonalUnitVariance(),
    shrinkage::Shrinkage) = LinearShrinkageEstimator(target, shrinkage)


function cov(X::AbstractMatrix{<:Real}, lse::LinearShrinkageEstimator;
             corrected::Bool=false, dims::Int=1)

    @assert dims ∈ [1, 2] "Argument dims can only be 1 or 2 (given: $dims)"

    Xc = (dims == 1) ? centercols(X) : centercols(transpose(X))
    # sample covariance of size (p x p)
    n, p = size(Xc)
    S    = cov(Xc, Simple(corrected=corrected))
    return linear_shrinkage(lse.target, Xc, S, lse.shrinkage, n, p, corrected)
end

##############################################################################
# Helper functions for Linear Shrinkage estimators

# M: square (symmetric)
rescale(M::AbstractMatrix, D::Diagonal) = D * M * D

# uncentered cov
uccov(X::AbstractMatrix) = (X'*X)/size(X, 1)

# sum of off diagonal elements, basically just a sum across elements
# of a square matrix with the possibility of removing the diagonal.
function sumij(S::AbstractMatrix; with_diag=false)
    acc = sum(S)
    with_diag || return acc - tr(S)
    return acc
end

# helper function ∑_{i≂̸j} f_ij
function sum_fij(Xc, S, n, p, corrected=false)
    tmat  = @inbounds [(Xc[t,:] * Xc[t,:]' - S) for t ∈ 1:n]
    dS    = diag(S)
    sdS   = sqrt.(dS)
    tdiag = @inbounds [Diagonal(Xc[t,:].^2 - dS) for t ∈ 1:n]
    # estimator for cov(s_ii, s_ij) --> row scaling
    ϑ̂ᵢᵢ   = @inbounds sum(tdiag[t] * tmat[t] for t ∈ 1:n) / n
    # estimator for cov(s_jj, s_ij) --> column scaling
    ϑ̂ⱼⱼ   = @inbounds sum(tmat[t] * tdiag[t] for t ∈ 1:n) / n
    ∑fij  = zero(eltype(Xc))
    @inbounds for i ∈ 1:p, j ∈ 1:p
        (j == i) && continue
        αᵢⱼ   = sdS[j]/sdS[i]
        ∑fij += ϑ̂ᵢᵢ[i,j]*αᵢⱼ + ϑ̂ⱼⱼ[i,j]/αᵢⱼ
    end
    ∑fij /= (2.0 * ifelse(corrected, n-1, n))
    return ∑fij
end
##############################################################################

## TARGET A

function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, I, λ)
end

function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F   = I
    κ   = ifelse(corrected, n-1, n)
    γ   = κ/n
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS̄² = γ^2 * sumij(S.^2, with_diag=true)
        λ   = (sumij(uccov(Xc²), with_diag=true) - ΣS̄²) / κ
        λ  /=  (ΣS̄² - 2tr(S) + p)
    else
        error("Unsupported shrinkage method for target DiagonalUnitVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET B

target_B(S::AbstractMatrix, p::Int) = tr(S)/p * I

function linear_shrinkage(::DiagonalCommonVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, target_B(S, p), λ)
end

function linear_shrinkage(::DiagonalCommonVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F   = target_B(S, p)
    κ   = ifelse(corrected, n-1, n)
    γ   = κ/n
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        v   = F.λ # tr(S)/p
        ΣS̄² = γ^2 * sumij(S.^2, with_diag=true)
        λ   = (sumij(uccov(Xc²), with_diag=true) - ΣS̄²) / κ
        λ  /=  (ΣS̄² - p*v^2)
    elseif λ == :rblw
        # https://arxiv.org/pdf/0907.4698.pdf equations 17, 19
        trS² = sum(S.^2)
        tr²S = tr(S)^2
        # note: using corrected or uncorrected S does not change λ
        λ = ((n-2)/n * trS² + tr²S)/((n+2) * (trS² - tr²S/p))
    elseif λ == :oas
        # https://arxiv.org/pdf/0907.4698.pdf equation 23
        trS² = sum(S.^2)
        tr²S = tr(S)^2
        # note: using corrected or uncorrected S does not change λ
        λ = ((1.0-2.0/p) * trS² + tr²S)/((n+1.0-2.0/p) * (trS² - tr²S/p))
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET D

target_D(S::AbstractMatrix) = Diagonal(S)

function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, target_D(S), λ)
end

function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F   = target_D(S)
    κ   = ifelse(corrected, n-1, n)
    γ   = κ/n
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS̄² = γ^2 * sumij(S.^2)
        λ   = (sumij(uccov(Xc²)) - ΣS̄²) / (κ * ΣS̄²)
    elseif λ == :ss
        # use the standardised data matrix
        D   = Diagonal(1.0 ./ sum(Xc², dims=1)[:])
        ΣS̄² = γ^2 * sumij(rescale(S.^2, D))
        λ   = (sumij(rescale(uccov(Xc²), D)) - ΣS̄²) / κ
        λ  /= ΣS̄²
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET C

function target_C(S::AbstractMatrix, p::Int)
    v = tr(S)/p
    c = sumij(S; with_diag=false) / (p*(p - 1))
    F = c * ones(p, p)
    F -= Diagonal(F)
    F += v * I
    return F, v, c
end

function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    F, _, _ = target_C(S, p)
    return linshrink(S, F, λ)
end

function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F, v, c = target_C(S, p)
    κ   = ifelse(corrected, n-1, n)
    γ   = κ/n
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS̄² = γ^2 * sumij(S.^2, with_diag=true)
        λ   = (sumij(uccov(Xc²), with_diag=true) - ΣS̄²) / κ
        λ  /=  (ΣS̄² - p*(p-1)*c^2 - p*v^2)
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET E

function target_E(S::AbstractMatrix)
    d = diag(S)
    return sqrt.(d*d')
end

function linear_shrinkage(::PerfectPositiveCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, target_E(S), λ)
end

function linear_shrinkage(::PerfectPositiveCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F   = target_E(S)
    κ   = ifelse(corrected, n-1, n)
    γ   = κ/n
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS̄² = γ^2 * sumij(S.^2)
        λ   = (sumij(uccov(Xc²)) - ΣS̄²) / κ
        λ  -= sum_fij(Xc, S, n, p, corrected)
        λ  /= sumij((S - F).^2)
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET F

function target_F(S::AbstractMatrix, p::Int)
    s  = sqrt.(diag(S))
    s_ = @inbounds [s[i]*s[j] for i ∈ 1:p, j ∈ 1:p]
    r̄  = (sum(S ./ s_) - p)/(p * (p - 1))
    F_ = r̄ * s_
    F  = F_ + (Diagonal(s_) - Diagonal(F_))
    return F, r̄
end

function linear_shrinkage(::ConstantCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    F, _ = target_F(S, p)
    return linshrink(S, F, λ)
end

function linear_shrinkage(::ConstantCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F, r̄ = target_F(S, p)
    κ    = ifelse(corrected, n-1, n)
    γ    = κ/n
    Xc²  = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS̄² = γ^2 * sumij(S.^2)
        λ   = (sumij(uccov(Xc²)) - ΣS̄²) / κ
        λ  -= r̄ * sum_fij(Xc, S, n, p, corrected)
        λ  /= sumij((S - F).^2)
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end
