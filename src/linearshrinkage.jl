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


function cov(X::AbstractMatrix{<:Real}, lse::LinearShrinkageEstimator;
             corrected::Bool=false, dims::Int=1)

    @assert dims ∈ [1, 2] "Argument dims can only be 1 or 2 (given: $dims)"

    Xc = (dims == 1) ? centercols(X) : centercols(transpose(X))
    # sample covariance of size (p x p)
    n, p = size(Xc)
    S    = cov(Xc, Simple(corrected=corrected))
    return linear_shrinkage(lse.target, Xc, S, lse.shrinkage, n, p, corrected)
end

##################################################

# helper function
function sum_var_sij(Z::AbstractMatrix, CovZ::AbstractMatrix, n::Int,
                     corrected=false; with_diag=true)
    scale = ifelse(corrected, n-1, n)
    Z²    = Z.^2
    π̂mat  = ((Z²'*Z²)/n - CovZ.^2) / scale
    with_diag && return sum(π̂mat)
    return sum(π̂mat) - sum(diag(π̂mat))
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

# helper function for schafer-strimmer, compute the standardised data matrix
# and the corresponding uncorrected covariance.
function std_cov(Xc::AbstractMatrix, dS::AbstractVector)
    Xs = Xc * Diagonal(1.0 ./ sqrt.(dS))
    R  = cov(Xs, Simple())
    return Xs, R
end

## TARGET A

function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, I, λ)
end

function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F = I
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ = sum_var_sij(Xc, S, n, corrected) / sum((S - F).^2)
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

    F = target_B(S, p)
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ = sum_var_sij(Xc, S, n, corrected) / sum((S - F).^2)
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

target_D(S::AbstractMatrix) = Diagonal(diag(S))

function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, target_D(S), λ)
end

function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F = target_D(S)
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ = sum_var_sij(Xc, S, n, corrected; with_diag=false) / sum((S - F).^2)
    elseif λ == :ss
        # use the standardised data matrix and its covariance
        Xs, R = std_cov(Xc, diag(S))
        λ  = sum_var_sij(Xs, R, n, corrected; with_diag=false)
        λ /= sum((R - target_D(R)).^2)
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end

## TARGET C

function target_C(S::AbstractMatrix, p::Int)
    v = tr(S)/p
    # average of off-diagonal terms
    c = sum(S)/(p * (p - 1)) - v / (p - 1)
    # target: off diag terms = average of s_{ij} for i≂̸j
    F = c * ones(p, p)
    # target: diag terms = average of s_{ii}
    F -= Diagonal(diag(F))
    F += v * I
    return F
end

function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int,
                          corrected::Bool)

    return linshrink(S, target_C(S, p), λ)
end

function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int,
                          corrected::Bool)

    F = target_C(S, p)
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ = sum_var_sij(Xc, S, n, corrected) / sum((S - F).^2)
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

    F = target_E(S)
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ  = sum_var_sij(Xc, S, n, corrected; with_diag=false)
        λ -= sum_fij(Xc, S, n, p, corrected)
        λ /= sum((S - F).^2)
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
    F  = F_ + Diagonal(diag(s_) .- diag(F_))
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
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        λ  = sum_var_sij(Xc, S, n, corrected; with_diag=false)
        λ -= r̄ * sum_fij(Xc, S, n, p, corrected)
        λ /= sum((S - F).^2)
    else
        error("Unsupported shrinkage method for target DiagonalCommonVariance.")
    end
    λ = clamp(λ, 0.0, 1.0)
    return linshrink(S, F, λ)
end
