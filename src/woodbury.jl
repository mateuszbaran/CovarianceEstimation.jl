# Covariance estimation for high-dimensional data
# Uses a covariance model of the form `Σ = σ²I + U * Λ * U'`, where `U` is a low-rank matrix of eigenvectors and
# `Λ` is a diagonal matrix capturing the excess width along the dimensions in `U` compared to isotropic.

# If you're curious about dispatch and inferrability, see the "Implementation note" at the top of src/loss.jl.
# We employ the same de-specialization trick even for `σ²` in `WoodburyEstimator`, as we'll adopt the eltype of
# `Λ` anyway.

"""
    WoodburyEstimator(loss::LossFunction, rank::Integer;
                      σ²::Union{Real,Nothing}=nothing, corrected::Bool=false)

Specify that covariance matrices should be estimated using a "spiked" covariance model

    Σ = σ²I + U * Λ * U'

`loss` is either a [`NormLossCov`](@ref) or [`StatLossCov`](@ref) object, which specifies the
loss function for which the estimated covariance will be optimal. `rank` is the maximum
number of eigenvalues `Λ` to retain in the model. Optionally, one may specify `σ²` directly,
or it can be estimated from the data matrix (`σ²=nothing`). Set `corrected=true` to use
the unbiased estimator of the variance.
"""
struct WoodburyEstimator{L<:LossFunction} <: CovarianceEstimator
    loss::L
    rank::Int
    σ²::Union{Real,Nothing} # common diagonal variance, `nothing` indicates unknown
    corrected::Bool
end
WoodburyEstimator(loss::LossFunction, rank::Integer; σ²::Union{Real,Nothing}=nothing, corrected::Bool=false) =
    WoodburyEstimator(loss, rank, σ², corrected)

"""
    cov(estimator::WoodburyEstimator, X::AbstractMatrix, weights::FrequencyWeights...; dims::Int=1, mean=nothing, UsV=nothing)

Estimate the covariance matrix from the data matrix `X` using a "spiked" covariance model

    Σ = σ²I + U * Λ * U',

where `U` is a low-rank matrix of eigenvectors and `Λ` is a diagonal matrix.

Reference:
    Donoho, D.L., Gavish, M. and Johnstone, I.M., 2018.
    Optimal shrinkage of eigenvalues in the spiked covariance model. Annals of statistics, 46(4), p.1742.

When `σ²` is not supplied in `estimator`, it is calculated from the residuals `X - X̂`, where `X̂` is the
low-rank approximation of `X` used to generate `U` and `Λ`.

If `X` is too large to manipulate in memory, you can pass `UsV = (U, s, V)` (a truncated SVD of `X - mean(X; dims)`)
and then `X` will only be used compute the dimensionality and number of observations. This requires that you
specify `estimator.σ²`.
"""
function cov(estimator::WoodburyEstimator,  X::AbstractMatrix{<:Real}, weights::FrequencyWeights...;
             dims::Int=1, mean=nothing, UsV = nothing)
    # Argument validation
    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    p = size(X, 3 - dims)
    p >= estimator.rank || throw(ArgumentError("Argument rank (got $(estimator.rank)) must be less than the number of observations (size(X, dims)=$(size(X, dims)))"))
    wn = totalweight(size(X, dims), weights...)

    local ΔX
    U, s, V = if UsV === nothing
        if mean === nothing
            mean = Statistics.mean(X, weights...; dims=dims)
        end
        # Compute the low-rank approximation of the centered data matrix
        ΔX = weightedX(X .- mean, weights...; dims=dims)
        tsvd(ΔX, estimator.rank)
    else
        UsV
    end

    T = eltype(s)
    σ² = estimator.σ²
    σ² = if σ² === nothing
        ΔΔX = ΔX - U*Diagonal(s)*V'
        # The number of degrees of freedom is (number of observations minus the rank)*dimensionality
        nσ = (totalweight(size(X, dims), weights...) - estimator.rank) * size(X, 3-dims)
        sum(abs2, ΔΔX) / (nσ - estimator.corrected)
    else
        T(σ²)
    end::T     # fix inferrability (see note at top of file)

    # Ratio of dimensionality to number of observations (the principal parameter in Random Matrix Theory)
    γ = p / wn

    # Implement the optimal shrinkage algorithm
    λ_shrunk = shrink.(Ref(estimator.loss), s.^2 ./ wn, σ², γ)
    keep = (!iszero).(λ_shrunk)

    # Return the shrunk covariance matrix as a WoodburyMatrix
    return SymWoodbury(σ² * I(p), dims == 1 ? V[:, keep] : U[:, keep], Diagonal(λ_shrunk[keep]))
end

function shrink(loss::LossFunction, λ::Real, σ²::Real, γ::Real)
    # Implement the procedure on Donoho et al. (2018), p. 1758
    # We return the difference from σ², since that's already contained in the diagonal term
    λu = λ / σ²
    λ₊ = (1 + sqrt(γ))^2
    λu < λ₊ && return zero(σ²)
    # Calculate the "de-biased" eigenvalue ℓ (Eq. 1.10)
    λ′ = λu + 1 - γ
    ℓ = (λ′ + sqrt(λ′^2 - 4λu)) / 2
    # Calculate the cosine (Eq. 1.6)
    c = sqrt((1 - γ / (ℓ - 1)^2) / (1 + γ / (ℓ - 1)^2))
    # Calculate the sine
    s = sqrt(1 - c^2)
    # Apply the shrinker
    return σ² * (shrinker(loss, ℓ, c, s) - 1)
end
