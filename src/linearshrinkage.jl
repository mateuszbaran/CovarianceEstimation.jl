const Shrinkage = Union{Symbol, Real}
abstract type LinearShrinkageTarget end

"""
    linear_shrinkage(S, F, ρ)

Given a sample covariance `S` and a matching target `F`, return an estimator
corresponding to the James-Stein estimator corresponding to a shrinkage of
intensity `ρ` between `S` and `F`.
"""
linear_shrinkage(S::AbstractMatrix, F::LinearShrinkageTarget, ρ::Real) =
    (1.0 - ρ) * S + ρ * F


struct LinearShrinkageEstimator{T<:LinearShrinkageTarget, S<:Shrinkage} <: CovarianceEstimator
    target::T
    shrinkage::S
end

function LinearShrinkageEstimator(t::LinearShrinkageTarget, s::Real)
    @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
    return LinearShrinkageEstimator(t, s)
end

function LinearShrinkageEstimator(t::LinearShrinkageTarget, s::Symbol=:optimal)
    @assert s ∈ [:optimal] "Shrinkage setting not supported"
    return LinearShrinkageEstimator(t, s)
end


function cov(X::AbstractMatrix{<:Real}, lse::LinearShrinkageEstimator;
             dims::Int=1)

    Xc = copy(X)
    if dims == 2
        Xc = transpose(Xc)
    elseif dims != 1
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    centercols!(Xc)
    # sample covariance of size (p x p)
    n, p = size(Xc)
    S    = (Xc'*Xc)/n

    F, ρ = targetandshrinkage(lse.target, S, Xc)
    ρ    = ifelse(lse.shrinkage == :optimal, ρ : ls.shrinkage)
    return shrink(S, F, ρ̂)
end
