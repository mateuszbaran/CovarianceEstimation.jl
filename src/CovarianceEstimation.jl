module CovarianceEstimation

using Statistics
import Statistics: cov
using LinearAlgebra

"""
    CovarianceEstimator

Basic type for all covariance estimators.
"""
abstract type CovarianceEstimator end

abstract type TargetType end

struct LinearShrinkage{TT<:TargetType, S<:Union{Symbol, Real}} <: CovarianceEstimator
    targetType::TT
    shrinkage::S
    function LinearShrinkage(targetType::TT, s::T) where {TT<:TargetType, T<:Real}
        @assert 0 ≤ s ≤ 1 "Shrinkage value should be between 0 and 1"
        new{TT, T}(targetType, s)
    end
    function LinearShrinkage(targetType::TT, s::Symbol=:optimal) where TT<:TargetType
        @assert s ∈ [:optimal] "Shrinkage setting not supported"
        new{TT, Symbol}(targetType, s)
    end
end

function cov(X::AbstractMatrix{<:Real}, ls::LinearShrinkage; dims::Int=1)
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

    F, ρ = targetandshrinkage(ls.targetType, S, Xc)
    ρ̂    = ls.shrinkage == :optimal ? ρ : ls.shrinkage
    return shrink(S, F, ρ̂)
end

export cov
export CovarianceEstimator, Simple,
    LedoitWolf, RaoBlackwellLedoitWolf,
    OracleApproximatingShrinkage, LW, RBLW, OAS

include("utils.jl")
include("basicmethods.jl")
include("ledoitwolf.jl")
include("chen.jl")

const LW = LinearShrinkage(LedoitWolf(), :optimal)
const RBLW = LinearShrinkage(RaoBlackwellLedoitWolf(), :optimal)
const OAS = LinearShrinkage(OracleApproximatingShrinkage(), :optimal)


end # module
