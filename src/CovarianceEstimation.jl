module CovarianceEstimation

using Statistics
import Statistics: cov
using LinearAlgebra

"""
    CovarianceEstimator

Basic type for all covariance estimators.
"""
abstract type CovarianceEstimator end

export cov
export CovarianceEstimator, Simple, Corrected,
    LedoitWolf, RaoBlackwellLedoitWolf,
    OracleApproximatingShrinkage

include("utils.jl")
include("basicmethods.jl")
include("ledoitwolf.jl")
include("chen.jl")

end # module
