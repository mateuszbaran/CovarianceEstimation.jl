module CovarianceEstimation

using Statistics
import Statistics: cov

"""
    CovarianceEstimator

Basic type for all covariance estimators.
"""
abstract type CovarianceEstimator end

export cov
export CovarianceEstimator, SimpleCovariance, CorrectedCovariance,
    LedoitWolfCovariance, RBLWCovariance, OASCovariance

include("basicmethods.jl")
include("ledoitwolf.jl")
include("chen.jl")

end # module
