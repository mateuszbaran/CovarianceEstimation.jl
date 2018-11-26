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
    LedoitWolfCovariance

include("basicmethods.jl")
include("ledoitwolf.jl")

end # module
