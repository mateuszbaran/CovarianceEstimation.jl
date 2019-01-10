module CovarianceEstimation

using Statistics
using StatsBase
using LinearAlgebra
import Statistics: cov

"""
    CovarianceEstimator

Basic type for all covariance estimators.
"""
abstract type CovarianceEstimator end

export cov
export CovarianceEstimator, Simple,
    LinearShrinkage,
    # Targets for linear shrinkage
    DiagonalUnitVariance, DiagonalCommonVariance, DiagonalUnequalVariance,
    CommonCovariance, PerfectPositiveCorrelation, ConstantCorrelation,
    # Eigendecomposition-based methods
    AnalyticalNonlinearShrinkage


include("utils.jl")
include("simplecov.jl")
include("linearshrinkage.jl")
include("nonlinearshrinkage.jl")

end # module
