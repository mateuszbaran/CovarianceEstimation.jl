module CovarianceEstimation

using Statistics
using StatsBase
using LinearAlgebra
import StatsBase: cov

export cov
export CovarianceEstimator, SimpleCovariance,
    LinearShrinkage,
    # Targets for linear shrinkage
    DiagonalUnitVariance, DiagonalCommonVariance, DiagonalUnequalVariance,
    CommonCovariance, PerfectPositiveCorrelation, ConstantCorrelation,
    # Eigendecomposition-based methods
    AnalyticalNonlinearShrinkage


include("utils.jl")
include("linearshrinkage.jl")
include("nonlinearshrinkage.jl")

end # module
