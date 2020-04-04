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
    AnalyticalNonlinearShrinkage,
    TylerMEstimator, NormalizedRegularizedTylerMEstimator


include("utils.jl")
include("linearshrinkage.jl")
include("nonlinearshrinkage.jl")
include("m-estimators.jl")

end # module
