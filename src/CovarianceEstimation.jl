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
    # Biweight midcovariance
    BiweightMidcovariance


include("utils.jl")
include("biweight.jl")
include("linearshrinkage.jl")
include("nonlinearshrinkage.jl")

end # module
