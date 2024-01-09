module CovarianceEstimation

using Statistics
using StatsBase
using LinearAlgebra
import StatsBase: cov
using WoodburyMatrices
using TSVD

export cov
export CovarianceEstimator, SimpleCovariance,
    LinearShrinkage,
    # Targets for linear shrinkage
    DiagonalUnitVariance, DiagonalCommonVariance, DiagonalUnequalVariance,
    CommonCovariance, PerfectPositiveCorrelation, ConstantCorrelation,
    # Eigendecomposition-based methods
    AnalyticalNonlinearShrinkage,
    # Biweight midcovariance
    BiweightMidcovariance,
    # Woodbury-based methods
    WoodburyEstimator,
    # Loss functions
    NormLossCov, StatLossCov


include("utils.jl")
include("loss.jl")
include("biweight.jl")
include("linearshrinkage.jl")
include("nonlinearshrinkage.jl")
include("woodbury.jl")

end # module
