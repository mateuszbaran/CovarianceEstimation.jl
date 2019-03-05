using BenchmarkTools
using Random
using CovarianceEstimation

const SUITE = BenchmarkGroup()

Random.seed!(1234)

const X10_100  = randn(15, 100)
const X100_100  = randn(100, 100)
const X100_1000  = randn(100, 1000)
const X100_10000  = randn(100, 10000)
const X1000_100  = randn(1000, 100)
const X1000_10000  = randn(1000, 10000)

matrices = [X10_100, X100_100, X100_1000, X100_10000, X1000_10000]
estimators = Dict(
    "Simple" => (ce = SimpleCovariance(), maxfeatures = 10000),
    "Linear constant correlation" => (ce = LinearShrinkage(ConstantCorrelation()), maxfeatures = 1000),
    "Linear diagonal unit variance" => (ce = LinearShrinkage(DiagonalUnitVariance()), maxfeatures = 1000),
    "Linear diagonal common variance" => (ce = LinearShrinkage(DiagonalCommonVariance()), maxfeatures = 1000),
    "Linear common covariance" => (ce = LinearShrinkage(CommonCovariance()), maxfeatures = 1000),
    "Linear diagonal unequal variance" => (ce = LinearShrinkage(DiagonalUnequalVariance()), maxfeatures = 1000),
    "Linear perfect positive correlation" => (ce = LinearShrinkage(PerfectPositiveCorrelation()), maxfeatures = 1000),
    "Linear RBLW" => (ce = LinearShrinkage(DiagonalCommonVariance(), :rblw), maxfeatures = 10000),
    "Linear OAS" => (ce = LinearShrinkage(DiagonalCommonVariance(), :oas), maxfeatures = 10000),
    "Analytical nonlinear shrinkage" => (ce = AnalyticalNonlinearShrinkage(), maxfeatures = 1000)
)

for (k, v) in estimators
    SUITE[k] = BenchmarkGroup()
    for m in matrices
        if size(m, 2) ≤ v[:maxfeatures]
            SUITE[k]["Size $(size(m))"] = @benchmarkable cov($v[:ce], $m)
        end
    end
end
