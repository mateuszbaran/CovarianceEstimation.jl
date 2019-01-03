| Status | Coverage | Docs |
| :----: | :----: | :----: |
| [![Build Status](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl.svg?branch=master)](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl) [![Build status](https://ci.appveyor.com/api/projects/status/7riq3mtk8wy6k3yl?svg=true)](https://ci.appveyor.com/project/mateuszbaran/covarianceestimation-jl) | [ ![codecov.io](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl/coverage.svg?branch=master)](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl?branch=master) | [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mateuszbaran.github.io/CovarianceEstimation.jl/dev) |

# CovarianceEstimation.jl
Lightweight robust covariance estimation in Julia.
The package is currently unregistered but can be installed with `Pkg` using

```julia-repl
] add https://github.com/mateuszbaran/CovarianceEstimation.jl
```

## Quick start

```julia
using CovarianceEstimation

X = randn(5, 7)

S_uncorrected  = cov(X, Simple())
S_corrected    = cov(X, Simple(corrected=true))

# using linear shrinkage with different targets
LSE = LinearShrinkageEstimator
# - Ledoit-Wolf target + shrinkage
method = LSE(ConstantCorrelation())
S_ledoitwolf = cov(X, method)
# - Chen target + shrinkage (using the more verbose call)
method = LSE(target=DiagonalCommonVariance(), shrinkage=:rblw)
S_chen_rblw = cov(X, method)
method = LSE(target=DiagonalCommonVariance(), shrinkage=:oas)
S_chen_oas = cov(X, method)

# a pre-defined shrinkage can be used as well
method = LinearShrinkageEstimator(DiagonalUnitVariance(), 0.5)
# using a given shrinkage
S_05 = cov(X, method)
```

## Currently supported algorithms

* `Simple`: basic corrected and uncorrected sample covariance (via the `Statistics` package),
* `LinearShrinkageEstimator`: James-Stein type estimator of the form `(1-λ)S+λF` where `S` is the uncorrected simple estimator, `F` is a target and `λ∈[0,1]` a shrinkage intensity.
  - common targets are implemented following the taxonomy given in [**1**] along with Ledoit-Wolf optimal shrinkage intensities [**2**].
  - in the case of the `DiagonalCommonVariance` target, a Rao-Blackwellised intensity and Oracle-Approximating intensity are also supported (see [**3**]).
  - **Note**: `S` is symmetric semi-positive definite so that if the `F` is symmetric positive definite and provided `λ` is non-zero, the estimator obtained after shrinkage is also symmetric positive definite. For the diagonal targets `DiagonalUnitVariance`, `DiagonalCommonVariance` and `DiagonalUnequalVariance` the target is necessarily SPD.

## Comparison to existing libraries

Rough benchmarks are run over random matrices of various sizes (`40x20, 20x40, 400x200, 200x400`).
These benchmarks should (as usual) be taken with a pinch of salt but essentially a significant speedup should be expected for a standard problem.

* **Sklearn** (implements `DiagonalCommonVariance` target with `oas` and `lw` shrinkage)
  - average speedup: `5x`
* **Corpcor** (implements `DiagonalUnequalVariance` target with `ss` shrinkage)
  - average speedup: `22x`
* **Ledoit-Wolfe 1** (implements `ConstantCorrelation` target with `lw` shrinkage, we used Octave for the comparison)
  - average speedup: `12x`
* **Ledoit-Wolfe 2** (implements `Nonlinear shrinkage`)
  - average speedup: `25x`


## References

* [**1**] J. Schäfer and K. Strimmer, *[A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf)*, Statistical Applications in Genetics and Molecular Biology, 2005.
* [**2**] O. Ledoit and M. Wolf, *[Honey, I Shrunk the Sample Covariance Matrix](http://www.ledoit.net/honey.pdf)*, The Journal of Portfolio Management, 2004.
* [**3**] Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero, *[Shrinkage Algorithms for MMSE Covariance Estimation](https://arxiv.org/pdf/0907.4698.pdf)*, IEEE Transactions on Signal Processing, 2010.
