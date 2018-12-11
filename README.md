| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl.svg?branch=master)](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl) [![Build status](https://ci.appveyor.com/api/projects/status/7riq3mtk8wy6k3yl?svg=true)](https://ci.appveyor.com/project/mateuszbaran/covarianceestimation-jl) | [ ![codecov.io](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl/coverage.svg?branch=master)](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl?branch=master) |

# CovarianceEstimation.jl
Lightweight covariance estimation in Julia.
The package is currently unregistered but can be installed with `Pkg` using

```julia-repl
] add https://github.com/mateuszbaran/CovarianceEstimation.jl
```

## Quick start

```julia
using CovarianceEstimation
using Random

X = randn(5, 3)

S_uncorrected  = cov(Simple(), X)
S_corrected    = cov(Corrected(), X)

# using optimal shrinkage
S_ledoitwolf   = cov(LedoitWolf(), X)
S_rbledoitwolf = cov(RaoBlackwellLedoitWolf(), X)
S_oracleapprox = cov(OracleApproximatingShrinkage(), X)

# using a given shrinkage
S_ledoitwolf_05 = cov(LedoitWolf(0.5), X)
```

## Currently supported algorithms

* Basic corrected and uncorrected sample covariance (via the `Statistics` package),
* Ledoit-Wolf shrinkage [**1**],
* Rao-Blackwell Ledoit-Wolf shrinkage and Oracle Approximating shrinkage [**2**].

## References

* [**1**] O. Ledoit and M. Wolf, *[Honey, I Shrunk the Sample Covariance Matrix](http://www.ledoit.net/honey.pdf)*, The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
* [**2**] Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero, *[Shrinkage Algorithms for MMSE Covariance Estimation](https://arxiv.org/pdf/0907.4698.pdf)*, IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
