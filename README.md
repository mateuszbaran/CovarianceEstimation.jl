| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl.svg?branch=master)](https://travis-ci.com/mateuszbaran/CovarianceEstimation.jl) [![Build Status](https://ci.appveyor.com/api/projects/status/fsut3j3onulvws1w?svg=true)](https://ci.appveyor.com/project/nalimilan/statsbase-jl) | [ ![codecov.io](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl/coverage.svg?branch=master)](http://codecov.io/github/mateuszbaran/CovarianceEstimation.jl?branch=master) |

# CovarianceEstimation.jl
Lightweight covariance estimation in Julia

Supported algorithms:
* Ledoit-Wolf covariance shrinkage

  O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,” The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.


* Rao-Blackwell theorem modified Ledoit-Wolf shrinkage and Oracle Approximating shrinkage

  Y. Chen, A. Wiesel, Y. C. Eldar, and A. O. Hero, “Shrinkage Algorithms for MMSE Covariance Estimation,” IEEE Transactions on Signal Processing, vol. 58, no. 10, pp. 5016–5029, Oct. 2010.
