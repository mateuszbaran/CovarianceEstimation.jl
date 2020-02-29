# [M-estimators](@id mestim)

Example:
```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1)
X = rand(10, 10) * randn(10, 100) # a wide matrix
tyler = CovarianceEstimation.tme(X; verbose = true) # the shape matrix is 10x10
nrtmeRMT = CovarianceEstimation.nrtme(X; verbose = true)
nrtmeLW = CovarianceEstimation.nrtme(X; reg = :LW, verbose = true)
```
