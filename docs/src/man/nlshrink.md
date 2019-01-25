# [Nonlinear shrinkage estimators](@id nlshrink)

Nonlinear shrinkage estimators correspond to covariance estimators based on the eigendecomposition of the sample matrix:

```
F = eigen(X)
# ... (transformation of eigenvalues)
F.U*(d̃ .* F.U') # d̃ is a vector of transformed eigenvalues
```

Currently, only the analytical nonlinear shrinkage ([`AnalyticalNonlinearShrinkage`](@ref)) method is implemented.
