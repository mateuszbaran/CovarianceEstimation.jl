# Methods

## Notations

In these docs, $X$ denotes an $n\times p$ data matrix describing $n$ observations with $p$ features (variables) and possibly $p > n$.
For most estimators, the matrix is assumed to have entries in $\mathbb R$.

!!! note

    While in the docs we assume that the rows of $X$ correspond to observations and the columns to features, in the code you can specify this using the `dims` keyword with `dims=1` being the default whereas `dims=2` will take rows as features and columns as observations.

We will write $X_c$ the *centered* matrix i.e. where each column sums up to zero.
We will also write $X_s$ the *standardised* matrix i.e. where each column not only sums up to zero but is scaled to have sample variance one.
Finally, we will write $S$ the standard sample covariance estimator (see below) of size $p\times p$ and $D$ the diagonal matrix matching the diagonal of $S$.

Note that the sample variance and covariance can be corrected or uncorrected (and this can be specified in the code).
In order to avoid having to specify this everywhere in the document, it is useful to introduce a last symbol: $\kappa$ which is either set to $n$ (uncorrected case) or $(n-1)$ (corrected case).

With these notations we can write:

```math
\begin{eqnarray}
    S &=& \kappa^{-1}X_c^T X_c \label{simple-covariance}\\
    D_{ij} &=& S_{ij}\mathbf 1_{[i=j]}\\
    X_s &=& X_c D^{-1/2}
\end{eqnarray}
```

## Simple estimator

The standard covariance estimator is easily obtained via \eqref{simple-covariance}.
It can be specified with the constructor `SimpleCovariance` which can take a named argument `corrected` (either `false` (default) or `true`).

```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1)
n, p = 5, 7
X = randn(n, p)
# corrected covariance
S = cov(SimpleCovariance(corrected=true), X)
# we can also manually compute it and compare
Xc = (X .- sum(X, dims=1)/n) # centering
κ = n-1 # correction factor
S ≈ (Xc'*Xc)/κ
```

## Linear shrinkage estimators

Linear shrinkage estimators correspond to covariance estimators of the form

```math
\hat\Sigma = (1-\lambda)S + \lambda F
```

where $F$ is a *target* matrix of appropriate dimensions, $\lambda\in[0,1]$ is a shrinkage intensity and $S$ is the sample covariance estimator.
There are several standard targets that can be used, a simple example being the identity matrix.

The shrinkage intensity $\lambda$ can be specified manually or computed automatically.
Depending on the target, different approaches are implemented to compute a good intensity such as, for example, the Ledoit-Wolf optimal intensity (which is the default intensity if you don't specify it).

You can read more on the targets that can be used and the corresponding automatic intensities [here](@ref lshrink).

Here is an example using the identity matrix as a target and automatic shrinkage intensity (Ledoit-Wolfe):

```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1)
n, p = 2, 3
X = randn(n, p)
target = DiagonalUnitVariance()
shrinkage = :lw # Ledoit-Wolf optimal shrinkage
method = LinearShrinkage(target, shrinkage)
cov(method, X)
```

You can also specify the intensity manually:

```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1) # hide
n, p = 2, 3 # hide
X = randn(n, p) # hide
target = DiagonalUnitVariance() # hide
shrinkage = 0.8
method2 = LinearShrinkage(target, shrinkage)
cov(method2, X)
```

[Read more on linear shrinkage estimators...](@ref lshrink)

## Nonlinear shrinkage estimators

[Read more on nonlinear shrinkage estimators...](@ref nlshrink)

## Comparing estimators

You may want to look at our simple [comparison of covariance estimators](@ref msecomp) which compares the MSE of the various estimators in a range of situations.
Long story short, the `LinearShrinkageEstimator` with `DiagonalUnequalVariance` target performs well in the case $n<p$ though most other estimators don't fare too badly in comparison.
In the case $n>p$, the nonlinear shrinkage method does very well (though it is more expensive to compute).
