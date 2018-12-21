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
    X_s &=& X_c D^{-1}
\end{eqnarray}
```

## Simple estimator

The standard covariance estimator is easily obtained via \eqref{simple-covariance}.
It can be specified with the constructor `Simple` which can take a named argument `corrected` (either `false` (default) or `true`).

```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1)
n, p = 5, 7
X = randn(n, p)
# corrected covariance
S = cov(X, Simple(corrected=true))
# we can also manually compute it and compare
Xc = (X .- sum(X, dims=1)/n) # centering
κ = n-1
S ≈ (Xc'*Xc)/κ
```

## Linear shrinkage estimators

Linear shrinkage estimators correspond to covariance estimators of the form

```math
C = (1-\lambda)S + \lambda F
```

where $F$ is a *target* matrix of appropriate dimensions and $\lambda\in[0,1]$ is a shrinkage intensity.
There are several standard targets that can be used, a simple example being the identity matrix.

The shrinkage intensity $\lambda$ can be specified manually or computed automatically.
Depending on the target, different approaches are implemented to compute a good intensity such as, for example, the Ledoit-Wolfe optimal intensity (which is the default intensity if you don't specify it).

You can read more on the targets that can be used and the corresponding automatic intensities [here](@ref lshrink).

Here is an example using the identity matrix as a target and automatic shrinkage intensity (Ledoit-Wolfe):

```@example
using CovarianceEstimation # hide
using Random # hide
Random.seed!(1)
n, p = 2, 3
X = randn(n, p)
target = DiagonalUnitVariance()
shrinkage = :lw # Ledoit-Wolfe optimal shrinkage
method = LinearShrinkageEstimator(target, shrinkage)
cov(X, method)
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
method2 = LinearShrinkageEstimator(target, shrinkage)
cov(X, method2)
```

## Nonlinear shrinkage estimators

[More on linear shrinkage estimators...](@ref nlshrink)
