# [Linear shrinkage estimators](@id lshrink)

Linear shrinkage estimators correspond to covariance estimators of the form

```math
\hat\Sigma = (1-\lambda)S + \lambda F
```

where $F$ is a *target* matrix of appropriate dimensions, $\lambda\in[0,1]$ is a shrinkage intensity and $S$ is the sample covariance estimator.

## Targets and intensities

There are several standard targets that can be used (we follow here [Schaffer & Strimmer 2005](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf)):

1. `DiagonalUnitVariance` where $F=I$ the identity matrix,
1. `DiagonalCommonVariance` where $F=vI$
1. `DiagonalUnequalVariance` where $F=\mathrm{diag}(S)$
1. `CommonCovariance` where $F_{ii}=v$ and $F_{ij}=c$
1. `PerfectPositiveCorrelation` where $F_{ii}=S_{ii}$ and $F_{ij}=\sqrt{S_{ii}S_{jj}}$
1. `ConstantCorrelation` where $F_{ii}=S_{ii}$ and $F_{ij}=\overline{r}\sqrt{S_{ii}S_{jj}}$

where $ v = \mathrm{tr}(S)/p $ is the average variance, $c = \sum_{i\neq j} S_{ij}/(p*(p-1))$ is the average of off-diagonal terms of $S$ and $\overline{r}$ is the average of sample correlations.

For each of these targets, an optimal shrinkage intensity $\lambda^\star$ can be computed.
A standard approach is to apply the Ledoit-Wolfe formula (`shrinkage=:lw`) though there are some variants that can be applied too.
See [Ledoit & Wolfe 2004](http://www.ledoit.net/honey.pdf).

Notably, Schaffer & Strimmer's variant (`shrinkage=:ss`) will ensure that the $\lambda^\star$ computed is the same for $X_c$ (the centered data matrix) as for $X_s$ (the standardised data matrix).
See [Schaffer & Strimmer 2005](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf).

Chen's variant includes a Rao-Blackwellised estimator (`shrinkage=:rblw`) and an Oracle-Approximating one (`shrinkage=:oas`) for the `DiagonalCommonVariance` target.
See [Chen, Wiesel, Eldar & Hero 2010](https://arxiv.org/pdf/0907.4698.pdf).

## What to use

In general, all linear shrinkage estimators offer roughly similar accuracies though the `DiagonalCommonVariance` target and the `CommonCovariance` with `:ss` shrinkage seem to perform marginally better in our experiments.
