# [Linear shrinkage estimators](@id lshrink)

Linear shrinkage estimators correspond to covariance estimators of the form

```math
\hat\Sigma = (1-\lambda)S + \lambda F
```

where $F$ is a *target* matrix of appropriate dimensions, $\lambda\in[0,1]$ is a shrinkage intensity and $S$ is the sample covariance estimator (corrected or uncorrected depending on the `corrected` keyword).

## Targets and intensities

There are several standard target matrices (denoted by $F$) that can be used (we follow here the notations and naming conventions of [Schaffer & Strimmer 2005](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf)):


| Target name | $F_{ii}$     | $F_{ij}$ ($i\neq j$) | Comment   |
| ------------ | ------------ | ------------         | --------- |
| `DiagonalUnitVariance`| $1$ | $0$ | $F = \mathbf I$ |
| `DiagonalCommonVariance`| $v$ | 0 | $F = v\mathbf I$ |
| `DiagonalUnequalVariance` | $S_{ii}$ | 0 | $F = \mathrm{diag}(S)$, very common |
| `CommonCovariance` | $v$ | $c$ | |
| `PerfectPositiveCorrelation` | $S_{ii}$ | $\sqrt{S_{ii}S_{jj}}$ | |
| `ConstantCorrelation` | $S_{ii}$ | $\overline{r}\sqrt{S_{ii}S_{jj}}$ | used in [Ledoit & Wolf 2004](http://www.ledoit.net/honey.pdf) |


where $ v = \mathrm{tr}(S)/p $ is the average variance, $c = \sum_{i\neq j} S_{ij}/(p(p-1))$ is the average of off-diagonal terms of $S$ and $\overline{r}$ is the average of sample correlations (see [Schaffer & Strimmer 2005](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf)).

For each of these targets, an optimal shrinkage intensity $\lambda^\star$ can be computed.
A standard approach is to apply the Ledoit-Wolf formula (`shrinkage=:lw`, see [Ledoit & Wolf 2004](http://www.ledoit.net/honey.pdf)) though there are some variants that can be applied too.
Notably, Schaffer & Strimmer's variant (`shrinkage=:ss`) will ensure that the $\lambda^\star$ computed is the same for $X_c$ (the centered data matrix) as for $X_s$ (the standardised data matrix).
See [Schaffer & Strimmer 2005](http://strimmerlab.org/publications/journals/shrinkcov2005.pdf).

Chen's variant includes a Rao-Blackwellised estimator (`shrinkage=:rblw`) and an Oracle-Approximating one (`shrinkage=:oas`) for the `DiagonalCommonVariance` target.
See [Chen, Wiesel, Eldar & Hero 2010](https://arxiv.org/pdf/0907.4698.pdf).
