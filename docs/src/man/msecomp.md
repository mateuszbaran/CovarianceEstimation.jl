# [MSE Comparison](@id msecomp)

Below are results obtained with a variety of data matrices of dimensions $n\times p$.
For each pair of dimension, 50 covariance matrices are generated with associated sample data matrices.
The covariance obtained with the different estimators are then compared to the ground-truth and the MSE is reported.

| Abbreviation | Method |
| ---- | ---- |
| `anshrink` | analytical nonlinear shrinkage |
| `ccor` | LSE with constant correlation target |
| `ccov` | LSE with constant covariance target |
| `d1v`  | LSE with identity target |
| `dcv`  | LSE with diagonal common variance target |
| `duv`  | LSE with diagonal unequal variance target |
| `ppc`  | LSE with perfect positive correlation target |
| `s`  | Simple estimator (baseline) |
| `_lw`  | uses ledoit-wolf shrinkage |
| `_ss`  | uses schaffer-strimmer shrinkage |
| `_oas`  | uses oracle approximating shrinkage |
| `_rblw`  | uses rao-blackwellised ledoit-wolf shrinkage |


## Fat matrices

![](../assets/mse_comp/bm_15x20.png)
![](../assets/mse_comp/bm_20x30.png)
![](../assets/mse_comp/bm_20x50.png)
![](../assets/mse_comp/bm_100x200.png)

## Tall matrices

![](../assets/mse_comp/bm_20x15.png)
![](../assets/mse_comp/bm_30x20.png)
![](../assets/mse_comp/bm_50x20.png)
![](../assets/mse_comp/bm_200x100.png)
