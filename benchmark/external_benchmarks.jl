# ==========================================================================
# This script compares the facilities offered by CovarianceEstimation.jl
# to existing packages and compares the performances for a set of reasonably
# large matrices.
#
# Libraries considered:
#
# - [python] sklearn
# - [R] corpcor
# - [Matlab/Octave] Ledoit & Wolf's implementations
#
# ==========================================================================

using CovarianceEstimation
using Random:seed!
using DelimitedFiles
using Statistics
using LinearAlgebra

seed!(65489) # random seed for reproducibility

p = [20, 40, 200, 400]
n = [40, 20, 400, 200]
L = [randn(pᵢ, pᵢ) for pᵢ ∈ p]
C = [Lᵢ' * Lᵢ for Lᵢ ∈ L]

X = [randn(n[i], p[i]) * L[i] for i ∈ eachindex(p)]

for Cᵢ ∈ C
    pᵢ = size(Cᵢ, 1)
    writedlm("benchmark/tmpmat/C_$(pᵢ).csv", Cᵢ)
end

for Xᵢ ∈ X
    nᵢ, pᵢ = size(Xᵢ)
    writedlm("benchmark/tmpmat/X_$(nᵢ)x$(pᵢ).csv", Xᵢ)
    # octave is unstable for eigenvalue decomposition so we do everything here
    F = eigen(Xᵢ' * Xᵢ / size(Xᵢ, 1))
    writedlm("benchmark/tmpmat/S_lambda_$(nᵢ)x$(pᵢ).csv", F.values)
    writedlm("benchmark/tmpmat/S_u_$(nᵢ)x$(pᵢ).csv", F.vectors)
end

# =========
# SKLEARN
# =========
oas = LinearShrinkage(
        target=DiagonalCommonVariance(),
        shrinkage=:oas)
lw = LinearShrinkage(
        target=DiagonalCommonVariance(),
        shrinkage=:lw)
times = zeros(length(p))
res = zeros(length(p))
for i ∈ eachindex(p)
    Xᵢ = X[i]
    Cᵢ = C[i]
    start = time()
    C_oas = cov(oas, Xᵢ)
    times[i] = time() - start
    res[i] = norm(C_oas - Cᵢ)
end
res_ours = res
times_ours = times
res_sklearn = [  51.05167435  214.21678806 1602.04944525 6533.94732887][:]
times_sklearn = [0.00026703 0.00023913 0.00225997 0.00262403][:]

δ = mean(res_ours - res_sklearn)
speedup = mean(times_sklearn ./ times_ours)

for i ∈ eachindex(p)
    Xᵢ = X[i]
    Cᵢ = C[i]
    start = time()
    C_lw = cov(lw, Xᵢ)
    times[i] = time() - start
    res[i] = norm(C_lw - Cᵢ)
end
res_ours = res
times_ours = times
res_sklearn = [  50.95012365  215.37414448 1602.01890586 6534.65450682][:]
times_sklearn = [0.00076485 0.00042701 0.00249386 0.00562191][:]

δ = mean(res_ours - res_sklearn)
speedup = mean(times_sklearn ./ times_ours)


# =========
# CORPCOR
# =========
ss = LinearShrinkage(
        target=DiagonalUnequalVariance(),
        shrinkage=:ss)
times = zeros(length(p))
res = zeros(length(p))
for i ∈ eachindex(p)
    Xᵢ = X[i]
    Cᵢ = C[i]
    start = time()
    C_ss = cov(ss, Xᵢ)
    times[i] = time() - start
    res[i] = norm(C_ss - Cᵢ)
end
res_ours = res
times_ours = times
res_corpcor = [  51.39624  219.00562 1606.20378 6562.85887][:]
times_corpcor = [0.001765013, 0.001574039, 0.1169169, 0.122371]

δ = mean(res_ours .- res_corpcor)

# roughly 22x speedup
speedup = mean(times_corpcor ./ times_ours)


# ===============
# LEDOIT WOLF 1
# ===============
lw = LinearShrinkage(
        target=ConstantCorrelation(),
        shrinkage=:lw)
times = zeros(length(p))
res = zeros(length(p))
for i ∈ eachindex(p)
    Xᵢ = X[i]
    Cᵢ = C[i]
    start = time()
    C_lw = cov(lw, Xᵢ)
    times[i] = time() - start
    res[i] = norm(C_lw - Cᵢ)
end
res_ours = res
times_ours = times
res_lw1 = [  51.503, 221.721, 1606.199, 6563.666]
times_lw1 = [0.0017309, 0.0015059, 0.0107820, 0.0208941]

δ = mean(res_ours .- res_lw1)

# roughly 10x speedup
speedup = mean(times_lw1 ./ times_ours)


# ===============
# LEDOIT WOLF 2
# ===============
lw = AnalyticalNonlinearShrinkage()
times = zeros(length(p))
res = zeros(length(p))
for i ∈ [1, 3] # fat matrix don't work reliably atm, see #38
    Xᵢ = X[i]
    Cᵢ = C[i]
    start = time()
    C_lw = cov(lw, Xᵢ)
    times[i] = time() - start
    res[i] = norm(C_lw - Cᵢ)
end
res_ours = res
times_ours = times
res_lw2 = [ 49.10317, 6985.10751]
times_lw2 = [0.015798, 0.007848]

# to our advantage
δ = mean(res_ours[[1, 3]] .- res_lw2)

# roughly 25x speedup
speedup = mean(times_lw2 ./ times_ours[[1, 3]])
