using CovarianceEstimation
using Statistics
using LinearAlgebra
using Random:seed!
using DelimitedFiles
using ProgressMeter
using PyPlot

const CE = CovarianceEstimation

seed!(6312)

# dimensions

const np = (
    # fat matrices
    (15, 20), (20, 30), (20, 50), (100, 200),
    # tall matrices
    (20, 15), (30, 20), (50, 20), (200, 100),
    )

const LSE = LinearShrinkageEstimator

const estimators = Dict(
    # BASELINE ESTIMATORS
    "s_uncorr" => Simple(),
    "s_corr"   => Simple(corrected=true),
    # LINEAR SHRINKAGE, BASED ON UNCORRECTED SAMPLE COV
    "d1v_lw_uncorr"   => LSE(target=DiagonalUnitVariance(),
                              shrinkage=:lw),
    "d1v_ss_uncorr"   => LSE(target=DiagonalUnitVariance(),
                              shrinkage=:ss),
    "dcv_lw_uncorr"   => LSE(target=DiagonalCommonVariance(),
                              shrinkage=:lw),
    "dcv_ss_uncorr"   => LSE(target=DiagonalCommonVariance(),
                              shrinkage=:ss),
    "dcv_rblw_uncorr" => LSE(target=DiagonalCommonVariance(),
                              shrinkage=:rblw),
    "dcv_oas_uncorr"  => LSE(target=DiagonalCommonVariance(),
                              shrinkage=:oas),
    "duv_lw_uncorr"   => LSE(target=DiagonalUnequalVariance(),
                              shrinkage=:lw),
    "duv_ss_uncorr"   => LSE(target=DiagonalUnequalVariance(),
                              shrinkage=:ss),
    "ccov_lw_uncorr"  => LSE(target=CommonCovariance(),
                              shrinkage=:lw),
    "ccov_ss_uncorr"  => LSE(target=CommonCovariance(),
                              shrinkage=:ss),
    "ppc_lw_uncorr"   => LSE(target=PerfectPositiveCorrelation(),
                              shrinkage=:lw),
    "ppc_ss_uncorr"   => LSE(target=PerfectPositiveCorrelation(),
                              shrinkage=:ss),
    "ccor_lw_uncorr"  => LSE(target=ConstantCorrelation(),
                              shrinkage=:lw),
    "ccor_ss_uncorr"  => LSE(target=ConstantCorrelation(),
                              shrinkage=:ss),
    # LINEAR SHRINKAGE, BASED ON CORRECTED SAMPLE COV
    "d1v_lw_corr"   => LSE(target=DiagonalUnitVariance(),
                            shrinkage=:lw, corrected=true),
    "d1v_ss_corr"   => LSE(target=DiagonalUnitVariance(),
                            shrinkage=:ss, corrected=true),
    "dcv_lw_corr"   => LSE(target=DiagonalCommonVariance(),
                            shrinkage=:lw, corrected=true),
    "dcv_ss_corr"   => LSE(target=DiagonalCommonVariance(),
                            shrinkage=:ss, corrected=true),
    "dcv_rblw_corr" => LSE(target=DiagonalCommonVariance(),
                            shrinkage=:rblw, corrected=true),
    "dcv_oas_corr"  => LSE(target=DiagonalCommonVariance(),
                            shrinkage=:oas, corrected=true),
    "duv_lw_corr"   => LSE(target=DiagonalUnequalVariance(),
                            shrinkage=:lw, corrected=true),
    "duv_ss_corr"   => LSE(target=DiagonalUnequalVariance(),
                            shrinkage=:ss, corrected=true),
    "ccov_lw_corr"  => LSE(target=CommonCovariance(),
                            shrinkage=:lw, corrected=true),
    "ccov_ss_corr"  => LSE(target=CommonCovariance(),
                            shrinkage=:ss, corrected=true),
    "ppc_lw_corr"   => LSE(target=PerfectPositiveCorrelation(),
                            shrinkage=:lw, corrected=true),
    "ppc_ss_corr"   => LSE(target=PerfectPositiveCorrelation(),
                            shrinkage=:ss, corrected=true),
    "ccor_lw_corr"  => LSE(target=ConstantCorrelation(),
                            shrinkage=:lw, corrected=true),
    "ccor_ss_corr"  => LSE(target=ConstantCorrelation(),
                            shrinkage=:ss, corrected=true),
    # NONLINEAR SHRINKAGE, BASED ON UNCORRECTED SAMPLE COV
    "anshrink_uncorr" => AnalyticalNonlinearShrinkage(),
    # NONLINEAR SHRINKAGE, BASED ON CORRECTED SAMPLE COV
    "anshrink_uncorr" => AnalyticalNonlinearShrinkage()
    )

# warmup case
begin
    println("Warmup...")
    X = randn(15, 5)
    cov(Simple(), X)
    cov(Simple(corrected=true), X)
    @showprogress for (_, estimator) ∈ estimators
        cov(estimator, X)
    end
end

# number of time each experiment is ran
nrounds = 50

results = Dict{String, Float64}()

for (n, p) ∈ np
    # Generate a colouring matrix (cholesky of true cov)
    # the first randn allows for different amplitudes, the last one allows
    # for the various dimensions to have their variance on different scales
    L = randn() .* randn(p, p) .* randn(p)
    Σ = L * L'
    p2 = p^2
    println("Dims: $n x $p...")
    @showprogress for round ∈ 1:nrounds
        # generate a sample matrix
        X = randn(n, p) * L
        # estimate the covariance
        for (name, estimator) ∈ estimators
            results[name * "_$(n)x$(p)_$round"] =
                norm(cov(estimator, X) - Σ) / p2
        end
    end
end

reskeys = collect(keys(results))

results2 = Dict{String, NTuple{2, Float64}}()

# print results for each size and every estimator (cumbersome)
for (name, estimator) ∈ estimators
    for (n, p) ∈ np
        # find all the results that match
        matchkeys = filter(k -> occursin("$(name)_$(n)x$(p)_", k), reskeys)
        r = [results[key] for key ∈ matchkeys]
        med, iqr = mean(r), -(quantile(r, [0.75, 0.25])...)
        results2["$(name)_$(n)x$(p)_"] = (med, iqr)
    end
end

res2keys = collect(keys(results2))

# print best estimator per size
for (n, p) ∈ np
    matchkeys = filter(k -> occursin("_$(n)x$(p)_", k), res2keys)
    # dirty way to find the min
    name_best = ""
    val_best = Inf
    iqr_best = Inf
    for key ∈ matchkeys
        med, iqr = results2[key]
        if med < val_best
            name_best = key
            val_best = med
            iqr_best = iqr
        end
    end
    println("Size $(n)x$(p), best: $name_best ($val_best; $iqr_best)")
end


for (n, p) ∈ np
    figure()
    rkeys = filter(κ->occursin("$(n)x$(p)_", κ), res2keys)
    data = Vector{Vector{Float64}}(undef, length(rkeys))
    for (i, k) ∈ enumerate(rkeys)
        # find the results
        matchkeys = filter(κ -> occursin(k, κ), reskeys)
        data[i] = [results[κ] for κ ∈ matchkeys]
    end
    perm = sortperm(rkeys)
    boxplot(data[perm])
    title("Size: $(n)x$(p)")
    xticks(1:length(rkeys), rkeys[perm], rotation=90)
    tight_layout()
    savefig("docs/src/assets/mse_comp/bm_$(n)x$(p).png")
end
