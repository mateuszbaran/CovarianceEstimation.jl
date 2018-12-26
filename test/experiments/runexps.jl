using CovarianceEstimation
using Statistics
using LinearAlgebra
using Random:seed!
using DelimitedFiles
using ProgressMeter

const CE = CovarianceEstimation

seed!(6312)

# dimensions

np = (
    # fat matrices
    (15, 20), (20, 30), (20, 50), (100, 1000),
    # tall matrices
    (20, 15), (30, 20), (50, 20), (1000, 100),
    )

const LSE = LinearShrinkageEstimator

const estimators = Dict(
    "d1v_lw"        => LSE(target=DiagonalUnitVariance(),      shrinkage=:lw),
    "d1v_ss"        => LSE(target=DiagonalUnitVariance(),      shrinkage=:ss),
    "dcv_lw"        => LSE(target=DiagonalCommonVariance(),    shrinkage=:lw),
    "dcv_ss"        => LSE(target=DiagonalCommonVariance(),    shrinkage=:ss),
    "dcv_rblw"      => LSE(target=DiagonalCommonVariance(),    shrinkage=:rblw),
    "dcv_oas"       => LSE(target=DiagonalCommonVariance(),    shrinkage=:oas),
    "duv_lw"        => LSE(target=DiagonalUnequalVariance(),   shrinkage=:lw),
    "duv_ss"        => LSE(target=DiagonalUnequalVariance(),   shrinkage=:ss),
    "ccov_lw"       => LSE(target=CommonCovariance(),          shrinkage=:lw),
    "ccov_ss"       => LSE(target=CommonCovariance(),          shrinkage=:ss),
    "ppc_lw"        => LSE(target=PerfectPositiveCorrelation(),shrinkage=:lw),
    "ppc_ss"        => LSE(target=PerfectPositiveCorrelation(),shrinkage=:ss),
    "ccor_lw"       => LSE(target=ConstantCorrelation(),       shrinkage=:lw),
    "ccor_ss"       => LSE(target=ConstantCorrelation(),       shrinkage=:ss),
    "anshrink"      => AnalyticalNonlinearShrinkage()
    )

# warmup case
begin
    println("Warmup...")
    X = randn(15, 5)
    cov(X, Simple())
    cov(X, Simple(corrected=true))
    @showprogress for (_, estimator) ∈ estimators
        cov(X, estimator)
    end
end

# number of time each experiment is ran
nrounds = 50

results = Dict{String, Float64}()

@showprogress for (n, p) ∈ np
    # Generate a colouring matrix (cholesky of true cov)
    L = randn() .* randn(p, p)
    Σ = L * L'
    p2 = p^2
    for round ∈ 1:nrounds
        # generate a sample matrix
        X = randn(n, p) * L
        # baseline estimators
        results["s_uncorr_$(n)x$(p)_$round"] =
            norm(cov(X, Simple()) - Σ) / p2
        results["s_corr_$(n)x$(p)_$round"] =
            norm(cov(X, Simple(true)) - Σ) / p2
        # estimate the covariance
        for (name, estimator) ∈ estimators
            if estimator isa LSE
                results[name * "_uncorr_$(n)x$(p)_$round"] =
                    norm(cov(X, estimator, corrected=false) - Σ) / p2
                results[name * "_corr_$(n)x$(p)_$round"] =
                    norm(cov(X, estimator, corrected=true) - Σ) / p2
            else
                # XXX this will be fixed with #37
                results[name * "_uncorr_$(n)x$(p)_$round"] =
                    norm(cov(X, estimator) - Σ) / p2
            end
        end
    end
end

reskeys = collect(keys(results))

results2 = Dict{String, NTuple{2, Float64}}()

# print results for each size and every estimator (cumbersome)
for (name, estimator) ∈ estimators
    if estimator isa LSE
        println("Res for $(name)_uncorr")
        for (n, p) ∈ np
            # find all the results that match
            matchkeys = filter(k -> occursin("$(name)_uncorr_$(n)x$(p)_", k),
                               reskeys)
            r = [results[key] for key ∈ matchkeys]
            μ = mean(r)
            σ = std(r)
            results2["$(name)_uncorr_$(n)x$(p)_"] = (μ, σ)
            println("$(n)x$(p) -- $μ ($σ)")
        end
        println("Res for $(name)_corr")
        for (n, p) ∈ np
            # find all the results that match
            matchkeys = filter(k -> occursin("$(name)_corr_$(n)x$(p)_", k),
                               reskeys)
            r = [results[key] for key ∈ matchkeys]
            μ = mean(r)
            σ = std(r)
            results2["$(name)_corr_$(n)x$(p)_"] = (μ, σ)
            println("$(n)x$(p) -- $μ ($σ)")
        end
    else
        println("Res for $(name)_uncorr")
        for (n, p) ∈ np
            # find all the results that match
            matchkeys = filter(k -> occursin("$(name)_uncorr_$(n)x$(p)_", k),
                               reskeys)
            r = [results[key] for key ∈ matchkeys]
            μ = mean(r)
            σ = std(r)
            results2["$(name)_uncorr_$(n)x$(p)_"] = (μ, σ)
            println("$(n)x$(p) -- $μ ($σ)")
        end
    end
end

res2keys = collect(keys(results2))

# print best estimator per size
for (n, p) ∈ np
    matchkeys = filter(k -> occursin("_$(n)x$(p)_", k), res2keys)
    # dirty way to find the min
    name_best = ""
    val_best = Inf
    sig_best = Inf
    for key ∈ matchkeys
        μ, σ = results2[key]
        if μ < val_best
            name_best = key
            val_best = μ
            sig_best = σ
        end
    end
    println("Size $(n)x$(p), best: $name_best ($val_best; $sig_best)")
end
