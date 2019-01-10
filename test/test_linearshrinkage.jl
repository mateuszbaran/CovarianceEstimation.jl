@testset "LinShrink: target F with LW (ref⭒) " begin
    lw = LinearShrinkageEstimator(ConstantCorrelation())
    testTransposition(lw, X)
    testUncorrelated(lw)
    testTranslation(lw, X)
    for X̂ ∈ test_matrices
        ref_results = matlab_ledoitwolf_covcor(X̂)
        lwfixed = LinearShrinkageEstimator(ConstantCorrelation(), ref_results["shrinkage"])
        @test cov(X̂, lw) ≈ ref_results["lwcov"]
        @test cov(X̂, lwfixed) ≈ ref_results["lwcov"]
    end
end


@testset "LinShrink: target D with SS (ref⭒) " begin
    p = ifelse(endswith(pwd(), "CovarianceEstimation"), "test", "")
    ## R Script used to compare:
    # require(corpcor)
    # tm1 = read.table("20x100.csv")
    # tm2 = read.table("100x20.csv")
    # tm3 = read.table("50x50.csv")
    # c1 = cov.shrink(tm1, lambda.var=0.0)
    # c2 = cov.shrink(tm2, lambda.var=0.0)
    # c3 = cov.shrink(tm3, lambda.var=0.0)

    ss = LinearShrinkageEstimator(target=DiagonalUnequalVariance(),
                                  shrinkage=:ss; corrected=true)
    test_mat1 = readdlm(joinpath(p, "test_matrices/20x100.csv"))
    ref_cov1  = readdlm(joinpath(p, "test_matrices/20x100_corpcor.csv"))
    test_mat2 = readdlm(joinpath(p, "test_matrices/100x20.csv"))
    ref_cov2  = readdlm(joinpath(p, "test_matrices/100x20_corpcor.csv"))
    test_mat3 = readdlm(joinpath(p, "test_matrices/50x50.csv"))
    ref_cov3  = readdlm(joinpath(p, "test_matrices/50x50_corpcor.csv"))
    @test cov(test_mat1, ss) ≈ ref_cov1
    @test cov(test_mat2, ss) ≈ ref_cov2
    @test cov(test_mat3, ss) ≈ ref_cov3
end


@testset "LinShrink: target ABCDE with LW    " begin
    # TARGET A
    lwa = LinearShrinkageEstimator(DiagonalUnitVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(X̂, Simple())
        Xtmp = centercols(X̂)
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum((S-Diagonal(S)).^2) + sum((diag(S).-1).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(X̂, lwa) ≈ (1.0-shrinkage) * S + shrinkage * I
        lwafixed = LinearShrinkageEstimator(DiagonalUnitVariance(), shrinkage)
        @test cov(X̂, lwafixed) ≈ (1.0 - shrinkage) * S + shrinkage * I
    end
    # TARGET B
    lwb = LinearShrinkageEstimator(DiagonalCommonVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(X̂, Simple())
        Xtmp = centercols(X̂)
        v = tr(S)/p
        F = v * I
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum((S-Diagonal(S)).^2) + sum((diag(S).-v).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(X̂, lwb) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwbfixed = LinearShrinkageEstimator(DiagonalCommonVariance(), shrinkage)
        @test cov(X̂, lwbfixed) ≈ (1.0 - shrinkage) * S + shrinkage * F
    end
    # TARGET C
    lwc = LinearShrinkageEstimator(CommonCovariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(X̂, Simple())
        Xtmp = centercols(X̂)
        v = tr(S)/p
        c = sum(S-Diagonal(S))/(p*(p-1))
        F = v * I + c * (ones(p, p) - I)
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum(((S-Diagonal(S)) - c*(ones(p, p)-I)).^2) + sum((diag(S) .- v).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(X̂, lwc) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwcfixed = LinearShrinkageEstimator(CommonCovariance(), shrinkage)
        @test cov(X̂, lwcfixed) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
    # TARGET D
    lwd = LinearShrinkageEstimator(DiagonalUnequalVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(X̂, Simple())
        Xtmp = centercols(X̂)
        F = Diagonal(S)
        shrinkage  = sum_var_sij(Xtmp, S, n, false; with_diag=false)
        shrinkage /= sum((S-Diagonal(S)).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(X̂, lwd) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwdfixed = LinearShrinkageEstimator(DiagonalUnequalVariance(), shrinkage)
        @test cov(X̂, lwdfixed) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
    # TARGET E
    lwe = LinearShrinkageEstimator(PerfectPositiveCorrelation())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(X̂, Simple())
        Xtmp = centercols(X̂)
        d = diag(S)
        F = sqrt.(d*d')
        shrinkage  = sum_var_sij(Xtmp, S, n; with_diag=false)
        shrinkage -= CE.sum_fij(Xtmp, S, n, n)
        shrinkage /= sum((S - F).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(X̂, lwe) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwefixed = LinearShrinkageEstimator(PerfectPositiveCorrelation(), shrinkage)
        @test cov(X̂, lwefixed) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
end


@testset "LinShrink: target B with RBLW+OAS  " begin
    rblw = LinearShrinkageEstimator(DiagonalCommonVariance(), :rblw)
    testTransposition(rblw, X)
    testUncorrelated(rblw)
    testTranslation(rblw, X)

    oas = LinearShrinkageEstimator(DiagonalCommonVariance(), :oas)
    testTransposition(oas, X)
    testUncorrelated(oas)
    testTranslation(oas, X)

    for X̂ ∈ test_matrices
        Ŝ_rblw = cov(X̂, rblw)
        Ŝ_oas  = cov(X̂, oas)

        X̂ = centercols(X̂)
        n, p = size(X̂)
        Ŝ    = cov(X̂, Simple())

        F_ref = tr(Ŝ)/p * I
        # https://arxiv.org/pdf/0907.4698.pdf eq 17
        λ_rblw_ref = ((n-2)/n*tr(Ŝ^2)+tr(Ŝ)^2)/((n+2)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        λ_rblw_ref = clamp(λ_rblw_ref, 0.0, 1.0)
        # https://arxiv.org/pdf/0907.4698.pdf eq 23
        λ_oas_ref = ((1-2/p)*tr(Ŝ^2)+tr(Ŝ)^2)/((n+1-2/p)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        λ_oas_ref = clamp(λ_oas_ref, 0.0, 1.0)

        @test cov(X̂, rblw) ≈ CE.linshrink(F_ref, Ŝ, λ_rblw_ref)
        @test cov(X̂, oas) ≈ CE.linshrink(F_ref, Ŝ, λ_oas_ref)

        rblw_fixed = LinearShrinkageEstimator(DiagonalCommonVariance(), λ_rblw_ref)
        oas_fixed = LinearShrinkageEstimator(DiagonalCommonVariance(), λ_oas_ref)
        @test cov(X̂, rblw_fixed) ≈ CE.linshrink(F_ref, Ŝ, λ_rblw_ref)
        @test cov(X̂, oas_fixed) ≈ CE.linshrink(F_ref, Ŝ, λ_oas_ref)
    end
end


@testset "LinShrink: all targets, std (SS)   " begin
    for target ∈ [
            DiagonalUnitVariance(),
            DiagonalCommonVariance(),
            DiagonalUnequalVariance(),
            CommonCovariance(),
            PerfectPositiveCorrelation(),
            # ConstantCorrelation()
            ]
        #=
        :ss assumes that the shrinkage should be the same if the data
        matrix is standardised. This is what is tested.
        =#
        for c ∈ [false, true], X̂ ∈ test_matrices
            Xcs = X̂
            Xcs = centercols(Xcs)
            n = size(Xcs, 1)
            for i ∈ 1:size(Xcs, 2)
                Xcs[:, i] ./= sqrt(var(Xcs[:, i], corrected=c))
            end
            LSEss = LinearShrinkageEstimator(target, :ss, corrected=c)
            LSElw = LinearShrinkageEstimator(target, :lw, corrected=c)
            @test cov(Xcs, LSEss) ≈ cov(Xcs, LSElw)
        end
    end
end
