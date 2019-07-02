@testset "LinShrink: target F with LW (ref⭒) " begin
    lw = LinearShrinkage(ConstantCorrelation())
    testTransposition(lw, X)
    testUncorrelated(lw)
    testTranslation(lw, X)
    testDims(lw, X)
    for X̂ ∈ test_matrices
        ref_results = matlab_ledoitwolf_covcor(X̂)
        lwfixed = LinearShrinkage(ConstantCorrelation(), ref_results["shrinkage"])
        @test cov(lw, X̂) ≈ ref_results["lwcov"]
        @test cov(lwfixed, X̂) ≈ ref_results["lwcov"]
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

    ss = LinearShrinkage(target=DiagonalUnequalVariance(),
                                  shrinkage=:ss; corrected=true)
    test_mat1 = readdlm(joinpath(p, "test_matrices/20x100.csv"))
    ref_cov1  = readdlm(joinpath(p, "test_matrices/20x100_corpcor.csv"))
    test_mat2 = readdlm(joinpath(p, "test_matrices/100x20.csv"))
    ref_cov2  = readdlm(joinpath(p, "test_matrices/100x20_corpcor.csv"))
    test_mat3 = readdlm(joinpath(p, "test_matrices/50x50.csv"))
    ref_cov3  = readdlm(joinpath(p, "test_matrices/50x50_corpcor.csv"))
    @test cov(ss, test_mat1) ≈ ref_cov1
    @test cov(ss, test_mat2) ≈ ref_cov2
    @test cov(ss, test_mat3) ≈ ref_cov3
end


@testset "LinShrink: target ABCDE with LW    " begin
    # TARGET A
    lwa = LinearShrinkage(DiagonalUnitVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(SimpleCovariance(), X̂)
        Xtmp = centercols(X̂)
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum((S-Diagonal(S)).^2) + sum((diag(S).-1).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(lwa, X̂) ≈ (1.0-shrinkage) * S + shrinkage * I
        lwafixed = LinearShrinkage(DiagonalUnitVariance(), shrinkage)
        @test cov(lwafixed, X̂) ≈ (1.0 - shrinkage) * S + shrinkage * I
    end
    # TARGET B
    lwb = LinearShrinkage(DiagonalCommonVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(SimpleCovariance(), X̂)
        Xtmp = centercols(X̂)
        v = tr(S)/p
        F = v * I
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum((S-Diagonal(S)).^2) + sum((diag(S).-v).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(lwb, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwbfixed = LinearShrinkage(DiagonalCommonVariance(), shrinkage)
        @test cov(lwbfixed, X̂) ≈ (1.0 - shrinkage) * S + shrinkage * F
    end
    # TARGET C
    lwc = LinearShrinkage(CommonCovariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(SimpleCovariance(), X̂)
        Xtmp = centercols(X̂)
        v = tr(S)/p
        c = sum(S-Diagonal(S))/(p*(p-1))
        F = v * I + c * (ones(p, p) - I)
        shrinkage  = sum_var_sij(Xtmp, S, n)
        shrinkage /= sum(((S-Diagonal(S)) - c*(ones(p, p)-I)).^2) + sum((diag(S) .- v).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(lwc, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwcfixed = LinearShrinkage(CommonCovariance(), shrinkage)
        @test cov(lwcfixed, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
    # TARGET D
    lwd = LinearShrinkage(DiagonalUnequalVariance())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(SimpleCovariance(), X̂)
        Xtmp = centercols(X̂)
        F = Diagonal(S)
        shrinkage  = sum_var_sij(Xtmp, S, n, false; with_diag=false)
        shrinkage /= sum((S-Diagonal(S)).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(lwd, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwdfixed = LinearShrinkage(DiagonalUnequalVariance(), shrinkage)
        @test cov(lwdfixed, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
    # TARGET E
    lwe = LinearShrinkage(PerfectPositiveCorrelation())
    for X̂ ∈ test_matrices
        n, p = size(X̂)
        S = cov(SimpleCovariance(), X̂)
        Xtmp = centercols(X̂)
        d = diag(S)
        F = sqrt.(d*d')
        shrinkage  = sum_var_sij(Xtmp, S, n; with_diag=false)
        shrinkage -= CE.sum_fij(Xtmp, S, n, n)
        shrinkage /= sum((S - F).^2)
        shrinkage = clamp(shrinkage, 0.0, 1.0)
        @test cov(lwe, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
        lwefixed = LinearShrinkage(PerfectPositiveCorrelation(), shrinkage)
        @test cov(lwefixed, X̂) ≈ (1.0-shrinkage) * S + shrinkage * F
    end
end


@testset "LinShrink: target B with RBLW+OAS  " begin
    rblw = LinearShrinkage(DiagonalCommonVariance(), :rblw)
    testTransposition(rblw, X)
    testUncorrelated(rblw)
    testTranslation(rblw, X)
    testDims(rblw, X)

    oas = LinearShrinkage(DiagonalCommonVariance(), :oas)
    testTransposition(oas, X)
    testUncorrelated(oas)
    testTranslation(oas, X)
    testDims(oas, X)

    for X̂ ∈ test_matrices
        Ŝ_rblw = cov(rblw, X̂)
        Ŝ_oas  = cov(oas, X̂)

        X̂ = centercols(X̂)
        n, p = size(X̂)
        Ŝ    = cov(SimpleCovariance(), X̂)

        F_ref = tr(Ŝ)/p * I
        # https://arxiv.org/pdf/0907.4698.pdf eq 17
        λ_rblw_ref = ((n-2)/n*tr(Ŝ^2)+tr(Ŝ)^2)/((n+2)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        λ_rblw_ref = clamp(λ_rblw_ref, 0.0, 1.0)
        # https://arxiv.org/pdf/0907.4698.pdf eq 23
        λ_oas_ref = ((1-2/p)*tr(Ŝ^2)+tr(Ŝ)^2)/((n+1-2/p)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        λ_oas_ref = clamp(λ_oas_ref, 0.0, 1.0)

        @test cov(rblw, X̂) ≈ CE.linshrink(F_ref, Ŝ, λ_rblw_ref)
        @test cov(oas, X̂) ≈ CE.linshrink(F_ref, Ŝ, λ_oas_ref)

        rblw_fixed = LinearShrinkage(DiagonalCommonVariance(), λ_rblw_ref)
        oas_fixed = LinearShrinkage(DiagonalCommonVariance(), λ_oas_ref)
        @test cov(rblw_fixed, X̂) ≈ CE.linshrink(F_ref, Ŝ, λ_rblw_ref)
        @test cov(oas_fixed, X̂) ≈ CE.linshrink(F_ref, Ŝ, λ_oas_ref)
    end
end


@testset "LinShrink: all targets, std (SS)   " begin
    for target ∈ [
            DiagonalUnitVariance(),
            DiagonalCommonVariance(),
            DiagonalUnequalVariance(),
            CommonCovariance(),
            PerfectPositiveCorrelation(),
            ConstantCorrelation()
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
            LSEss = LinearShrinkage(target, :ss, corrected=c)
            LSElw = LinearShrinkage(target, :lw, corrected=c)
            @test cov(LSEss, Xcs) ≈ cov(LSElw, Xcs)
            # weights are not currently exported but it seems to be a good
            # idea to ensure proper errors
            fw1 = FrequencyWeights(rand(1:10, size(Xcs, 1)))
            @test_throws ErrorException cov(LSEss, Xcs, fw1)
            @test_throws ErrorException cov(LSEss, Xcs, fw1; dims=2)
            @test_throws ErrorException cov(LSEss, Xcs, fw1; mean=nothing)
        end
    end
end

@testset "LinShrink: all targets, mean arg   " begin
    for target ∈ [
            DiagonalUnitVariance(),
            DiagonalCommonVariance(),
            DiagonalUnequalVariance(),
            CommonCovariance(),
            PerfectPositiveCorrelation(),
            ConstantCorrelation()
            ]
        for c ∈ [true, false], dim ∈ [1, 2], X̂ ∈ test_matrices[1:3]
            Xcs = X̂
            Xcs = centercols(Xcs)
            n = size(Xcs, dim)
            LSEss = LinearShrinkage(target, :ss, corrected=c)
            LSElw = LinearShrinkage(target, :lw, corrected=c)
            m = mean(cov(LSEss, Xcs, dims=dim, mean=nothing))
            @test cov(LSEss, Xcs, dims=dim, mean=mean(Xcs, dims=dim)) ≈ cov(LSEss, Xcs, dims=dim, mean=nothing)
            @test cov(LSElw, Xcs, dims=dim, mean=mean(Xcs, dims=dim)) ≈ cov(LSElw, Xcs, dims=dim, mean=nothing)
            if dim == 2
                meanvec = vec(mean(Xcs, dims=dim))
                @test cov(LSEss, Xcs, dims=dim, mean=meanvec) ≈ cov(LSEss, Xcs, dims=dim, mean=nothing)
                @test cov(LSElw, Xcs, dims=dim, mean=meanvec) ≈ cov(LSElw, Xcs, dims=dim, mean=nothing)
            end
        end
    end
end
