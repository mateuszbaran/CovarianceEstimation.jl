@testset "Analytical Nonlin shrinkage (ref⭒) " begin
    ANS = AnalyticalNonlinearShrinkage()
    ANS_alpha = AnalyticalNonlinearShrinkage(; alpha=1.0)

    for X̂ ∈ test_matrices
        size(X̂, 1) < 12 && continue
        ref_result = matlab_ledoitwolf_analytical_shrinkage(X̂)
        c = cov(ANS, X̂); @test c ≈ ref_result
        @test issymmetric(c)
        size(X̂, 1) < 24 && continue
        n2 = size(X̂, 1) ÷ 2
        w = FrequencyWeights(vcat(ones(n2), zeros(size(X̂, 1) - n2)))
        ref_result = matlab_ledoitwolf_analytical_shrinkage(X̂[1:n2, :])
        c = cov(ANS, X̂, w); @test c ≈ ref_result
        # Weight types besides FrequencyWeights are not supported
        aw1 = AnalyticWeights(rand(size(test_matrices[1], 1)))
        @test_throws ErrorException cov(ANS, test_matrices[1], aw1)
        @test_throws ErrorException cov(ANS, test_matrices[1], aw1; dims=2)
        @test_throws ErrorException cov(ANS, test_matrices[1], aw1; mean=nothing)

        # test alpha
        c = cov(ANS_alpha, X̂)
        @test minimum(eigvals(c)) ≈ ANS_alpha.alpha
    end
end
