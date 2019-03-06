@testset "Analytical Nonlin shrinkage (ref⭒) " begin
    ANS = AnalyticalNonlinearShrinkage()
    for X̂ ∈ test_matrices
        size(X̂, 1) < 12 && continue
        ref_result = matlab_ledoitwolf_analytical_shrinkage(X̂)
        @test cov(ANS, X̂) ≈ ref_result
    end

    # weights are not currently exported but it seems to be a good
    # idea to ensure proper errors
    fw1 = FrequencyWeights(rand(1:10, size(test_matrices[1], 1)))
    @test_throws ErrorException cov(ANS, test_matrices[1], fw1)
    @test_throws ErrorException cov(ANS, test_matrices[1], fw1; dims=2)
    @test_throws ErrorException cov(ANS, test_matrices[1], fw1; mean=nothing)
end
