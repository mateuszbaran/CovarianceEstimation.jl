@testset "Analytical Nonlin shrinkage (ref⭒) " begin
    ANS = AnalyticalNonlinearShrinkage()
    for X̂ ∈ test_matrices
        size(X̂, 1) < 12 && continue
        ref_result = matlab_ledoitwolf_analytical_shrinkage(X̂)
        @test cov(X̂, ANS) ≈ ref_result
    end
end
