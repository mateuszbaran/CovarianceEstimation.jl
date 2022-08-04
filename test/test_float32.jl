@testset "Float32 matrices" begin
    # linear shrinkages
    for X in test_matrices
        x = convert(Matrix{Float32}, X)
        for target in [
            DiagonalUnitVariance(),
            DiagonalCommonVariance(),
            DiagonalUnequalVariance(),
            CommonCovariance(),
            PerfectPositiveCorrelation(),
            ConstantCorrelation()
            ]
            for shrinkage in [:lw, :ss]
                @test eltype(cov(LinearShrinkage(target, shrinkage), x)) == Float32
            end
        end
    end

    # nonlinear shrinkages
    ANS = AnalyticalNonlinearShrinkage()
    for X in test_matrices
        size(X, 1) < 12 && continue
        x = convert(Matrix{Float32}, X)
        @test eltype(cov(ANS, x)) == Float32
    end
end
