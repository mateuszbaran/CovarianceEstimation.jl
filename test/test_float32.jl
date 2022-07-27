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
end
