function _test_cov_vec(ce::CovarianceEstimator, X::AbstractMatrix{<:Real})
    c = cov(ce, X)
    @inbounds for i in 1:size(X, 2)
        @test cov(ce, @view(X[:, i])) ≈ c[i, i] atol = 1e-12 rtol = 1e-16
    end
end

function _test_cov_vecvec(ce::CovarianceEstimator, X::AbstractMatrix{<:Real})
    c = cov(ce, X)
    @inbounds for i in 1:size(X, 2)
        for j in 1:size(X, 2)
            cij = cov(ce, @view(X[:, i]), @view(X[:, j]))
            @test cij ≈ c[i, j] atol = 1e-12 rtol = 1e-16
        end
    end
end

@testset "BiweightMidcovariance: c=9.0, modify_sample_size=false" begin
    ce = BiweightMidcovariance()

    test_mat1 = _get_ref("20x100")
    test_mat2 = _get_ref("100x20")
    test_mat3 = _get_ref("50x50")

    testTransposition(ce, test_mat1)
    testUncorrelated(ce)
    testTranslation(ce, test_mat1)

    # Test that other dispatches are consistent with the main one we are testing.
    _test_cov_vec(ce, test_mat1)
    _test_cov_vecvec(ce, test_mat1)

    c1 = cov(ce, test_mat1)
    # TODO compare against references
end
