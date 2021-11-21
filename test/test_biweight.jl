@testset "BiweightMidcovariance: c=9.0, modify_sample_size=false" begin
    ce = BiweightMidcovariance()

    test_mat1 = _get_ref("20x100")
    test_mat2 = _get_ref("100x20")
    test_mat3 = _get_ref("50x50")

    testTransposition(ce, test_mat1)
    testUncorrelated(ce)
    testTranslation(ce, test_mat1)

    c1 = cov(ce, test_mat1)
    # TODO compare against references
end
