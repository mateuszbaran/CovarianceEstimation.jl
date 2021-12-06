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

@testset "BiweightMidcovariance" begin
    for c in (9, 5)
        for modify_sample_size in (true, false)
            @testset "basic properties c=$c, modify_sample_size=$modify_sample_size" begin
                ce = BiweightMidcovariance(; c=c, modify_sample_size=modify_sample_size)
                for X in test_matrices
                    testTransposition(ce, X)
                    testDims(ce, X)
                    testUncorrelated(ce)
                    testTranslation(ce, X)

                    # Test that other dispatches are consistent with the one tested above.
                    _test_cov_vec(ce, X)
                    _test_cov_vecvec(ce, X)
                end
            end
        end
    end

    @testset "Z should give zeros, not NaNs" begin
        # This covers a difference discovered against astropy.
        # It was due to not implementing the special case that covariance is defined to be 0
        # when MAD = 0.
        @test cov(BiweightMidcovariance(), Z) == zeros(Float64, 3, 3)
    end

    @testset "astropy reference (c=9.0)" begin
        # These tests are against references generated with astropy.
        #
        # References generated with:
        #   * System: Linux 5.4.0-90-generic, x86_64, Intel 8700K
        #   * Software: Python 3.7.6, astropy 4.3.1, numpy 1.21.2
        #
        # import os
        # import numpy
        # from astropy.stats import biweight_midcovariance
        # prefix = os.path.join(
        #     os.path.expanduser("~"), ".julia", "dev", "CovarianceEstimation", "test",
        #     "test_matrices"
        # )
        # for sample in ["20x100", "50x50", "100x20", "50x50tdist3"]:
        #     for modify_sample_size in [False, True]:
        #         path = os.path.join(prefix, f"{sample}.csv")
        #         X = numpy.genfromtxt(path).T
        #         c = biweight_midcovariance(X, modify_sample_size=modify_sample_size)
        #         mss = "_mss" if modify_sample_size else ""
        #         out_path = os.path.join(prefix, f"{sample}_bwcov{mss}.csv")
        #         with open(out_path, "w") as file_:
        #             for i in range(c.shape[0]):
        #                 file_.write(" ".join(str(el) for el in c[i, :]))
        #                 file_.write("\n")
        names = ["20x100", "50x50", "100x20", "50x50tdist3"]
        _test_refs(BiweightMidcovariance(), names, "bwcov")
        _test_refs(BiweightMidcovariance(; modify_sample_size=true), names, "bwcov_mss")
    end
end
