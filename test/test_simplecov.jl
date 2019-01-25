@testset "Simple covariance                  " begin
    sc = Simple()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = false)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = false)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = false)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = false)
    testTransposition(sc, X)
    testUncorrelated(sc)
    testTranslation(sc, X)

    sc = Simple(corrected=true)
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = true)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = true)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = true)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = true)
    testTransposition(sc, X)
    testUncorrelated(sc)
    testTranslation(sc, X)

    fw1 = FrequencyWeights(rand(1:10, size(X, 1)))
    sfwc = Simple(fw1)
    @test cov(sfwc, X) ≈ cov(X, fw1; corrected = sfwc.corrected)

    fw2 = FrequencyWeights(rand(1:10, size(X, 2)))
    sfwc = Simple(fw2)
    @test cov(sfwc, X, dims=2) ≈ cov(X, fw2, 2; corrected = sfwc.corrected)

    sfwc = Simple(fw1, corrected = true)
    @test cov(sfwc, X) ≈ cov(X, fw1, corrected = sfwc.corrected)
end
