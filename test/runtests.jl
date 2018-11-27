using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test
using Random

Random.seed!(1234)

X = randn(3, 8)

function testTransposition(ce::CovarianceEstimator)
    @test cov(ce, X; dims=1) ≈ cov(ce, transpose(X); dims=2)
    @test cov(ce, X; dims=2) ≈ cov(ce, transpose(X); dims=1)

    @test_throws ArgumentError cov(ce, X, dims=0)
    # broken?
    # @test_throws ArgumentError cov(ce, X, dims=3)
end

function testUncorrelated(ce::CovarianceEstimator; corrected = false)
    X2s = [[1 0; 0 1; -1 0; 0 -1],
           [2 0; 0 2; -1 0; 0 -1]]

    for X2 ∈ X2s
        @test cov(ce, X2) ≈ cov(X2, corrected = corrected)
    end
end

@testset "Simple covariance" begin
    sc = SimpleCovariance()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = false)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = false)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = false)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = false)
    testTransposition(sc)
    testUncorrelated(sc; corrected = false)
end

@testset "Corrected covariance" begin
    sc = CorrectedCovariance()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = true)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = true)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = true)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = true)
    testTransposition(sc)
    testUncorrelated(sc; corrected = true)
end

@testset "Ledoit-Wolf covariance shrinkage" begin
    lwc = LedoitWolfCovariance()
    Z = [2. -1 -1; -1 2 -1; 2 -1 -1]
    C = cov(Z; dims=2)
    F, r̄ = CovarianceEstimation.ledoitwolfshrinkagetarget(C)
    @test F ≈ Matrix(3.0I, 3, 3)
    @test r̄ ≈ 0.0
    (N, Tnum) = size(Z)
    πmatrix = mean([(Z[i,t]*Z[j,t]-C[i,j])^2 for i in 1:N, j in 1:N] for t in 1:Tnum)
    # the following assumes that r̄ ≈ 0.0
    πhat = sum(πmatrix)
    ρhat = sum(diag(πmatrix))
    γhat = sum((F - C).^2)
    κhat = (πhat - ρhat)/γhat
    δstar = clamp(κhat/Tnum, 0.0, 1.0)
    shrunkcov = (1-δstar)*C + δstar*F
    @test cov(lwc, Z; dims=2) ≈ shrunkcov
    @test cov(lwc, Z, δstar; dims=2) ≈ shrunkcov

    testTransposition(lwc)
    testUncorrelated(lwc; corrected = true)
end

@testset "Chen covariance shrinkage" begin
    Z = [2. -1 -1; -1 2 -1; 2 -1 -1]
    C = cov(Z; dims=2)
    F = CovarianceEstimation.chenshrinkagetarget(Z)
    @test F ≈ Matrix(3.0I, 3, 3)
    ρhatZrblw = 99/135
    shrunkcov = (1-ρhatZrblw)*C + ρhatZrblw*F
    rblwc = RBLWCovariance()
    @test cov(rblwc, Z, ρhatZrblw; dims=2) ≈ shrunkcov
    @test cov(rblwc, Z; dims=2) ≈ shrunkcov
    testTransposition(rblwc)
    # broken
    # testUncorrelated(rblwc; corrected = true)

    ρhatZoas = 1
    shrunkcov = (1-ρhatZoas)*C + ρhatZoas*F
    oasc = OASCovariance()
    @test cov(oasc, Z, ρhatZoas; dims=2) ≈ shrunkcov
    @test cov(oasc, Z; dims=2) ≈ shrunkcov
    testTransposition(oasc)
    # broken
    # testUncorrelated(oasc; corrected = true)
end
