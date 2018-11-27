using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test
using Random

Random.seed!(1234)

X = randn(3, 8)

@testset "Simple covariance" begin
    sc = SimpleCovariance()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = false)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = false)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = false)
end

@testset "Corrected covariance" begin
    sc = CorrectedCovariance()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = true)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = true)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = true)
end

@testset "Ledoit-Wolf covariance shrinkage" begin
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
    @test cov(LedoitWolfCovariance(), Z) ≈ shrunkcov
    @test cov(LedoitWolfCovariance(), Z, δstar) ≈ shrunkcov
end

@testset "Chen covariance shrinkage" begin
    Z = [2. -1 -1; -1 2 -1; 2 -1 -1]
    C = cov(Z; dims=2)
    F = CovarianceEstimation.chenshrinkagetarget(Z)
    @test F ≈ Matrix(3.0I, 3, 3)
    ρhatZrblw = 99/135
    shrunkcov = (1-ρhatZrblw)*C + ρhatZrblw*F
    @test cov(RBLWCovariance(), Z, ρhatZrblw) ≈ shrunkcov
    @test cov(RBLWCovariance(), Z) ≈ shrunkcov

    ρhatZoas = 1
    shrunkcov = (1-ρhatZoas)*C + ρhatZoas*F
    @test cov(OASCovariance(), Z, ρhatZoas) ≈ shrunkcov
    @test cov(OASCovariance(), Z) ≈ shrunkcov
end
