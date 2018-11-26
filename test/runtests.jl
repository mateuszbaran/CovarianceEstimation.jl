using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test

X = rand(3, 8)

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
    C = Z * transpose(Z)
    F, r̄ = CovarianceEstimation.ledoitwolfshrinkagetarget(C)
    @test F ≈ Matrix(6.0I, 3, 3)
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
