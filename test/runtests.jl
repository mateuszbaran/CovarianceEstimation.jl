using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test
using Random

include("reference_ledoitwolf.jl")

Random.seed!(1234)

const X = randn(3, 8)
const Z = [2 -1 2; -1 2 -1; -1 -1 -1]
const X2s = [[1 0; 0 1; -1 0; 0 -1], [1 0; 0 1; 0 -1; -1 0]]

function testTransposition(ce::CovarianceEstimator)
    @test cov(ce, X; dims=1) ≈ cov(ce, transpose(X); dims=2)
    @test cov(ce, X; dims=2) ≈ cov(ce, transpose(X); dims=1)

    @test_throws ArgumentError cov(ce, X, dims=0)
    # XXX broken?
    # @test_throws ArgumentError cov(ce, X, dims=3)
end

function testUncorrelated(ce::CovarianceEstimator)
    for X2 ∈ X2s
        @test isdiag(cov(ce, X2))
    end
end

function testTranslation(ce::CovarianceEstimator)
    C1 = cov(ce, X)
    C2 = cov(ce, X .+ randn(1, 8))
    @test C1 ≈ C2
    C1t = cov(ce, X')
    C2t = cov(ce, X' .+ randn(1, 3))
    @test C1t ≈ C2t
end

@testset "Simple covariance                  " begin
    sc = Simple()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = false)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = false)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = false)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = false)
    testTransposition(sc)
    testUncorrelated(sc)
    testTranslation(sc)
end

@testset "Corrected covariance               " begin
    sc = Corrected()
    @test cov(sc, X; dims=1) ≈ cov(X; dims=1, corrected = true)
    @test cov(sc, X; dims=2) ≈ cov(X; dims=2, corrected = true)
    @test cov(sc, X[1,:], X[2,:]) ≈ cov(X[1,:], X[2,:]; corrected = true)
    @test cov(sc, X[1,:]) ≈ cov(X[1,:]; corrected = true)
    testTransposition(sc)
    testUncorrelated(sc)
    testTranslation(sc)
end

@testset "Ledoit-Wolf covariance shrinkage   " begin
    lw = LedoitWolf()
    testTransposition(lw)
    testUncorrelated(lw)
    testTranslation(lw)
    for X̂ ∈ [X, Z]
        ref_results = matlab_ledoitwolf_covcor(X̂)
        C = cov(Simple(), X̂; dims=1)
        sdC = sqrt.(diag(C))
        F, r̄ = CovarianceEstimation.lw_shrinkagetarget(C, sdC)
        @test r̄ ≈ ref_results["r̄"]
        @test F ≈ ref_results["F"]
        shrinkage = CovarianceEstimation.lw_optimalshrinkage(X̂, C, sdC, F, r̄)
        @test shrinkage ≈ ref_results["shrinkage"]
        @test cov(lw, X̂) ≈ ref_results["lwcov"]
    end
end

@testset "Chen covariance shrinkage          " begin
    Z = [2. -1 -1; -1 2 -1; 2 -1 -1]
    C = cov(Z; dims=2)
    F = CovarianceEstimation.rblw_shrinkagetarget(Z)
    @test F ≈ Matrix(3.0I, 3, 3)
    ρhatZrblw = 99/135
    shrunkcov = (1-ρhatZrblw)*C + ρhatZrblw*F
    rblwc = RaoBlackwellLedoitWolf()
    rblwcOptim = RaoBlackwellLedoitWolf(ρhatZrblw)
    @test cov(rblwcOptim, Z; dims=2) ≈ shrunkcov
    @test cov(rblwc, Z; dims=2) ≈ shrunkcov
    testTransposition(rblwc)
    testUncorrelated(rblwc)
    testTranslation(rblwc)

    testTransposition(rblwcOptim)

    ρhatZoas = 0.6
    Cdim1 = cov(Z; dims=1)
    Fdim1 = CovarianceEstimation.rblw_shrinkagetarget(Z, dims=1)
    shrunkcov = (1-ρhatZoas)*Cdim1 + ρhatZoas*Fdim1
    oasc = OracleApproximatingShrinkage()
    oascOptim = OracleApproximatingShrinkage(ρhatZoas)
    @test cov(oascOptim, Z; dims=1) ≈ shrunkcov
    @test cov(oasc, Z; dims=1) ≈ shrunkcov
    testTransposition(oasc)
    testUncorrelated(oasc)
    testTranslation(oasc)

    testTransposition(oascOptim)
end
