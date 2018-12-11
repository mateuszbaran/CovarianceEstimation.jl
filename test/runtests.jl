using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test
using Random

const CE = CovarianceEstimation

include("reference_ledoitwolf.jl")

Random.seed!(1234)

const X = randn(3, 8)
const Z = [2 -1 2; -1 2 -1; -1 -1 -1]
const X2s = [[1 0; 0 1; -1 0; 0 -1], [1 0; 0 1; 0 -1; -1 0]]

const test_matrices = [X, Z]

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
    for X̂ ∈ test_matrices
        ref_results = matlab_ledoitwolf_covcor(X̂)
        # center columns
        X̂c = X̂
        n, p = size(X̂c)
        μ = mean(X̂c, dims=1)
        @inbounds for i ∈ 1:n, j ∈ 1:p
            X̂c[i, j] -= μ[j]
        end
        # compute the different elements and check they match the reference
        Ŝ    = cov(Simple(), X̂)
        V̂    = sqrt.(diag(Ŝ))
        F, r̄ = CE.lw_shrinkagetarget(Ŝ, V̂, p)
        @test r̄ ≈ ref_results["r̄"]
        @test F ≈ ref_results["F"]
        shrinkage = CE.lw_optimalshrinkage(X̂c, Ŝ, V̂, F, r̄, n, p)
        @test shrinkage ≈ ref_results["shrinkage"]
        @test cov(lw, X̂) ≈ ref_results["lwcov"]
    end
end

@testset "Chen covariance shrinkage          " begin
    rblw = RaoBlackwellLedoitWolf()
    testTransposition(rblw)
    testUncorrelated(rblw)
    testTranslation(rblw)

    oas = OracleApproximatingShrinkage()
    testTransposition(oas)
    testUncorrelated(oas)
    testTranslation(oas)

    for X̂ ∈ test_matrices
        Ŝ_rblw = cov(rblw, X̂)
        Ŝ_oas  = cov(oas, X̂)

        CE.centercols!(X̂)
        n, p = size(X̂)
        Ŝ    = cov(Simple(), X̂)

        F = CE.rblw_shrinkagetarget(Ŝ, p)

        ρ_rblw = CE.rblw_optimalshrinkage(Ŝ, n, p)
        ρ_oas  = CE.oas_optimalshrinkage(Ŝ, n, p)

        @test F ≈ tr(Ŝ)/p * I
        # https://arxiv.org/pdf/0907.4698.pdf eq 17
        ρ_rblw_ref = ((n-2)/n*tr(Ŝ^2)+tr(Ŝ)^2)/((n+2)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        ρ_rblw_ref = min(ρ_rblw_ref, 1)
        @test ρ_rblw ≈ ρ_rblw_ref
        # https://arxiv.org/pdf/0907.4698.pdf eq 23
        ρ_oas_ref = ((1-2/p)*tr(Ŝ^2)+tr(Ŝ)^2)/((n+1-2/p)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        ρ_oas_ref = min(ρ_oas_ref, 1)
        @test ρ_oas ≈ ρ_oas_ref

        @test Ŝ_rblw ≈ CE.shrink(Ŝ, F, ρ_rblw_ref)
        @test Ŝ_oas ≈ CE.shrink(Ŝ, F, ρ_oas_ref)
    end
end
