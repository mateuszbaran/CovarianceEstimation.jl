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
    @test cov(X, ce; dims=1) ≈ cov(transpose(X), ce; dims=2)
    @test cov(X, ce; dims=2) ≈ cov(transpose(X), ce; dims=1)

    @test_throws ArgumentError cov(X, ce, dims=0)
    # XXX broken?
    # @test_throws ArgumentError cov(ce, X, dims=3)
end

function testUncorrelated(ce::CovarianceEstimator)
    for X2 ∈ X2s
        @test isdiag(cov(X2, ce))
    end
end

function testTranslation(ce::CovarianceEstimator)
    C1 = cov(X, ce)
    C2 = cov(X .+ randn(1, 8), ce)
    @test C1 ≈ C2
    C1t = cov(X', ce)
    C2t = cov(X' .+ randn(1, 3), ce)
    @test C1t ≈ C2t
end

@testset "Simple covariance                  " begin
    sc = Simple()
    @test cov(X, sc; dims=1) ≈ cov(X; dims=1, corrected = false)
    @test cov(X, sc; dims=2) ≈ cov(X; dims=2, corrected = false)
    @test cov(X[1,:], X[2,:], sc) ≈ cov(X[1,:], X[2,:]; corrected = false)
    @test cov(X[1,:], sc) ≈ cov(X[1,:]; corrected = false)
    testTransposition(sc)
    testUncorrelated(sc)
    testTranslation(sc)
end

@testset "Corrected covariance               " begin
    sc = Simple(corrected=true)
    @test cov(X, sc; dims=1) ≈ cov(X; dims=1, corrected = true)
    @test cov(X, sc; dims=2) ≈ cov(X; dims=2, corrected = true)
    @test cov(X[1,:], X[2,:], sc) ≈ cov(X[1,:], X[2,:]; corrected = true)
    @test cov(X[1,:], sc) ≈ cov(X[1,:]; corrected = true)
    testTransposition(sc)
    testUncorrelated(sc)
    testTranslation(sc)
end

@testset "Ledoit-Wolf covariance shrinkage   " begin
    lw = LinearShrinkageEstimator(ConstantCorrelation())
    testTransposition(lw)
    testUncorrelated(lw)
    testTranslation(lw)
    for X̂ ∈ test_matrices
        ref_results = matlab_ledoitwolf_covcor(X̂)
        @test cov(X̂, lw) ≈ ref_results["lwcov"]
    end
end

@testset "Chen covariance shrinkage          " begin
    rblw = LinearShrinkageEstimator(DiagonalCommonVariance(), :rblw)
    testTransposition(rblw)
    testUncorrelated(rblw)
    testTranslation(rblw)

    oas = LinearShrinkageEstimator(DiagonalCommonVariance(), :oas)
    testTransposition(oas)
    testUncorrelated(oas)
    testTranslation(oas)

    for X̂ ∈ test_matrices
        Ŝ_rblw = cov(X̂, rblw)
        Ŝ_oas  = cov(X̂, oas)

        CE.centercols!(X̂)
        n, p = size(X̂)
        Ŝ    = cov(X̂, Simple())

        F_ref = tr(Ŝ)/p * I
        # https://arxiv.org/pdf/0907.4698.pdf eq 17
        ρ_rblw_ref = ((n-2)/n*tr(Ŝ^2)+tr(Ŝ)^2)/((n+2)*(tr(Ŝ^2)-tr(Ŝ)^2/p))
        # https://arxiv.org/pdf/0907.4698.pdf eq 23
        ρ_oas_ref = ((1-2/p)*tr(Ŝ^2)+tr(Ŝ)^2)/((n+1-2/p)*(tr(Ŝ^2)-tr(Ŝ)^2/p))

        @test cov(X̂, rblw) ≈ CE.linshrink(Ŝ, F_ref, ρ_rblw_ref)
        #@test cov(X̂, oas) ≈ CE.linshrink(Ŝ, F_ref, ρ_oas_ref)
    end
end
