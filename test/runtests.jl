using CovarianceEstimation
using Statistics
using LinearAlgebra
using Test
using Random
using DelimitedFiles
using StatsBase


include("ref_lw_lshrink.jl")
include("ref_lw_nlshrink.jl")
include("legacy.jl")
Random.seed!(1234)


const CE = CovarianceEstimation

const X   = randn(3, 8)
const Z   = [2 -1 2; -1 2 -1; -1 -1 -1]
const X2s = [[1 0; 0 1; -1 0; 0 -1], [1 0; 0 1; 0 -1; -1 0]]
const Xa  = randn(5, 5)
const Xb  = randn(8, 3)
const Xc  = randn(15, 20)
const Xd  = randn(30, 50)
const Xe  = randn(50, 30)

const test_matrices = [X, Xa, Xb, Xc, Xd, Xe, Z]


function testTransposition(ce::CovarianceEstimator, X)
    @test cov(ce, X; dims=1) ≈ cov(ce, transpose(X); dims=2)
    @test cov(ce, X; dims=2) ≈ cov(ce, transpose(X); dims=1)
end

function testDims(ce::CovarianceEstimator, X)
    @test_throws ArgumentError cov(ce, X, dims=0)
    @test_throws ArgumentError cov(ce, X, dims=3)
end

function testUncorrelated(ce::CovarianceEstimator)
    for X2 ∈ X2s
        @test isdiag(cov(ce, X2))
    end

end

function testTranslation(ce::CovarianceEstimator, X)
    C1 = cov(ce, X)
    C2 = cov(ce, X .+ randn(1, size(X, 2)))
    @test C1 ≈ C2 atol = 1e-12 rtol = 1e-16
    C1t = cov(ce, X')
    C2t = cov(ce, X' .+ randn(1, size(X, 1)))
    @test C1t ≈ C2t atol = 1e-12 rtol = 1e-16
end

include("test_linearshrinkage.jl")
include("test_nonlinearshrinkage.jl")
