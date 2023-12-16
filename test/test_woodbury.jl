using CovarianceEstimation
using Test
using Statistics
using StatsBase
using WoodburyMatrices

function wapprox(W1::SymWoodbury, W2::SymWoodbury; atolU=1e-6, kwargs...)
    isapprox(W1.A, W2.A; kwargs...) && isapprox(W1.D, W2.D; kwargs...) || return false
    # For B, each column can have sign flipped without changing the matrix
    for (col1, col2) in zip(eachcol(W1.B), eachcol(W2.B))
        isapprox(col1, col2; atol=atolU, kwargs...) || isapprox(col1, -col2; atol=atolU, kwargs...) || return false
    end
    return true
end

function woodbury_test(ce::WoodburyEstimator, X̂)
    c = cov(ce, X̂)
    @test issymmetric(c)
    @test wapprox(cov(ce, X̂'; dims=2), c)
    n2 = size(X̂, 1) ÷ 2
    w = FrequencyWeights(vcat(ones(n2), zeros(size(X̂, 1) - n2)))
    chalf = cov(ce, X̂, w);
    @test wapprox(chalf, cov(ce, X̂[1:n2, :]))
    # Weight types besides FrequencyWeights are not supported
    aw1 = AnalyticWeights(rand(size(test_matrices[1], 1)))
    @test_throws Exception cov(ce, test_matrices[1], aw1)
    @test_throws Exception cov(ce, test_matrices[1], aw1; dims=2)
    @test_throws Exception cov(ce, test_matrices[1], aw1; mean=nothing)
end

@testset "Woodbury shrinkage" begin
    @test_throws ArgumentError NormLoss(:not_a_norm, 1)
    @test_throws ArgumentError NormLoss(:L2, 0)
    for norm in (:L1, :L2, :Linf)
        for pivotidx in 1:7
            if pivotidx ∈ (5, 7) || (norm == :L1 && pivotidx ∈ (3, 4))
                X̂ = test_matrices[argmax([size(X̂, 1) for X̂ ∈ test_matrices])]
                @test_throws ArgumentError cov(WoodburyEstimator(NormLoss(norm, pivotidx), 2), X̂)
                continue
            end
            loss = NormLoss(norm, pivotidx)
            for X̂ ∈ test_matrices
                size(X̂, 1) < 6 && continue
                woodbury_test(WoodburyEstimator(loss, 2), X̂)
            end
        end
    end
    @test_throws ArgumentError StatLoss(:oops)
    for mode in (:st, :ent, :div, :aff, :fre)
        loss = StatLoss(mode)
        for X̂ ∈ test_matrices
            size(X̂, 1) < 6 && continue
            woodbury_test(WoodburyEstimator(loss, 2), X̂)
        end
    end
    # Take a case where we have some idea of what to expect: make a low-rank modification to a noise matrix
    Y = Xe
    v = zeros(size(Y, 2))
    v[1] = 100
    Y = [Y; v'; -v']    # avoid perturbing the mean
    ce = WoodburyEstimator(StatLoss(:ent), 2; σ²=1)
    C = cov(ce, Y)
    @test Matrix(C) ≈ cov(Y) rtol=0.1
    @test C.D[1,1] ≈ 400 rtol=0.1
    @test size(C.D) == (1, 1) || abs(C.D[2,2]) < 1
    @test all(==(1), diag(C.A))
    @test abs(C.B[1,1]) > 0.99
    @test wapprox(cov(ce, similar(Y); UsV = CovarianceEstimation.tsvd(Y .- mean(Y; dims=1), 2)), C)
end
