using LinearAlgebra

struct LedoitWolfCovariance <: CovarianceEstimator
end

function ledoitwolfshrinkagetarget(C::DenseMatrix{<:Real})
    N = size(C, 1)
    Cs = [sqrt(C[i,i]*C[j,j]) for i in 1:N, j in 1:N]
    r = C ./ Cs
    r̄ = (sum(r)-N)/(N*(N-1))
    Finterm = Cs .* r̄
    F = Finterm - Diagonal(diag(Finterm)) + Diagonal(diag(Cs))
    F, r̄
end

"""
    ldacalccov(LedoitWolfCovariance(), X, shrinkage)

Calculates shrunk covariance matrix for centered data `X` with shrinkage
parameter `shrinkage` and Ledoit-Wolf shrinkage target.

Implements shrinkage target and optimal shrinkage according to
O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,”
The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
"""
function cov(::LedoitWolfCovariance, X::DenseMatrix{T}, shrinkage::Number) where T<:Real
    (shrinkage ≥ 0 && shrinkage ≤ 1) || throw(ArgumentError("Shinkage must be in [0,1] (given shrinkage: $shrinkage)"))
    C = X * transpose(X)
    F, r̄ = ledoitwolfshrinkagetarget(C)
    (1-shrinkage)*C + shrinkage*F
end

"""
    ldacalccov(::LedoitWolfCovariance, X::DenseMatrix)

Calculates shrunk covariance matrix for centered data `X` with
Ledoit-Wolf optimal shrinkage.

Implements shrinkage target and optimal shrinkage according to
O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,”
The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
"""
function cov(::LedoitWolfCovariance, X::DenseMatrix{T}) where T<:Real
    C = X * transpose(X)
    Tnum = size(X, 2)
    N = size(X, 1)
    F, r̄ = ledoitwolfshrinkagetarget(C)
    πmatrix = mean([(X[i,t]*X[j,t]-C[i,j])^2 for i in 1:N, j in 1:N] for t in 1:Tnum)
    πhat = sum(πmatrix)
    ϑhatii = mean([(X[i,t]^2 - C[i,i])*(X[i,t]*X[j,t] - C[i,j]) for i in 1:N, j in 1:N] for t in 1:Tnum)
    ϑhatjj = mean([(X[i,t]^2 - C[j,j])*(X[i,t]*X[j,t] - C[i,j]) for i in 1:N, j in 1:N] for t in 1:Tnum)
    ρhatpart2 = zero(T)
    #TODO: inbounds/simd?
    for i in 1:N
        for j in 1:N
            ρhatpart2 += sqrt(C[j,j]/C[i,i])*ϑhatii[i,j] + sqrt(C[i,i]/C[j,j])*ϑhatjj[i,j]
        end
    end
    ρhat = sum(diag(πmatrix)) + (r̄/2)*ρhatpart2
    γhat = sum((F - C).^2)
    κhat = (πhat - ρhat)/γhat
    δstar = clamp(κhat/Tnum, 0.0, 1.0)
    (1-δstar)*C + δstar*F
end
