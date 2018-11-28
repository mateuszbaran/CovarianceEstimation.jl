using LinearAlgebra

struct LedoitWolfCovariance <: CovarianceEstimator
end

function ledoitwolfshrinkagetarget(C::AbstractMatrix{<:Real})
    N = size(C, 1)
    Cs = [sqrt(C[i,i]*C[j,j]) for i in 1:N, j in 1:N]
    r = C ./ Cs
    r̄ = (sum(r)-N)/(N*(N-1))
    Finterm = Cs .* r̄
    F = Finterm + Diagonal(diag(Cs) .- diag(Finterm))
    F, r̄
end

"""
    cov(LedoitWolfCovariance(), X::AbstractMatrix, shrinkage; dims::Int=1)

Calculates shrunk covariance matrix for data in `X` with shrinkage
parameter `shrinkage` and Ledoit-Wolf shrinkage target.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,”
The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
"""
function cov(::LedoitWolfCovariance, X::AbstractMatrix{T}, shrinkage::Number; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end

    (0 ≤ shrinkage ≤ 1) || throw(ArgumentError("Shinkage must be in [0,1] (given shrinkage: $shrinkage)"))
    C = cov(Xint; dims=2)
    F, r̄ = ledoitwolfshrinkagetarget(C)
    (1-shrinkage)*C + shrinkage*F
end

"""
    cov(::LedoitWolfCovariance, X::AbstractMatrix; dims::Int=1)

Calculates shrunk covariance matrix for data in `X` with Ledoit-Wolf
optimal shrinkage.

# Arguments
- `dims::Int`: the dimension along which the variables are organized.
When `dims = 1`, the variables are considered columns with observations
in rows; when `dims = 2`, variables are in rows with observations in columns.

Implements shrinkage target and optimal shrinkage according to
O. Ledoit and M. Wolf, “Honey, I Shrunk the Sample Covariance Matrix,”
The Journal of Portfolio Management, vol. 30, no. 4, pp. 110–119, Jul. 2004.
"""
function cov(::LedoitWolfCovariance, X::AbstractMatrix{T}; dims::Int=1) where T<:Real
    if dims == 1
        Xint = transpose(X)
    elseif dims == 2
        Xint = X
    else
        throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))
    end
    Xint = Xint .- mean(Xint; dims=2)
    C = cov(Xint; dims=2)
    Tnum = size(Xint, 2)
    N = size(Xint, 1)
    F, r̄ = ledoitwolfshrinkagetarget(C)
    πmatrix = mean([(Xint[i,t]*Xint[j,t]-C[i,j])^2 for i in 1:N, j in 1:N] for t in 1:Tnum)
    πhat = sum(πmatrix)
    ϑhatii = mean([(Xint[i,t]^2 - C[i,i])*(Xint[i,t]*Xint[j,t] - C[i,j]) for i in 1:N, j in 1:N] for t in 1:Tnum)
    ϑhatjj = mean([(Xint[i,t]^2 - C[j,j])*(Xint[i,t]*Xint[j,t] - C[i,j]) for i in 1:N, j in 1:N] for t in 1:Tnum)
    ρhatpart2 = zero(T)
    #TODO: inbounds/simd?
    sdC = sqrt.(C[i,i] for i ∈ 1:N)
    for i in 1:N
        for j in 1:N
            αij = sdC[j]/sdC[i]
            ρhatpart2 += ϑhatii[i,j]*αij + ϑhatjj[i,j]/αij
        end
    end
    ρhat = sum(diag(πmatrix)) + (r̄/2)*ρhatpart2
    γhat = sum((F - C).^2)
    κhat = (πhat - ρhat)/γhat
    if abs(γhat) < 1e-16
        # division by almost zero, so shrinkage doesn't really matter
        δstar = 0.0
    else
        δstar = clamp(κhat/Tnum, 0.0, 1.0)
    end
    (1-δstar)*C + δstar*F
end
