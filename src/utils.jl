const AM = AbstractMatrix

shrink(C::AM, F::Union{UniformScaling, AM}, ρ::Real) = (1.0 - ρ)*C + ρ*F

function centercols!(X::AM)
    # centering of the columns
    n, p = size(X)
    μ = mean(X, dims=1)
    @inbounds for i ∈ 1:n, j ∈ 1:p
        X[i, j] -= μ[j]
    end
    return nothing
end
