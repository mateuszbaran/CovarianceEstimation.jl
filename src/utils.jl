shrink(C::AbstractMatrix, F::Union{UniformScaling, AbstractMatrix}, ρ::Real) =
    (1.0 - ρ) * C + ρ * F

function centercols!(X::AbstractMatrix)
    # centering of the columns
    μ = mean(X, dims=1)
    X .-= μ
    return nothing
end
