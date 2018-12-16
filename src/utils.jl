function centercols!(X::AbstractMatrix)
    # centering of the columns
    μ = mean(X, dims=1)
    X .-= μ
    return nothing
end


linshrink(S::AbstractMatrix, F::Union{UniformScaling, AbstractMatrix},
          λ::Real) = (1.0 - λ) * S + λ * F
