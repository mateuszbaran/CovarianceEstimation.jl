centercols(X::AbstractMatrix) = (X .- mean(X, dims=1))


linshrink(S::AbstractMatrix, F::Union{UniformScaling, AbstractMatrix},
          λ::Real) = (1.0 - λ) * S + λ * F
