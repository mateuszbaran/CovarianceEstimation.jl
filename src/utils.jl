function centercols!(X::AbstractMatrix)
    # centering of the columns
    μ = mean(X, dims=1)
    X .-= μ
    return nothing
end


safe01clamp(x::Real) = isnan(x) ? 0.0 : clamp(x, 0.0, 1.0)


function linshrink(S::AbstractMatrix, F::Union{UniformScaling, AbstractMatrix},
                   ρ::Real)
    ρsafe = safe01clamp(ρ)
    return (1.0 - ρsafe) * S + ρsafe * F
end
