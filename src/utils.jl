function linshrink(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    return Symmetric((one(λ) .- λ).*S .+ λ.*F)
end

function linshrink(F::UniformScaling, S::AbstractMatrix, λ::Real)
    return Symmetric((one(λ) .- λ).*S + λ.*F)
end

function linshrink!(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    F .= (one(λ) .- λ).*S .+ λ.*F
    return Symmetric(F)
end
