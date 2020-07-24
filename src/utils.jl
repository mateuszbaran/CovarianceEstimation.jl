function linshrink(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    return Symmetric((1.0 .- λ).*S .+ λ.*F)
end

function linshrink(F::UniformScaling, S::AbstractMatrix, λ::Real)
    return Symmetric((1.0 .- λ).*S + λ.*F)
end

function linshrink!(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    F .= (1.0 .- λ).*S .+ λ.*F
    return Symmetric(F)
end
