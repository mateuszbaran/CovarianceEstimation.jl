function linshrink(F::Union{UniformScaling, AbstractMatrix},
                   S::AbstractMatrix, λ::Real)
    return Symmetric((1.0 - λ)*S + λ*F)
end


function linshrink!(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    F .*= λ
    F .+= (1.0 - λ)*S
    return Symmetric(F)
end
