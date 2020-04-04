function linshrink(F::Union{UniformScaling, AbstractMatrix},
                   S::AbstractMatrix, λ::Real)
    return Symmetric((1.0 - λ)*S + λ*F)
end

function linshrink!(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    F .*= λ
    F .+= (1.0 - λ)*S
    return Symmetric(F)
end

function subtract_mean!(Xc::AbstractMatrix, mean; dims::Int=1)
    if mean === nothing
        Xc .-= Statistics.mean(Xc, dims=1)
    elseif mean isa AbstractArray
        if dims == 1
            Xc .-= mean
        else
            Xc .-= mean'
        end
    end
end
