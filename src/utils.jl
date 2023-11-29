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

totalweight(n) = n
totalweight(_, weights) = sum(weights)

# Dividing by zero produces zero
guardeddiv(num, denom) = iszero(denom) ? zero(num)/oneunit(denom) : num/denom
diaginv(guard::Bool, num, v) = guard ? map(z -> guardeddiv(num, z), v) : num ./ v
