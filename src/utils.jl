"""
    centercols(X, μ; dims)

Internal function to center the columns of `X` given `μ` which can either be
    * `nothing` (in which case the mean is estimated and subtracted from `X`)
    * `0` (in which case nothing is done, the matrix is assumed to be centred)
    * a given vector (in which case that vector is subtracted from `X`)
"""
function centercols(X::AbstractMatrix, μ=nothing; dims=1)
    Xc = (dims == 1) ? copy(X) : transpose(X)
    if μ === nothing
        Xc .-= mean(X, dims=1)
    elseif iszero(μ)
        # assume X is centered already
    elseif μ isa AbstractVector
        @assert length(μ) == size(Xc, 2) "Dimension mismatch, given mean has length $(length(μ)) while p=$(size(Xc, 2))"
        Xc .-= μ
    else
        throw(ArgumentError("Incorrect assignment for `mean`, expected either 0, nothing or a vector, got $μ."))
    end
    Xc
end

function linshrink(F::Union{UniformScaling, AbstractMatrix},
                   S::AbstractMatrix, λ::Real)
    return (1.0 - λ)*S + λ*F
end


function linshrink!(F::AbstractMatrix, S::AbstractMatrix, λ::Real)
    F .*= λ
    F .+= (1.0 - λ) * S
    return F
end
