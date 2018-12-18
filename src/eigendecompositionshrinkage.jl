
struct AnalyticalNonlinearShrinkage <: CovarianceEstimator end

"""
    analytical_nonlinear_shrinkage(X)

Based on Matlab code in Olivier Ledoit and Michael Wolf. Analytical Nonlinear
Shrinkage of Large-Dimensional Covariance Matrices. (Nov 2018)
http://www.econ.uzh.ch/static/wp/econwp264.pdf
"""
function analytical_nonlinear_shrinkage(X::AbstractMatrix; decomp::Union{Eigen,Nothing} = nothing)
    n, p = size(X)
    if n < 12
        # explained in the paper
        throw(ArgumentError("The number of samples `n` must be at least 12 (given: $n)."))
    end
    sample = (X'*X)./n
    # sample eigenvalues sorted in ascending order and eigenvectors
    F = isa(decomp, Nothing) ? eigen(sample) : decomp
    perm = sortperm(F.values)
    λ    = F.values[perm]
    u    = F.vectors[:,perm]
    # compute analytical nonlinear shrinkage kernel formula
    λ    = λ[max(1,p-n+1):p]
    L    = repeat(λ, outer = (1, min(p, n)))
    # Equation (4.9)
    h    = n^(-1/3)
    H    = h*L'
    x    = (L-L')./H
    f̃    = (3/4/sqrt(5))*mean(max.(1 .-(x.^2)./5, 0)./H, dims=2)
    # Equation (4.7)
    Hftemp = (-3/10/π)*x+(3/4/sqrt(5)/π)*(1 .-(x.^2)./5).*log.(abs.((sqrt(5).-x)./(sqrt(5).+x)))
    # Equation (4.8)
    Hftemp[abs.(x).==sqrt(5)] .= (-3/10/π)*x[abs.(x).==sqrt(5)]
    Hf̃ = mean(Hftemp./H, dims=2)
    if p<=n
        # Equation (4.3)
        d̃ = λ./((π*(p/n)*λ.*f̃).^2 + (1-(p/n).-π*(p/n)*λ.*Hf̃).^2)
    else
        # Equation (C.8)
        Hf̃0 = (1/π)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2)*log((1+sqrt(5)*h)/(1-sqrt(5)*h)))*mean(1 ./λ)
        d̃0 = 1/(π*(p-n)/n*Hf̃0)
        # Equation (C.5)
        d̃1=λ./(π^2*(λ.^2).*(f̃.^2+Hf̃.^2))
        # Eq. (C.4)
        d̃ = vcat(d̃0*ones(p-n,1), d̃1)
    end
    # Equation (4.4)
    return u*Diagonal(vec(d̃'))*u'
end

function cov(X::AbstractMatrix{<:Real}, ::AnalyticalNonlinearShrinkage;
             dims::Int=1, decomp::Union{Eigen,Nothing}=nothing)

    @assert dims ∈ [1, 2] "Argument dims can only be 1 or 2 (given: $dims)"

    Xc = (dims == 1) ? centercols(X) : centercols(transpose(X))
    return analytical_nonlinear_shrinkage(Xc, decomp = decomp)
end
