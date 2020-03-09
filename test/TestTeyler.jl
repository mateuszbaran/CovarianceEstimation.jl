using LinearAlgebra


## Tyler M-Estimator fixed point algorithm (Tyler, 1987)
# `X` (the data) must be a wide matrix (for the sake of efficiency)
# `tol` is the stopping criterion
# `maxiter` is the maximum number of iterations allowed
# if `verbose`, information on convergence will be printed in the REPL.
function tme(   X::AbstractMatrix{T};
                tol::Real = real(T)(1e-6),
                maxiter::Int = 200,
                verbose::Bool = false) where {T<:Union{Real,Complex}}
    n, t = size(X)
    R = Matrix{T}(I, n, n)
    Rnew = Matrix{T}(undef, n, n)
    iter, 😋 = 1, false

    verbose && println("Iterating M-estimator fixed-point algorithm...")
    while true
        C = cholesky(R)
        fill!(Rnew, zero(T))
        @inbounds for i = 1:t
            @views v = C.L \ X[:, i]
            Rnew += (X[:, i] .* X[:, i]') ./ (v ⋅ v)
        end
        Rnew *= inv(tr(Rnew))
        conv = norm(Rnew - R) / norm(R)
        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (overRun = iter == maxiter) && @warn(
            "M-estimator reached the max number of iterations before convergence:",
            iter,
        )
        (😋 = conv <= tol) || overRun == true ? break : (iter += 1; R[:] = Rnew)
    end # while
    verbose && @info("Convergence has " * (😋 ? "" : "not ") * "been attained.\n\n")
    return Rnew
end

using MAT, Test

# change this to point to the directory where the .mat files are
dataDir="C:\\Users\\congedom\\Documents\\Code\\MATLAB\\Tayler"

files=["file3_100.mat", "file5_100.mat", "file10_100.mat"]
M=Vector{Matrix}(undef, 3)
myM=similar(M)

for (i, f) ∈ enumerate(files)
    file = joinpath(dataDir, f) # the \\ works only in
    vars = matread(file)
    X    = vars["X"] # data matrix
    M[i] = vars["M"] # benchmark estimator
    myM[i]=tme(X) # CovarianceEstimation.jl estimator
end

@test M[1]≈myM[1]
@test M[2]≈myM[2]
@test M[3]≈myM[3]
