function sum_var_sij(Z::AbstractMatrix, CovZ::AbstractMatrix, n::Int,
                     corrected=false; with_diag=true)
    scale = ifelse(corrected, n-1, n)
    Z²    = Z.^2
    π̂mat  = ((Z²'*Z²)/n - CovZ.^2) / scale
    with_diag && return sum(π̂mat)
    return sum(π̂mat) - sum(diag(π̂mat))
end
