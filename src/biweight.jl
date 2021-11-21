"""
    BiweightMidcovariance(c=9.0)

The biweight midcovariance is a covariance estimator that is resilient to outliers.

Consider two random variables ``x`` and ``y``, for which we have ``n`` observations
``\\{(x_i, y_i)\\}``. The biweight midcovariance is then defined to be:
```math
n\\cdot\\frac{
    \\sum_{|u_i|<1,|v_i|<1}(x_i - M_x)(1 - u_i^2)^2(y_i - M_y)(1 - v_i^2)^2
}{
    \\left(\\sum_{|u_i|<1}(1 - u_i^2)(1-5u_i^2)\\right)
    \\left(\\sum_{|v_i|<1}(1 - v_i^2)(1-5v_i^2)\\right)
}
```
where ``M_x`` and ``M_y`` are the medians of ``\\{x_i\\}`` and ``\\{y_i\\}`` respectively,
and
```math
\\begin{aligned}
u_i &= \\frac{x_i - M_x}{c \\cdot \\mathrm{MAD}_x} \\\\
v_i &= \\frac{y_i - M_y}{c \\cdot \\mathrm{MAD}_y},
\\end{aligned}
```
where ``\\mathrm{MAD}`` represents the median absolute deviation,
```math
\\begin{aligned}
\\mathrm{MAD}_x &= \\mathrm{median}(\\left\\{\\left|x_i - M_x\\right|\\right\\}) \\\\
\\mathrm{MAD}_y &= \\mathrm{median}(\\left\\{\\left|y_i - M_y\\right|\\right\\}).
\\end{aligned}
```

The parameter ``c`` is a tuning constant, for which the default is ``9.0``.
Larger values will reduce the number of outliers that are removed — i.e. reducing
robustness, but increasing sample efficiency.
"""
struct BiweightMidcovariance <: CovarianceEstimator
    c::Float64
    modify_sample_size::Bool
end
function BiweightMidcovariance(; c::Float64=9.0, modify_sample_size::Bool=false)
    return BiweightMidcovariance(c, modify_sample_size)
end

function covzm(ce::BiweightMidcovariance, x::AbstractVector{<:Real})
    MADx = median(abs.(x))
    numerator = zero(eltype(x))
    denominator = zero(eltype(x))
    count = 0
    for xi in x
        ui² = (xi / (ce.c * MADx))^2
        ui² >= 1 && continue
        count += 1
        numerator += xi^2 * (1 - ui²)^4
        denominator += (1 - ui²) * (1 - 5 * ui²)
    end
    n = ce.modify_sample_size ? count : length(x)
    return n * numerator / (denominator^2)
end

function cov(
    ce::BiweightMidcovariance,
    x::AbstractVector{<:Real};
    mean::Union{Nothing,<:Real}=nothing,
)
    Mx = isnothing(mean) ? median(x) : mean
    return covzm(ce, x .- Mx)
end

function covzm(
    ce::BiweightMidcovariance, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}
)
    MADx = median(abs.(x))
    MADy = median(abs.(y))

    # Promote types between x & y for numerator
    numerator = zero(promote_type(eltype(x), eltype(y)))
    denominator_x = zero(eltype(x))
    denominator_y = zero(eltype(y))
    count = 0
    for (xi, yi) in zip(x, y)
        ui² = (xi / (ce.c * MADx))^2
        vi² = (yi / (ce.c * MADy))^2

        if ui² < 1
            denominator_x += (1 - ui²) * (1 - 5 * ui²)
        end

        if vi² < 1
            denominator_y += (1 - vi²) * (1 - 5 * vi²)
        end

        if ui² < 1 && vi² < 1
            numerator += xi * (1 - ui²)^2 * yi * (1 - vi²)^2
            count += 1
        end
    end

    n = ce.modify_sample_size ? count : length(x)
    return n * numerator / (denominator_x * denominator_y)
end

function cov(
    ce::BiweightMidcovariance, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}
)
    return covzm(ce, x .- median(x), y .- median(y))
end

function cov(
    ce::BiweightMidcovariance,
    X::AbstractMatrix{<:Real};
    mean::Union{Nothing,AbstractVector{<:Real},AbstractMatrix{<:Real}}=nothing,
    dims::Int=1,
)
    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))

    # Deal with the `mean` argument -- if it isn't provided, then we centralise with the
    # median.
    MX = isnothing(mean) ? median(X; dims=dims) : mean
    X = X .- MX

    # Standardise the orientation of X, and ascertain the number of variables.
    X = dims == 2 ? transpose(X) : X
    p = size(X, 2)

    # Infer the output element type. This follows the approach used in the stdlib, in
    # Statistics.covzm.
    T = promote_type(typeof(one(eltype(X)) / 1), eltype(X))

    # Allocate a concrete output buffer, and populate it by iterating over all pairs.
    # Since different outliers may be removed
    result = Matrix{T}(undef, p, p)
    @inbounds for i in 1:p
        for j in i:p
            if i == j
                result[i, j] = covzm(ce, @view(X[:, i]))
            else
                cov_ij = covzm(ce, @view(X[:, i]), @view(X[:, j]))
                result[i, j] = cov_ij
                result[j, i] = cov_ij
            end
        end
    end

    return result
end
