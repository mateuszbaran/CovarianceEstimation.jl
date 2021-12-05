"""
    BiweightMidcovariance(; c=9.0, modify_sample_size=false)

The biweight midcovariance is a covariance estimator that is resilient to outliers.

The technique derives originally from astrophysics [1], and is implemented in the Python
module [Astropy](https://www.astropy.org/) [2], as well as in NIST's Dataplot [3].

Consider two random variables ``x`` and ``y``, for which we have ``n`` observations
``\\{(x_i, y_i)\\}``. The biweight midcovariance is then defined to be:
```math
n_s\\cdot\\frac{
    \\sum_{|u_i|<1,|v_i|<1}(x_i - M_x)(1 - u_i^2)^2(y_i - M_y)(1 - v_i^2)^2
}{
    \\left(\\sum_{|u_i|<1}(1 - u_i^2)(1-5u_i^2)\\right)
    \\left(\\sum_{|v_i|<1}(1 - v_i^2)(1-5v_i^2)\\right)
}
```
where ``n_s`` is the sample size, ``M_x`` and ``M_y`` are the medians of ``\\{x_i\\}`` and
``\\{y_i\\}`` respectively, and
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

If either ``\\mathrm{MAD}_x = 0`` or ``\\mathrm{MAD}_y = 0``, the pairwise covariance is
defined to be zero.

The parameter ``c`` is a tuning constant, for which the default is ``9.0``.
Larger values will reduce the number of outliers that are removed — i.e. reducing
robustness, but increasing sample efficiency.

# Fields
- `c::Float64`: The tuning constant corresponding to ``c`` above.
- `modify_sample_size::Bool`: If `false`, then we use a sample size ``n_s`` equal to the
    total number of observations ``n``. This is consistent with the standard definition of
    biweight midcovariance in the literature. Otherwise, we count only those elements which
    are not rejected as outliers in the numerator, i.e. those for which ``|u_i|<1``
    and ``|v_i|<1``.
    This follows the implementation in astropy [2].

# Complexity
- Space: ``O(p^2)``
- Time: ``O(np^2)``

# References
[1] Beers, Flynn, and Gebhardt (1990; AJ 100, 32)
"Measures of Location and Scale for Velocities in Clusters of Galaxies -- A Robust Approach"

[2] [Astropy biweight_midcovariance](https://docs.astropy.org/en/stable/api/astropy.stats.biweight_midcovariance.html)

[3] [NIST Dataplot biweight midcovariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm)
"""
struct BiweightMidcovariance <: CovarianceEstimator
    c::Float64
    modify_sample_size::Bool
end
function BiweightMidcovariance(; c::Real=9.0, modify_sample_size::Bool=false)
    c > 0 || throw(ArgumentError("c must be positive, got $c"))
    return BiweightMidcovariance(Float64(c), modify_sample_size)
end

function covzm(ce::BiweightMidcovariance, x::AbstractVector{<:Real})
    MADx = median!(abs.(x))
    numerator = zero(eltype(x))

    # If MADx is zero, return zero.
    iszero(MADx) && return zero(one(eltype(x)) / 1)

    denominator = zero(eltype(x))
    count = 0
    for xi in x
        ui² = (xi / (ce.c * MADx))^2
        # Benchmarking suggests that having this branch is similar or faster to a form
        # using `ifelse`, where we always compute the updates but avoid the branch.
        ui² >= 1 && continue
        count += 1
        numerator += xi^2 * (1 - ui²)^4
        denominator += (1 - ui²) * (1 - 5 * ui²)
    end
    n = ifelse(ce.modify_sample_size, count, length(x))
    return n * numerator / (denominator^2)
end

function cov(
    ce::BiweightMidcovariance,
    x::AbstractVector{<:Real};
    mean::Union{Nothing,<:Real}=nothing,
)
    Mx = mean === nothing ? median(x) : mean
    return covzm(ce, x .- Mx)
end

function covzm(
    ce::BiweightMidcovariance, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}
)
    MADx = median!(abs.(x))
    MADy = median!(abs.(y))

    # Promote types between x & y for numerator
    numerator = zero(promote_type(eltype(x), eltype(y)))

    # If either of the MADs are zero, return zero.
    #   `numerator` is always zero here, and of the correct type.
    iszero(MADx) && return numerator
    iszero(MADy) && return numerator

    denominator_x = zero(eltype(x))
    denominator_y = zero(eltype(y))
    count = 0
    for (xi, yi) in zip(x, y)
        ui² = (xi / (ce.c * MADx))^2
        vi² = (yi / (ce.c * MADy))^2

        # Compute the updates, and increment if desired.
        # This way should avoid branches, and benchmarks suggest it's 1-2% faster than
        # using conditionals.
        ax = (1 - ui²) * (1 - 5 * ui²)
        ay = (1 - vi²) * (1 - 5 * vi²)
        anum = xi * (1 - ui²)^2 * yi * (1 - vi²)^2

        denominator_x += ifelse(ui² < 1, ax, zero(denominator_x))
        denominator_y += ifelse(vi² < 1, ay, zero(denominator_y))
        i_count, i_numerator = ifelse(ui² < 1 && vi² < 1, (1, anum), (0, zero(numerator)))
        count += i_count
        numerator += i_numerator
    end

    n = ifelse(ce.modify_sample_size, count, length(x))
    return n * numerator / (denominator_x * denominator_y)
end

function cov(
    ce::BiweightMidcovariance, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}
)
    if size(x) != size(y)
        throw(ArgumentError("x & y have different sizes, $(size(x)) & $(size(y))"))
    end

    return covzm(ce, x .- median(x), y .- median(y))
end

function cov(
    ce::BiweightMidcovariance,
    X::AbstractMatrix{<:Real};
    mean::Union{Nothing,AbstractVector{<:Real},AbstractMatrix{<:Real}}=nothing,
    dims::Int=1,
)
    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))

    # Ascertain the number of observations `n` and number of variables `p`.
    n, p = dims == 1 ? size(X) : (size(X, 2), size(X, 1))

    # Temporary storage to use when computing medians.
    temp = Vector{eltype(X)}(undef, n)

    # If the `mean` argument isn't provided, then we centralise with the median.
    # After this statement, Xc is standardised to be of shape (n, p)
    Xc = if mean === nothing
        # The following is equivalent to
        #   X .- median(X; dims=dims)
        # However, this implementation significantly reduces the runtime of this function.
        # This is primarily achieved by reducing the number of allocations -- we allocate a
        # single temporary buffer which we allow `median!` to mutate, and then copy the next
        # section of data into it.
        X = dims == 2 ? transpose(X) : X
        Xc = copy(X)
        @inbounds for i in 1:p
            # Avoiding copyto!(, ... @view(X[:, i])), because this can allocate a view if
            # the call to copyto! isn't inlined.
            for j in 1:n
                temp[j] = X[j, i]
            end
            Xc[:, i] .-= median!(temp)
        end
        Xc
    else
        Xc = X .- mean
        # Standardise the orientation of Xc.
        dims == 2 ? transpose(Xc) : Xc
    end

    # Compute all the median absolute deviations.
    # To avoid extra allocations, use a temporary buffer `temp` which will be populated with
    # the absolute values of the observations of each variable, and which we then allow
    # median! to mutate.
    T = typeof(one(eltype(X)) / 1)  # This is the output type of median.
    MAD = Matrix{T}(undef, 1, p)
    @inbounds for i in 1:p
        # Avoiding map!(abs, temp, @view(Xc[:, i])), because this can allocate a view if
        # the call to map! isn't inlined.
        for j in 1:n
            temp[j] = abs(Xc[j, i])
        end
        MAD[1, i] = median!(temp)
    end

    u² = @. (Xc / (ce.c * MAD))^2

    onesubu² = max.((1 .- u²), 0)
    a = @. Xc * onesubu²^2
    numerator = a' * a

    a[:, :] .= @. onesubu² * (1 - 5 * u²)  # Re-use `a` memory to avoid an allocation.
    b = sum(a; dims=1)
    denominator = b' * b

    n = if ce.modify_sample_size
        mask = onesubu² .> 0
        mask' * mask
    else
        n
    end

    # Re-use the numerator memory for the result.
    result = numerator
    result .*= n ./ denominator

    # Whereever the MAD is zero, we should set the result to zero.
    mad_zero = dropdims(iszero.(MAD); dims=1)
    result[mad_zero, :] .= zero(eltype(result))
    result[:, mad_zero] .= zero(eltype(result))

    return result
end
