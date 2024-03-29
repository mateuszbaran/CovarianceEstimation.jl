const Shrinkage = Union{Symbol, Real}
abstract type LinearShrinkageTarget end

# Taxonomy from https://strimmerlab.github.io/publications/journals/shrinkcov2005.pdf
# page 13

"""
    DiagonalUnitVariance

Target for linear shrinkage: unit matrix.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=1`` if ``i=j`` and
* ``F_{ij}=0`` otherwise
"""
struct DiagonalUnitVariance <: LinearShrinkageTarget end

"""
    DiagonalCommonVariance

Target for linear shrinkage: unit matrix multiplied by average variance of
variables.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=v`` if ``i=j`` with ``v=\\mathrm{tr}(S)/p`` and
* ``F_{ij}=0`` otherwise
"""
struct DiagonalCommonVariance <: LinearShrinkageTarget end

"""
    DiagonalUnequalVariance

Target for linear shrinkage: diagonal of covariance matrix.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=s_{ij}`` if ``i=j`` and
* ``F_{ij}=0`` otherwise
"""
struct DiagonalUnequalVariance <: LinearShrinkageTarget end

"""
    CommonCovariance

Target for linear shrinkage: see `target_C`.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=v`` if ``i=j`` with ``v=\\mathrm{tr}(S)/p`` and
* ``F_{ij}=c`` with ``c=\\sum_{i\\neq j} S_{ij}/(p(p-1))`` otherwise
"""
struct CommonCovariance <: LinearShrinkageTarget end

"""
    PerfectPositiveCorrelation

Target for linear shrinkage: see `target_E`.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=S_{ij}`` if ``i=j`` and
* ``F_{ij}=\\sqrt{S_{ii}S_{jj}}`` otherwise
"""
struct PerfectPositiveCorrelation <: LinearShrinkageTarget end

"""
    ConstantCorrelation

Target for linear shrinkage: see `target_F`.
A subtype of `LinearShrinkageTarget` where

* ``F_{ij}=S_{ij}`` if ``i=j`` and
* ``F_{ij}=\\overline{r}\\sqrt{S_{ii}S_{jj}}`` otherwise where
  ``\\overline{r}`` is the average sample correlation
"""
struct ConstantCorrelation <: LinearShrinkageTarget end

"""
    LinearShrinkage(target, shrinkage; corrected=false, drop_var0=false)

Linear shrinkage estimator described by equation
``(1 - \\lambda) S + \\lambda F`` where ``S`` is standard covariance matrix,
``F`` is shrinkage target described by argument `target` and ``\\lambda`` is a
shrinkage parameter, either given explicitly in `shrinkage` or automatically
determined according to one of the supported methods.

The corrected estimator is used if `corrected` is true.
`drop_var0=true` drops the zero-variance variables from the computation of `\\lambda`.
"""
struct LinearShrinkage{T<:LinearShrinkageTarget, S<:Shrinkage} <: CovarianceEstimator
    target::T
    shrinkage::S
    corrected::Bool
    drop_var0::Bool

    function LinearShrinkage(t::TT, s::SS; corrected=false, drop_var0=false) where {TT<:LinearShrinkageTarget, SS<:Real}
        0 ≤ s ≤ 1 || throw(ArgumentError("Shrinkage value should be between 0 and 1. Got $s."))
        new{TT, SS}(t, s, corrected, drop_var0)
    end
    function LinearShrinkage(t::TT, s::Symbol=:auto; corrected=false, drop_var0=false) where TT <: LinearShrinkageTarget
        s ∈ (:auto, :lw, :ss, :rblw, :oas) || throw(ArgumentError("Shrinkage method $s not supported."))
        new{TT, Symbol}(t, s, corrected, drop_var0)
    end
end

LinearShrinkage(;
    target::LinearShrinkageTarget=DiagonalUnitVariance(),
    shrinkage::Shrinkage,
    corrected::Bool=false,
    drop_var0::Bool=false) = LinearShrinkage(target, shrinkage, corrected=corrected, drop_var0=drop_var0)

"""
    cov(lse::LinearShrinkage, X, [weights::FrequencyWeights]; dims=1, mean=nothing)

Linear shrinkage covariance estimator for matrix `X` along dimension `dims`.
Computed using the method described by `lse`.

Optionally provide `weights` associated with each observation in `X` (see `StatsBase.FrequencyWeights`).

!!! note
    Theoretical guidance for the use of weights in shrinkage estimation seems sparse.
    `FrequencyWeights` have a straightforward implementation, but support for other `AbstractWeight` subtypes
    awaits analytical justification.
"""
function cov(lse::LinearShrinkage, X::AbstractMatrix{<:Real}, weights::FrequencyWeights...;
             dims::Int=1, mean=nothing)

    dims ∈ (1, 2) || throw(ArgumentError("Argument dims can only be 1 or 2 (given: $dims)"))

    Xc   = (dims == 1) ? copy(X) : copy(transpose(X))
    n, p = size(Xc)
    # sample covariance of size (p x p)
    S = cov(SimpleCovariance(corrected=lse.corrected), X, weights...; dims=dims, mean=mean)
    pvar = p - (lse.drop_var0 ? sum(iszero, diag(S)) : 0)

    # NOTE: don't need to check if mean is proper as this is already done above
    if mean === nothing
        Xc .-= Statistics.mean(Xc, weights...; dims=1)
    elseif mean isa AbstractArray
        if dims == 1
            Xc .-= mean
        else
            Xc .-= mean'
        end
    end

    return linear_shrinkage(lse.target, Xc, S, lse.shrinkage, n, p, pvar, lse.corrected, weights...)
end

##############################################################################
# Helper functions for Linear Shrinkage estimators

"""
    rescale(M, d)

Internal function to scale the rows and the columns of a square matrix `M`
according to the elements in `d`.
This is useful when dealing with the standardised data matrix which can be
written `Xs=Xc*D` where `Xc` is the centered data matrix and `D=Diagonal(d)`
so that its simple covariance is `D*S*D` where `S` is the simple covariance
of `Xc`. Such `D*M*D` terms appear often in the computations of optimal
shrinkage λ.

* Space complexity: ``O(p^2)``
* Time complexity: ``O(p^2)``
"""
rescale(M::AbstractMatrix, d::AbstractVector) = M .* d .* d'


"""
    rescale!(M, d)

Same as `rescale` but in place (no allocation).
"""
rescale!(M::AbstractMatrix, d::AbstractVector) = (M .*= d .* d')


"""
    uccov(X)

Internal function to compute `X*X'/n` where `n` is the number of rows of `X`.
This corresponds to the uncorrected covariance of `X` if `X` is centered.
This operation appears often in the computations of optimal shrinkage λ.

* Space complexity: ``O(p^2)``
* Time complexity: ``O(2np^2)``
"""
uccov(X::AbstractMatrix) = (X' * X) / size(X, 1)
uccov(X::AbstractMatrix, weights::FrequencyWeights) = (X' * (weights .* X)) / sum(weights)


"""
    sumij(S)

Internal function to compute the sum of elements of a square matrix `S`.
A keyword `with_diag` can be passed to indicate whether to include or not the
diagonal of `S` in the sum. Both cases happen often in the computations of
optimal shrinkage λ.

* Space complexity: ``O(1)``
* Time complexity: ``O(p^2)``
"""
function sumij(S::AbstractMatrix; with_diag=false)
    acc = sum(S)
    with_diag || return acc - tr(S)
    return acc
end


"""
    sumij2(S)

Internal function identical to `sumij` except that it passes the function
`abs2` to the sum so that it is the sum of the elements of `S` squared
which is computed. This is significantly more efficient than using
`sumij(S.^2)` for large matrices as it allocates very little.

* Space complexity: ``O(1)``
* Time complexity: ``O(2p^2)``
"""
function sumij2(S::AbstractMatrix; with_diag=false)
    acc = sum(abs2, S)
    with_diag || return acc - sum(abs2, diag(S))
    return acc
end


"""
    sum_fij(Xc, S, n, κ)

Internal function corresponding to ``∑_{i≂̸j}f_{ij}`` that appears in
https://strimmerlab.github.io/publications/journals/shrinkcov2005.pdf p.11.

* Space complexity: ``O(np + 2p^2)``
* Time complexity: ``O(2np^2)``
"""
function sum_fij(Xc, S, n, κ)
    sd  = sqrt.(diag(S))
    sdinv = map(z -> guardeddiv(1, z), sd)
    M   = ((Xc.^3)' * Xc) .* sdinv
    M .-= κ .* S .* sd
    M .*= sd'
    return sumij(M) / (n * κ)
end
function sum_fij(Xc, S, n, κ, weights)
    sd  = sqrt.(diag(S))
    sdinv = map(z -> guardeddiv(1, z), sd)
    M   = ((Xc.^3)' * (weights .* Xc)) .* sdinv
    M .-= κ .* S .* sd
    M .*= sd'
    return sumij(M) / (sum(weights) * κ)
end
##############################################################################

"""
    linear_shrinkage(target, Xc, S, λ, n, p, pvar, corrected, [weights])

Performs linear shrinkage with target of type `target` for data matrix `Xc`
of size `n` by `p` with covariance matrix `S` and shrinkage parameter `λ`.
Calculates corrected covariance if `corrected` is true.

`pvar == p` or `pvar = p - sum(iszero, diag(S))`, the number of non-zero
diagonal variances in `S`. The choice is controlled by `LinearShrinkage(...; drop_var0=true/false)`.
"""
linear_shrinkage

## TARGET A

function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    return linshrink(I, S, λ)
end

"""
    linear_shrinkage(::DiagonalUnitVariance, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a `DiagonalUnitVariance`.
"""
function linear_shrinkage(::DiagonalUnitVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F   = I
    T   = float(eltype(S))
    wn = totalweight(n, weights...)
    κ   = wn - Int(corrected)
    γ   = T(κ/wn)
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS² = sumij2(S, with_diag=true)
        λ   = sumij(uccov(Xc², weights...), with_diag=true) / γ^2 - ΣS²
        λ  /= κ * (ΣS² - 2tr(S) + pvar)
    elseif λ == :ss
        # use the standardised data matrix
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        S̄   = rescale(S, sqrt.(d)) # this has diagonal 1/κ
        ΣS̄² = sumij2(S̄, with_diag=true)
        λ   = sumij(rescale!(uccov(Xc², weights...), d), with_diag=true) / γ^2 - ΣS̄²
        λ  /= T(κ * ΣS̄² - pvar / κ)
    else
        throw(ArgumentError("Unsupported shrinkage method for target DiagonalUnitVariance: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink(F, S, λ)
end

## TARGET B

target_B(S::AbstractMatrix, pvar::Int) = tr(S)/pvar * I

function linear_shrinkage(::DiagonalCommonVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    return linshrink(target_B(S, pvar), S, λ)
end

"""
    linear_shrinkage(::DiagonalCommonVariance, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a `DiagonalCommonVariance`.
"""
function linear_shrinkage(::DiagonalCommonVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F   = target_B(S, pvar)
    T   = float(eltype(F))
    wn = totalweight(n, weights...)
    κ   = wn - Int(corrected)
    γ   = T(κ/wn)
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        v   = F.λ # tr(S)/p
        ΣS² = sumij2(S, with_diag=true)
        λ   = sumij(uccov(Xc², weights...), with_diag=true) / γ^2 - ΣS²
        λ  /= κ * (ΣS² - pvar*v^2)
    elseif λ == :ss
        # use the standardised data matrix
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        S̄   = rescale(S, sqrt.(d)) # this has diagonal 1/κ
        v̄   = κ # tr(S̄)/p
        ΣS̄² = sumij2(S̄, with_diag=true)
        λ   = sumij(rescale!(uccov(Xc², weights...), d), with_diag=true) / γ^2 - ΣS̄²
        λ  /= T(κ * ΣS̄² - pvar/κ)
    elseif λ == :rblw
        # https://arxiv.org/pdf/0907.4698.pdf equations 17, 19
        trS² = sum(abs2, S)
        tr²S = tr(S)^2
        # note: using corrected or uncorrected S does not change λ
        λ = T(((wn-2)/wn * trS² + tr²S) / ((wn+2) * (trS² - tr²S/pvar)))
    elseif λ == :oas
        # https://arxiv.org/pdf/0907.4698.pdf equation 23
        trS² = sum(abs2, S)
        tr²S = tr(S)^2
        # note: using corrected or uncorrected S does not change λ
        λ = ((one(T)-T(2.0)/pvar) * trS² + tr²S) / ((wn+one(T)-T(2.0)/pvar) * (trS² - tr²S/pvar))
    else
        throw(ArgumentError("Unsupported shrinkage method for target DiagonalCommonVariance: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink(F, S, λ)
end

## TARGET D

target_D(S::AbstractMatrix) = Diagonal(S)

function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    return linshrink(target_D(S), S, λ)
end

"""
    linear_shrinkage(::DiagonalUnequalVariance, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a `DiagonalUnequalVariance`.
"""
function linear_shrinkage(::DiagonalUnequalVariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F   = target_D(S)
    T   = float(eltype(F))
    wn = totalweight(n, weights...)
    κ   = wn - Int(corrected)
    γ   = T(κ / wn)
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS² = sumij2(S)
        λ   = sumij(uccov(Xc², weights...)) / γ^2 - ΣS²
        λ  /= κ * ΣS²
    elseif λ == :ss
        keep = diag(S) .> zero(T)
        Xc² = Xc²[:, keep]
        # use the standardised data matrix
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        ΣS̄² = sumij2(rescale(S[keep, keep], sqrt.(d)))
        λ   = sumij(rescale!(uccov(Xc², weights...), d)) / γ^2 - ΣS̄²
        λ  /= κ * ΣS̄²
    else
        throw(ArgumentError("Unsupported shrinkage method for target DiagonalUnequalVariance: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink(F, S, λ)
end

## TARGET C

function target_C(S::AbstractMatrix, p::Int, pvar::Int)
    v  = tr(S)/pvar
    c  = sumij(S; with_diag=false) / (pvar * (pvar - 1))
    F  = fill(c, (p, p))
    F -= Diagonal(F)
    F += v * I
    return F, v, c
end

function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F, _, _ = target_C(S, p, pvar)
    return linshrink!(F, S, λ)
end

"""
    linear_shrinkage(::CommonCovariance, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a `CommonCovariance`.
"""
function linear_shrinkage(::CommonCovariance, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F, v, c = target_C(S, p, pvar)
    T   = float(eltype(F))
    wn = totalweight(n, weights...)
    κ   = wn - Int(corrected)
    γ   = T(κ/wn)
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS² = sumij2(S, with_diag=true)
        λ   = sumij(uccov(Xc², weights...), with_diag=true) / γ^2 - ΣS²
        λ  /= κ * (ΣS² - pvar*(pvar-1)*c^2 - pvar*v^2)
    elseif λ == :ss
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        S̄   = rescale(S, sqrt.(d))
        ΣS̄² = sumij2(S̄, with_diag=true)
        λ   = sumij(rescale!(uccov(Xc², weights...), d), with_diag=true) / γ^2 - ΣS̄²
        λ  /= κ * ΣS̄² - pvar/κ - κ * sumij(S̄; with_diag=false)^2 / (pvar * (pvar - 1))
    else
        throw(ArgumentError("Unsupported shrinkage method for target CommonCovariance: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink!(F, S, λ)
end

## TARGET E

function target_E(S::AbstractMatrix)
    d = sqrt.(diag(S))
    return d*d'
end

function linear_shrinkage(::PerfectPositiveCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    return linshrink!(target_E(S), S, λ)
end

"""
    linear_shrinkage(::PerfectPositiveCorrelation, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a
`PerfectPositiveCorrelation`.
"""
function linear_shrinkage(::PerfectPositiveCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F   = target_E(S)
    T   = float(eltype(F))
    wn = totalweight(n, weights...)
    κ   = wn - Int(corrected)
    γ   = T(κ/wn)
    Xc² = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS² = sumij2(S)
        λ   = (sumij(uccov(Xc², weights...)) / γ^2 - ΣS²) / κ
        λ  -= sum_fij(Xc, S, n, κ, weights...)
        λ  /= sumij2(S - F)
    elseif λ == :ss
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        s   = sqrt.(d)
        S̄   = rescale(S, s)
        ΣS̄² = sumij2(S̄)
        λ   = (sumij(rescale!(uccov(Xc², weights...), d)) / γ^2 - ΣS̄²) / κ
        λ  -= sum_fij(Xc .* s', S̄, n, κ, weights...)
        F̄  = target_E(S̄)
        λ  /= sumij2(S̄ - F̄)
    else
        throw(ArgumentError("Unsupported shrinkage method for target PerfectPositiveCorrelation: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink!(F, S, λ)
end

## TARGET F

function target_F(S::AbstractMatrix, p::Int, pvar::Int)
    s  = sqrt.(diag(S))
    sinv = pvar < p ? map(z -> guardeddiv(1, z), s) : 1 ./ s
    sinv_ = sinv*sinv'
    r̄  = (sum(S .* sinv_) - pvar) / (pvar * (pvar - 1))
    F_ = r̄ * (s*s')
    F  = F_ + (Diagonal(S) - Diagonal(F_))
    return F, r̄
end

function linear_shrinkage(::ConstantCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Real, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F, _ = target_F(S, p, pvar)
    return linshrink!(F, S, λ)
end

"""
    linear_shrinkage(::ConstantCorrelation, Xc, S, λ, n, p, pvar, corrected)

Compute the shrinkage estimator where the target is a `ConstantCorrelation`.
"""
function linear_shrinkage(::ConstantCorrelation, Xc::AbstractMatrix,
                          S::AbstractMatrix, λ::Symbol, n::Int, p::Int, pvar::Int,
                          corrected::Bool, weights::FrequencyWeights...)

    F, r̄ = target_F(S, p, pvar)
    T    = float(eltype(F))
    wn = totalweight(n, weights...)
    κ    = wn - Int(corrected)
    γ    = T(κ/wn)
    Xc²  = Xc.^2
    # computing the shrinkage
    if λ ∈ [:auto, :lw]
        ΣS² = sumij2(S)
        λ   = (sumij(uccov(Xc², weights...)) / γ^2 - ΣS²) / κ
        λ  -= r̄ * sum_fij(Xc, S, n, κ, weights...)
        λ  /= sumij2(S - F)
    elseif λ == :ss
        d   = diaginv(pvar < p, oneunit(T), vec(sum(Xc², weights...; dims=1)))
        s    = sqrt.(d)
        S̄    = rescale(S, s)
        F̄, r̄ = target_F(S̄, p, pvar)
        ΣS̄²  = sumij2(S̄)
        λ    = (sumij(rescale!(uccov(Xc², weights...), d)) / γ^2 - ΣS̄²) / κ
        λ   -= r̄ * sum_fij(Xc .* s', S̄, n, κ, weights...)
        λ   /= sumij2(S̄ - F̄)
    else
        throw(ArgumentError("Unsupported shrinkage method for target ConstantCorrelation: $λ."))
    end
    λ = clamp(λ, zero(T), one(T))
    return linshrink!(F, S, λ)
end
