# Nonlinear shrinkage estimation can be described in terms of a loss function measuring the error between
# the target and the sample eigenvalues. In:
#    Donoho, D.L., Gavish, M. and Johnstone, I.M., 2018.
#    Optimal shrinkage of eigenvalues in the spiked covariance model. Annals of statistics, 46(4), p.1742.
# there is a systematic analysis of different loss functions; their classification scheme is encoded here.

# Implementation note:
# While we could parametrize all these different loss functions in the type system and use dispatch,
# that would induce needless specialization: every function that took a `LossFunction` would have to be
# specialized, even if the the argument encoding the loss function is merely "passed-through."
# So instead, we hide details from the type system, dividing the Donoho classification scheme into just two types,
# `NormLossCov` and `StatLossCov`, and use a fast runtime check to determine which shrinkage function to use.

abstract type LossFunction end

"""
    NormLossCov(norm::Symbol, pivotidx::Int)

Specify a loss function for which the estimated covariance will be optimal. `norm` is one of
`:L1`, `:L2`, or `:Linf`, and `pivotidx` is an integer from 1 to 7, as specified in Table 1 (p. 1755)
of Donoho et al. (2018). In the table below, `A` and `B` are the target and sample covariances,
respectively, and the loss function is the specified norm on the quantity in the `pivot` column:

| `pivotidx` | `pivot` | Notes |
|------------|---------|-------|
| 1          | `A - B`     | |
| 2          | `A⁻¹ - B⁻¹`     | |
| 3          | `A⁻¹ B - I`     | Not available for `:L1` |
| 4          | `B⁻¹ A - I`     | Not available for `:L1` |
| 5          | `A⁻¹ B + B⁻¹ A - 2I` | Not supported |
| 6          | `sqrt(A) \\ B / sqrt(A) - I` | |
| 7          | `log(sqrt(A) \\ B / sqrt(A))` | Not supported |

See also [`StatLossCov`](@ref).

Reference:
    Donoho, D.L., Gavish, M. and Johnstone, I.M., 2018.
    Optimal shrinkage of eigenvalues in the spiked covariance model. Annals of statistics, 46(4), p.1742.
"""
struct NormLossCov <: LossFunction
    # Lᴺᴷ where N is the norm and K is an integer (1 through 7) representing the pivot function
    norm::Symbol
    pivotidx::Int

    function NormLossCov(norm::Symbol, pivotidx::Int)
        norm ∈ (:L1, :L2, :Linf) || throw(ArgumentError("norm must be :L1, :L2, or :Linf"))
        1 <= pivotidx <= 7 || throw(ArgumentError("pivotidx must be from 1 to 7 (see Table 1 in Donoho et al. (2018))"))
        return new(norm, pivotidx)
    end
end

"""
    StatLossCov(mode::Symbol)

Specify a loss function for which the estimated covariance will be optimal. `mode` is one of
`:st`, `:ent`, `:div`, `:aff`, or `:fre`, as specified in Table 2 (p. 1757) of Donoho et al. (2018).
In the table below, `A` and `B` are the target and sample covariances, respectively:

| `mode` | loss |  Interpretation  |
|--------|---------|-----|
| `:st`  | `st(A, B) = tr(A⁻¹ B - I) - log(det(B)/det(A))` | Minimize KL-divergence between `N(0, A)` and `N(0, B)` where `N` is normal distribution |
| `:ent` | `st(B, A)` | Minimize errors in Mahalanobis distances |
| `:div` | `st(A, B) + st(B, A)` | |
| `:aff` | `0.5 * log(det(A + B) / (2 * sqrt(det(A*B))))` | Minimize Hellinger distance between `N(0, A)` and `N(0, B)` |
| `:fre` | `tr(A + B - 2sqrt(A*B))` | |
"""
struct StatLossCov <: LossFunction
    mode::Symbol

    function StatLossCov(mode::Symbol)
        statlosses = (:st, :ent, :div, :aff, :fre)

        mode ∈ statlosses || throw(ArgumentError("mode must be among $(statlosses)"))
        return new(mode)
    end
end


# Implement Table 2, Donoho et al. (2018), p. 1757

function shrinker(loss::NormLossCov, ℓ::Real, c::Real, s::Real)
    # See top of file for why these are branches rather than dispatch
    norm, pivotidx = loss.norm, loss.pivotidx
    pivotidx ∈ (5, 7) && throw(ArgumentError("Pivot index $(pivotidx) is not supported, see Table 2 in Donoho et al. 2018"))
    if norm == :L2         # Frobenius
        return  pivotidx == 1 ? ℓ * c^2 + s^2 :
                pivotidx == 2 ? ℓ / (c^2 + ℓ * s^2) :
                pivotidx == 3 ? (ℓ * c^2 + ℓ^2 * s^2) / (c^2 + ℓ^2 * s^2) :
                pivotidx == 4 ? (ℓ^2 * c^2 + s^2) / (ℓ * c^2 + s^2) :
                #= pivotidx == 6 =# 1 + (ℓ - 1) * c^2 / (c^2 + ℓ * s^2)^2
    elseif norm == :Linf   # Operator
        pivotidx ∈ (3, 4) && throw(ArgumentError("Pivot index $(pivotidx) is not supported for Linf norm, see Table 2 in Donoho et al. 2018"))
        return  pivotidx ∈ (1, 2) ? ℓ :
                #= pivotidx == 6 =# 1 + (ℓ - 1) / (c^2 + ℓ * s^2)
    elseif norm == :L1     # Nuclear
        val = pivotidx == 1 ? 1 + (ℓ - 1) * (1 - 2s^2) :
              pivotidx == 2 ? ℓ / (c^2 + (2ℓ-1)*s^2) :
              pivotidx == 3 ? ℓ / (c^2 + ℓ^2*s^2) :
              pivotidx == 4 ? (ℓ^2*c^2 + s^2) / ℓ :
              #= pivotidx == 6 =# (ℓ - (ℓ - 1)^2*c^2*s^2) / (c^2 + ℓ*s^2)^2
        return max(val, 1)
    end
    throw(ArgumentError("Norm $(norm) is not supported"))
end

function shrinker(loss::StatLossCov, ℓ::Real, c::Real, s::Real)
    mode = loss.mode
    if mode == :st
        return ℓ / (c^2 + ℓ * s^2)
    elseif mode == :ent
        return ℓ * c^2 + s^2
    elseif mode == :div
        return sqrt((ℓ^2 * c^2 + ℓ * s^2) / (c^2 + ℓ * s^2))
    elseif mode == :fre
        return (sqrt(ℓ) * c^2 + s^2)^2
    elseif mode == :aff
        return ((1 + c^2)*ℓ + s^2) / (1 + c^2 + ℓ * s^2)
    end
    throw(ArgumentError("Mode $(mode) is not supported"))
end
