export ConvergenceCondition,
    MaximumError, MaximumRelativeError, MaximumErrorReduction, MinimumRateOfConvergence, MultipleConditions

"""
    ConvergenceCondition

An abstract type for objects that can check whether a sequence of non-negative
scalar values `val[0], val[1], val[2], ...` has converged to some limit `L`,
given the errors `err[iter] = |val[iter] - L|`.

Every subtype of `ConvergenceCondition` must define
    `has_converged(::ConvergenceCondition, cache, val, err, iter)`.
The `cache`, which is set to `nothing` by default, may be used to store
information from previous iterations that is useful for determining convergence.
In order to have access to a `cache` of some particular type, a subtype of
`ConvergenceCondition` should define
    `cache_type(::ConvergenceCondition, ::Type{FT})`.
To specify on which iterations this cache should be updated, it should define
    `needs_cache_update(::ConvergenceCondition, iter)`.
To specify how the cache should be update on those iterations, it should define
    `updated_cache(::ConvergenceCondition, cache, val, err, iter)`.

Although `cache_type` can call `promote_type` to prevent potential type
instability errors, this should be avoided in order to ensure that users write
type-stable code.
"""
abstract type ConvergenceCondition end
cache_type(::ConvergenceCondition, ::Type{FT}) where {FT} = Nothing
needs_cache_update(::ConvergenceCondition, iter) = false

"""
    MaximumError(max_err)

Checks whether `err[iter] ≤ max_err`. Since `err[iter] ≥ 0`, this can only be
`true` if `max_err ≥ 0`.
"""
struct MaximumError{FT} <: ConvergenceCondition
    max_err::FT
end
has_converged((; max_err)::MaximumError, cache, val, err, iter) = err <= max_err

"""
    MaximumRelativeError(max_rel_err)

Checks whether `err[iter] ≤ max_rel_err * val[iter]`. Since `err[iter] ≥ 0` and
`val[iter] ≥ 0`, this can only be `true` if `max_rel_err ≥ 0`.
"""
struct MaximumRelativeError{FT} <: ConvergenceCondition
    max_rel_err::FT
end
has_converged((; max_rel_err)::MaximumRelativeError, cache, val, err, iter) = err <= max_rel_err * val

"""
    MaximumErrorReduction(max_reduction)

Checks whether `err[iter] ≤ max_reduction * err[0]` for all `iter ≥ 1`. Since
`err[iter] ≥ 0`, this can only be `true` if `max_reduction ≥ 0`. Also, it must
be the case that `max_reduction ≤ 1` in order for the sequence to not diverge
(i.e., to avoid `err[iter] > err[0]`).
"""
struct MaximumErrorReduction{FT} <: ConvergenceCondition
    max_reduction::FT
end
cache_type(::MaximumErrorReduction, ::Type{FT}) where {FT} = NamedTuple{(:max_err,), Tuple{FT}}
has_converged(::MaximumErrorReduction, cache, val, err, iter) = iter >= 1 && err <= cache.max_err
needs_cache_update(::MaximumErrorReduction, iter) = iter == 0
updated_cache((; max_reduction)::MaximumErrorReduction, cache, val, err, iter) = (; max_err = max_reduction * err)

"""
    MinimumRateOfConvergence(rate, order = 1)

Checks whether `err[iter] ≥ rate * err[iter - 1]^order` for all `iter ≥ 1`.
Since `err[iter] ≥ 0`, this can only be `true` if `rate ≥ 0`. Also, if
`order == 1`, it must be the case that `rate ≤ 1` in order for the sequence to
not diverge (i.e., to avoid `err[iter] > err[iter - 1]`). In addition, if
`err[iter] < 1` for all sufficiently large values of `iter`, it must be the case
that `order ≥ 1` for the sequence to not diverge.
"""
struct MinimumRateOfConvergence{FT, FT2} <: ConvergenceCondition
    rate::FT
    order::FT2
end
MinimumRateOfConvergence(rate) = MinimumRateOfConvergence(rate, 1)
cache_type(::MinimumRateOfConvergence, ::Type{FT}) where {FT} = NamedTuple{(:min_err,), Tuple{FT}}
has_converged(::MinimumRateOfConvergence, cache, val, err, iter) = iter >= 1 && err >= cache.min_err
needs_cache_update(::MinimumRateOfConvergence, iter) = true
updated_cache((; rate, order)::MinimumRateOfConvergence, cache, val, err, iter) = (; min_err = rate * err^order)

"""
    MultipleConditions(condition_combiner = all, conditions...)

Checks multiple `ConvergenceCondition`s, combining their results with either
`all` or `any`.
"""
struct MultipleConditions{CC, C} <: ConvergenceCondition
    condition_combiner::CC
    conditions::C
    function MultipleConditions(
        condition_combiner::Union{typeof(all), typeof(any)},
        conditions::Vararg{ConvergenceCondition},
    )
        return new{typeof(condition_combiner), typeof(conditions)}(condition_combiner, conditions)
    end
end

MultipleConditions(conditions::ConvergenceCondition...) = MultipleConditions(all, conditions)
cache_type((; conditions)::MultipleConditions, ::Type{FT}) where {FT} =
    Tuple{map(condition -> cache_type(condition, FT), conditions)...}
has_converged(alg::MultipleConditions, caches, val, err, iter) = alg.condition_combiner(
    ((condition, cache),) -> has_converged(condition, cache, val, err, iter),
    zip(alg.conditions, caches),
)
needs_cache_update((; conditions)::MultipleConditions, iter) =
    any(condition -> needs_cache_update(condition, iter), conditions)
updated_cache((; conditions)::MultipleConditions, caches, val, err, iter) = map(
    (condition, cache) ->
        needs_cache_update(condition, iter) ? updated_cache(condition, cache, val, err, iter) : cache,
    conditions,
    caches,
)
