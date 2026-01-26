export ConvergenceChecker

using LinearAlgebra: norm

"""
    ConvergenceChecker(;
        norm_condition,
        component_condition,
        condition_combiner,
        norm = LinearAlgebra.norm,
    )

Checks whether a sequence `val[0], val[1], val[2], ...` has converged to some
limit `L`, given the errors `err[iter] = val[iter] .- L`. This is done by
calling `is_converged!(::ConvergenceChecker, cache, val, err, iter)`, where
`val = val[iter]` and `err = err[iter]`. If the value of `L` is not known, `err`
can be an approximation of `err[iter]`. The `cache` for a `ConvergenceChecker`
can be obtained with `allocate_cache(::ConvergenceChecker, val_prototype)`,
where `val_prototype` is `similar` to `val` and `err`.

A `ConvergenceChecker` can perform two types of checks---it can check whether
`norm(val)` and `norm(err)` satisfy some `ConvergenceCondition`, and it can
check whether all the components of `abs.(val)` and `abs.(err)` individually
satisfy some `ConvergenceCondition`. These two checks can be combined with
either `&` or `|`. If one of the checks is not needed, the corresponding
`ConvergenceCondition` can be set to `nothing`.

Instead of `LinearAlgebra.norm`, `norm` can be set to anything that will convert
`val` and `err` to non-negative scalar values.
"""
Base.@kwdef struct ConvergenceChecker{
    NC <: Union{Nothing, ConvergenceCondition},
    CC <: Union{Nothing, ConvergenceCondition},
    C <: Union{typeof(&), typeof(|)},
    N,
}
    norm_condition::NC = nothing
    component_condition::CC = nothing
    condition_combiner::C = &
    norm::N = norm
end

function allocate_cache(alg::ConvergenceChecker, val_prototype)
    (; norm_condition, component_condition) = alg
    if isnothing(norm_condition) && isnothing(component_condition)
        error("ConvergenceChecker must have at least one ConvergenceCondition")
    end
    FT = eltype(val_prototype)
    return (;
        norm_cache = isnothing(norm_condition) ? nothing : Ref{cache_type(norm_condition, FT)}(),
        component_cache = isnothing(component_condition) ? nothing :
                          similar(val_prototype, cache_type(component_condition, FT)),
        component_bools = isnothing(component_condition) ? nothing : similar(val_prototype, Bool),
    )
end

function has_component_converged(alg, cache, val, err, iter)
    (; component_condition) = alg
    (; component_cache, component_bools) = cache
    # Caching abs.(val) and abs.(err) is probably not worth the overhead.
    @. component_bools = has_converged(
        (component_condition,), # wrap in a tuple to prevent iteration by .
        component_cache,
        abs(val),
        abs(err),
        iter,
    )
    return all(component_bools)
end

is_converged!(alg::Nothing, cache, val, err, iter) = false

function is_converged!(alg::ConvergenceChecker, cache, val, err, iter)
    (; norm_condition, component_condition, condition_combiner, norm) = alg
    (; norm_cache, component_cache) = cache
    if isnothing(norm_condition)
        converged = has_component_converged(alg, cache, val, err, iter)
    else
        norm_val = norm(val)
        norm_err = norm(err)
        converged = has_converged(norm_condition, norm_cache[], norm_val, norm_err, iter)
        if !isnothing(component_condition)
            if condition_combiner === &
                converged = converged && has_component_converged(alg, cache, val, err, iter)
            else # condition_combiner === |
                converged = converged || has_component_converged(alg, cache, val, err, iter)
            end
        end
    end
    if !converged # only update caches if they will be needed for future iters
        if !isnothing(norm_condition) && needs_cache_update(norm_condition, iter)
            norm_cache[] = updated_cache(norm_condition, norm_cache[], norm_val, norm_err, iter)
        end
        if !isnothing(component_condition) && needs_cache_update(component_condition, iter)
            @. component_cache = updated_cache((component_condition,), component_cache, abs(val), abs(err), iter)
        end
    end
    return converged
end
