abstract type SymbolicTypeTrait end

"""
    struct ScalarSymbolic <: SymbolicTypeTrait end

Trait indicating a type is a scalar symbolic variable.

See also: [`ArraySymbolic`](@ref), [`NotSymbolic`](@ref), [`symbolic_type`](@ref)
"""
struct ScalarSymbolic <: SymbolicTypeTrait end

"""
    struct ArraySymbolic <: SymbolicTypeTrait end

Trait indicating type is a symbolic array. Calling `collect` on a symbolic array must
return an `AbstractArray` containing `ScalarSymbolic` variables for each element in the
array, in the same shape as the represented array. For example, if `a` is a symbolic array
representing a 2x2 matrix, `collect(a)` must return a 2x2 array of scalar symbolic variables.

See also: [`ScalarSymbolic`](@ref), [`NotSymbolic`](@ref), [`symbolic_type`](@ref)
"""
struct ArraySymbolic <: SymbolicTypeTrait end

"""
    struct NotSymbolic <: SymbolicTypeTrait end

Trait indicating a type is not symbolic.

See also: [`ScalarSymbolic`](@ref), [`ArraySymbolic`](@ref), [`symbolic_type`](@ref)
"""
struct NotSymbolic <: SymbolicTypeTrait end

"""
    symbolic_type(x) = symbolic_type(typeof(x))
    symbolic_type(::Type)

Get the symbolic type trait of a type. Default to [`NotSymbolic`](@ref) for all types
except `Symbol` and `Expr`, both of which are [`ScalarSymbolic`](@ref).

See also: [`ScalarSymbolic`](@ref), [`ArraySymbolic`](@ref), [`NotSymbolic`](@ref)
"""
symbolic_type(x) = symbolic_type(typeof(x))
symbolic_type(::Type) = NotSymbolic()
symbolic_type(::Type{Symbol}) = ScalarSymbolic()
symbolic_type(::Type{Expr}) = ScalarSymbolic()

"""
    hasname(x)

Check whether the given symbolic variable (for which `symbolic_type(x) != NotSymbolic()`) has a valid name as per `getname`. Defaults to `true` for `x::Symbol`.
"""
function hasname end

hasname(::Symbol) = true
hasname(::Any) = false

"""
    getname(x)::Symbol

Get the name of a symbolic variable as a `Symbol`. Acts as the identity function for
`x::Symbol`.
"""
function getname end
getname(x::Symbol) = x

"""
    symbolic_evaluate(expr, syms::Dict; kwargs...)

Return the value of symbolic expression `expr` where the values of variables involved are
obtained from the dictionary `syms`. The keys of `syms` are symbolic variables (not
expressions of variables). The values of `syms` can be values or symbolic
expressions.

The returned value should either be a value or an expression involving symbolic variables
not present as keys in `syms`.

The function can take additional keyword arguments to control implementation-specific
behavior.

This is already implemented for 
`symbolic_evaluate(expr::Union{Symbol, Expr}, syms::Dict)`.
"""
function symbolic_evaluate(expr::Union{Symbol, Expr}, syms::Dict)
    while (newexpr = _symbolic_evaluate_helper(expr, syms)) != expr
        expr = newexpr
    end
    return try
        eval(expr)
    catch
        expr
    end
end

function _symbolic_evaluate_helper(expr, syms::Dict)
    if (res = get(syms, expr, nothing)) !== nothing
        return res
    end
    expr isa Expr || return expr

    newexpr = Expr(expr.head)
    sizehint!(newexpr.args, length(expr.args))
    for arg in expr.args
        push!(newexpr.args, _symbolic_evaluate_helper(arg, syms))
    end
    newexpr
end

############ IsTimeseriesTrait

abstract type IsTimeseriesTrait end

"""
    struct Timeseries <: IsTimeseriesTrait end

Trait indicating a type contains timeseries data. This affects the behaviour of
functions such as [`state_values`](@ref) and [`current_time`](@ref).

See also: [`NotTimeseries`](@ref), [`is_timeseries`](@ref)
"""
struct Timeseries <: IsTimeseriesTrait end

"""
    struct NotTimeseries <: IsTimeseriesTrait end

Trait indicating a type does not contain timeseries data. This affects the behaviour
of functions such as [`state_values`](@ref) and [`current_time`](@ref). Note that
if a type is `NotTimeseries` this only implies that it does not _store_ timeseries
data. It may still be time-dependent. For example, an `ODEProblem` only stores
the initial state of a system, so it is `NotTimeseries`, but still time-dependent.
This is the default trait variant for all types.

See also: [`Timeseries`](@ref), [`is_timeseries`](@ref).
"""
struct NotTimeseries <: IsTimeseriesTrait end

"""
    is_timeseries(x) = is_timeseries(typeof(x))
    is_timeseries(::Type)

Get the timeseries trait of a type. Defaults to [`NotTimeseries`](@ref) for all types.
A type for which `is_timeseries(T) == Timeseries()` may also have a parameter timeseries.
This is determined by the [`is_parameter_timeseries`](@ref) trait.

See also: [`Timeseries`](@ref), [`NotTimeseries`](@ref), [`is_parameter_timeseries`](@ref).
"""
function is_timeseries end

is_timeseries(x) = is_timeseries(typeof(x))
is_timeseries(::Type) = NotTimeseries()

"""
    is_parameter_timeseries(x) = is_parameter_timeseries(typeof(x))
    is_parameter_timeseries(::Type)

Get the parameter timeseries trait of a type. Defaults to [`NotTimeseries`](@ref) for all
types. A type for which `is_parameter_timeseries(T) == Timeseries()` must also have
`is_timeseries(T) == Timeseries()`.

See also: [`Timeseries`](@ref), [`NotTimeseries`](@ref), [`is_timeseries`](@ref).
"""
function is_parameter_timeseries end

is_parameter_timeseries(x) = is_parameter_timeseries(typeof(x))
is_parameter_timeseries(::Type) = NotTimeseries()
