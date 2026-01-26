"""
    symbolic_container(indp)

Using `indp`, return an object that implements the index provider interface. In case `indp`
itself implements the interface, `indp` can be returned as-is. All index provider interface
methods fall back to calling the same method on `symbolic_container(indp)`, so this may be
used for trivial implementations of the interface that forward all calls to another object.

Note that this method is optional. Thus the correct method to check for a fallback is:
```julia
hasmethod(symbolic_container, Tuple{typeof(indp)}) && symbolic_container(indp) != indp
```
"""
function symbolic_container end

"""
    is_variable(indp, sym)

Check whether the given `sym` is a variable in `indp`.
"""
is_variable(indp, sym) = is_variable(symbolic_container(indp), sym)

"""
    variable_index(indp, sym, [i])

Return the index of the given variable `sym` in `indp`, or `nothing` otherwise. If
[`constant_structure`](@ref) is `false`, this accepts the current time index as an
additional parameter `i`.
"""
variable_index(indp, sym) = variable_index(symbolic_container(indp), sym)
variable_index(indp, sym, i) = variable_index(symbolic_container(indp), sym, i)

"""
    variable_symbols(indp, [i])

Return a vector of the symbolic variables being solved for in the index provider `indp`.
If `constant_structure(sys) == false` this accepts an additional parameter indicating
the current time index. The returned vector should not be mutated.

For types that implement `Base.getindex` with symbolic indices using this interface,
the shorthand `valp[solvedvariables]` can be used as shorthand for
`valp[variable_symbols(sys)]`. See: [`solvedvariables`](@ref).
"""
variable_symbols(indp) = variable_symbols(symbolic_container(indp))
variable_symbols(indp, i) = variable_symbols(symbolic_container(indp), i)

"""
    is_parameter(indp, sym)

Check whether the given `sym` is a parameter in `indp`.
"""
is_parameter(indp, sym) = is_parameter(symbolic_container(indp), sym)

"""
    parameter_index(indp, sym)

Return the index of the given parameter `sym` in `indp`, or `nothing` otherwise.
"""
parameter_index(indp, sym) = parameter_index(symbolic_container(indp), sym)

"""
    is_timeseries_parameter(indp, sym)

Check whether the given `sym` is a timeseries parameter in `indp`.
"""
function is_timeseries_parameter(indp, sym)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        is_timeseries_parameter(sc, sym)
    else
        return false
    end
end

"""
    struct ParameterTimeseriesIndex
    function ParameterTimeseriesIndex(timeseries_idx, parameter_idx)

A struct storing the index of the timeseries of a timeseries parameter in a parameter
timeseries object. `timeseries_idx` refers to an index that identifies the timeseries
that the parameter belongs to. `parameter_idx` refers to the index of the parameter's
timeseries in that timeseries object. Note that `parameter_idx` may be different from
the object returned by [`parameter_index`](@ref) for a given parameter. The two fields in
this struct are `timeseries_idx` and `parameter_idx`.
"""
struct ParameterTimeseriesIndex{T, I}
    timeseries_idx::T
    parameter_idx::I
end

"""
    timeseries_parameter_index(indp, sym)

Return the index of timeseries parameter `sym` in `indp`. Must return this index as a
[`ParameterTimeseriesIndex`](@ref) object. Return `nothing` if `sym` is not a timeseries
parameter in `indp`. Defaults to returning `nothing`. Respects the
[`symbolic_container`](@ref) fallback for `indp` if present.
"""
function timeseries_parameter_index(indp, sym)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        timeseries_parameter_index(symbolic_container(indp), sym)
    else
        return nothing
    end
end

"""
    parameter_observed(indp, sym)

Return the observed function of `sym` in `indp`. This functions similarly to
[`observed`](@ref) except that `u` is not an argument of the returned function. For time-
dependent systems, the returned function must have the signature `(p, t) -> [values...]`.
For time-independent systems, the returned function must have the signature
`(p) -> [values...]`.

By default, this function returns `nothing`, indicating that the index provider does not
support generating parameter observed functions.
"""
function parameter_observed(indp, sym)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        return parameter_observed(symbolic_container(indp), sym)
    else
        return nothing
    end
end

"""
    struct ContinuousTimeseries end

A singleton struct corresponding to the timeseries index of the continuous timeseries.
"""
struct ContinuousTimeseries end

"""
    get_all_timeseries_indexes(indp, sym)

Return a `Set` of all unique timeseries indexes of variables in symbolic variable
`sym`. `sym` may be a symbolic variable or expression, an array of symbolics, an index,
or an array of indices. Continuous variables correspond to the
[`ContinuousTimeseries`](@ref) timeseries index. Non-timeseries parameters do not have a
timeseries index. Timeseries parameters have the same timeseries index as that returned by
[`timeseries_parameter_index`](@ref). Note that the independent variable corresponds to
the `ContinuousTimeseries` timeseries index.

Any ambiguities should be resolved in favor of variables. For example, if `1` could refer
to the variable at index `1` or parameter at index `1`, it should be interpreted as the
variable.

By default, this function returns `Set([ContinuousTimeseries()])`.
"""
function get_all_timeseries_indexes(indp, sym)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        return get_all_timeseries_indexes(symbolic_container(indp), sym)
    else
        return Set([ContinuousTimeseries()])
    end
end

"""
    parameter_symbols(indp)

Return a vector of the symbolic parameters of the given index provider `indp`. The returned
vector should not be mutated.
"""
parameter_symbols(indp) = parameter_symbols(symbolic_container(indp))

"""
    is_independent_variable(indp, sym)

Check whether the given `sym` is an independent variable in `indp`. The returned vector
should not be mutated.
"""
is_independent_variable(indp, sym) = is_independent_variable(symbolic_container(indp), sym)

"""
    independent_variable_symbols(indp)

Return a vector of the symbolic independent variables of the given index provider `indp`.
"""
independent_variable_symbols(indp) = independent_variable_symbols(symbolic_container(indp))

"""
    is_observed(indp, sym)

Check whether the given `sym` is an observed value in `indp`.
"""
is_observed(indp, sym) = is_observed(symbolic_container(indp), sym)

"""
    observed(indp, sym, [states])

Return the observed function of the given `sym` in `indp`. The returned function should
have the signature `(u, p) -> [values...]` where `u` and `p` is the current state and
parameter object, respectively. If `istimedependent(indp) == true`, the function should
accept the current time `t` as its third parameter. If `constant_structure(indp) == false`,
`observed` accepts a third parameter, which can either be a vector of symbols indicating
the order of states or a time index, which identifies the order of states. This function
does not need to be defined if [`is_observed`](@ref) always returns `false`. Thus,
it is mandatory to always check `is_observed` before using this function.

If `!is_markovian(indp)`, the returned function must have the signature
`(u, h, p, t) -> [values...]` where `h` is the history function, which can be called
to obtain past values of the state. The exact signature and semantics of `h` depend
on how it is used inside the returned function. `h` is obtained from a value
provider using [`get_history_function`](@ref).

See also: [`is_time_dependent`](@ref), [`is_markovian`](@ref), [`constant_structure`](@ref).
"""
observed(indp, sym) = observed(symbolic_container(indp), sym)
observed(indp, sym, states) = observed(symbolic_container(indp), sym, states)

"""
    supports_tuple_observed(indp)

Check if the given index provider supports generating observed functions for tuples of
symbolic variables. Falls back using `symbolic_container`, and returns `false` by
default.

See also: [`observed`](@ref), [`parameter_observed`](@ref), [`symbolic_container`](@ref).
"""
function supports_tuple_observed(indp)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) !== indp
        supports_tuple_observed(sc)
    else
        false
    end
end

"""
    is_time_dependent(indp)

Check if `indp` has time as (one of) its independent variables.
"""
is_time_dependent(indp) = is_time_dependent(symbolic_container(indp))

"""
    is_markovian(indp)

Check if an index provider represents a Markovian system. Markovian systems do not require
knowledge of past states to simulate. This function is only applicable to
index providers for which `is_time_dependent(indp)` returns `true`.

Non-Markovian index providers return [`observed`](@ref) functions with a different signature.
All value providers associated with a non-markovian index provider must implement
[`get_history_function`](@ref).

Returns `true` by default.
"""
function is_markovian(indp)
    if hasmethod(symbolic_container, Tuple{typeof(indp)})
        is_markovian(symbolic_container(indp))
    else
        true
    end
end

"""
    constant_structure(indp)

Check if `indp` has a constant structure. Constant structure index providers do not change
the number of variables or parameters over time.
"""
constant_structure(indp) = constant_structure(symbolic_container(indp))

"""
    all_variable_symbols(indp)

Return a vector of variable symbols in the system, including observed quantities.

For types that implement `Base.getindex` with symbolic indices using this interface,
The shorthand `sys[allvariables]` can be used as shorthand for
`valp[all_variable_symbols(indp)]`.

See: [`allvariables`](@ref).
"""
all_variable_symbols(indp) = all_variable_symbols(symbolic_container(indp))

"""
    all_symbols(indp)

Return an array of all symbols in the index provider. This includes parameters and
independent variables.
"""
all_symbols(indp) = all_symbols(symbolic_container(indp))

"""
    default_values(indp)

Return a dictionary mapping symbols in the index provider to their default value, if any.
This includes parameter symbols. The dictionary must be mutable.
"""
function default_values(indp)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        default_values(symbolic_container(indp))
    else
        Dict()
    end
end

struct SolvedVariables end

"""
    const solvedvariables = SolvedVariables()

This singleton is used as a shortcut to allow indexing of all solution variables
(excluding observed quantities). It has a [`symbolic_type`](@ref) of
[`ScalarSymbolic`](@ref). See: [`variable_symbols`](@ref).
"""
const solvedvariables = SolvedVariables()
symbolic_type(::Type{SolvedVariables}) = ScalarSymbolic()

struct AllVariables end

"""
    const allvariables = AllVariables()

This singleton is used as a shortcut to allow indexing of all solution variables
(including observed quantities). It has a [`symbolic_type`](@ref) of
[`ScalarSymbolic`](@ref). See [`all_variable_symbols`](@ref).
"""
const allvariables = AllVariables()
symbolic_type(::Type{AllVariables}) = ScalarSymbolic()
