###########
# Parameter Indexing
###########

"""
    parameter_values(valp)
    parameter_values(valp, i)

Return an indexable collection containing the value of each parameter in `valp`. The two-
argument version of this function returns the parameter value at index `i`. The
two-argument version of this function will default to returning
`parameter_values(valp)[i]`.

If this function is called with an `AbstractArray` or `Tuple`, it will return the same
array/tuple.
"""
function parameter_values end

parameter_values(arr::AbstractArray) = arr
parameter_values(arr::Tuple) = arr
parameter_values(arr::AbstractArray, i) = arr[i]
parameter_values(arr::Tuple, i) = arr[i]
parameter_values(prob, i) = parameter_values(parameter_values(prob), i)

"""
    get_parameter_timeseries_collection(valp)

Return the [`ParameterTimeseriesCollection`](@ref) contained in timeseries value provider
`valp`. Valid only for value providers where [`is_parameter_timeseries`](@ref) returns
[`Timeseries`](@ref).
"""
function get_parameter_timeseries_collection end

"""
    with_updated_parameter_timeseries_values(indp, params, args::Pair...)

Return an indexable collection containing the value of all parameters in `params`, with
parameters belonging to specific timeseries updated to different values. Each element in
`args...` contains the timeseries index as the first value, and the saved parameter values
in that partition. Not all parameter timeseries have to be updated using this method. If
an in-place update can be performed, it should be done and the modified `params` returned.
This method falls back on the basis of `symbolic_container(indp)`.

Note that here `params` is the parameter object.
"""
function with_updated_parameter_timeseries_values(indp, params, args...)
    return with_updated_parameter_timeseries_values(
        symbolic_container(indp), params, args...)
end

"""
    set_parameter!(valp, val, idx)

Set the parameter at index `idx` to `val` for value provider `valp`. This defaults to
modifying `parameter_values(valp)`. If any additional bookkeeping needs to be performed
or the default implementation does not work for a particular type, this method needs to
be defined to enable the proper functioning of [`setp`](@ref).

See: [`parameter_values`](@ref)
"""
function set_parameter! end

# Tuple only included for the error message
function set_parameter!(sys::Union{AbstractArray, Tuple}, val, idx)
    sys[idx] = val
end
set_parameter!(sys, val, idx) = set_parameter!(parameter_values(sys), val, idx)

"""
    finalize_parameters_hook!(valp, sym)

This is a callback run one for each call to the function returned by [`setp`](@ref)
which can be used to update internal data structures when parameters are modified.
This is in contrast to [`set_parameter!`](@ref) which is run once for each parameter
that is updated.
"""
finalize_parameters_hook!(valp, sym) = nothing

###########
# State Indexing
###########

"""
    state_values(valp)
    state_values(valp, i)

Return an indexable collection containing the values of all states in the value provider
`p`. If `is_timeseries(valp)` is [`Timeseries`](@ref), return a vector of arrays,
each of which contain the state values at the corresponding timestep. In this case, the
two-argument version of the function can also be implemented to efficiently return
the state values at timestep `i`. By default, the two-argument method calls
`state_values(valp)[i]`. If `i` consists of multiple indices (for example, `Colon`,
`AbstractArray{Int}`, `AbstractArray{Bool}`) specialized methods may be defined for
efficiency. By default, `state_values(valp, ::Colon) = state_values(valp)` to avoid
copying the timeseries.

If this function is called with an `AbstractArray`, it will return the same array.

See: [`is_timeseries`](@ref)
"""
function state_values end

state_values(arr::AbstractArray) = arr
state_values(arr, i) = state_values(arr)[i]
state_values(arr, ::Colon) = state_values(arr)

"""
    set_state!(valp, val, idx)

Set the state at index `idx` to `val` for value provider `valp`. This defaults to modifying
`state_values(valp)`. If any additional bookkeeping needs to be performed or the
default implementation does not work for a particular type, this method needs to be
defined to enable the proper functioning of [`setsym`](@ref).

See: [`state_values`](@ref)
"""
function set_state! end

"""
    current_time(valp)
    current_time(valp, i)

Return the current time in the value provider `valp`. If
`is_timeseries(valp)` is [`Timeseries`](@ref), return the vector of timesteps at which
the state value is saved. In this case, the two-argument version of the function can
also be implemented to efficiently return the time at timestep `i`. By default, the two-
argument method calls `current_time(p)[i]`. It is assumed that the timeseries is sorted
in increasing order.

If `i` consists of multiple indices (for example, `Colon`, `AbstractArray{Int}`,
`AbstractArray{Bool}`) specialized methods may be defined for efficiency. By default,
`current_time(valp, ::Colon) = current_time(valp)` to avoid copying the timeseries.

By default, the single-argument version acts as the identity function if
`valp isa AbstractVector`.

See: [`is_timeseries`](@ref)
"""
function current_time end

current_time(arr::AbstractVector) = arr
current_time(valp, i) = current_time(valp)[i]
current_time(valp, ::Colon) = current_time(valp)

"""
    get_history_function(valp)

Return the history function for a value provider. This is required for all value providers
associated with an index provider `indp` for which `!is_markovian(indp)`.

See also: [`is_markovian`](@ref).
"""
function get_history_function end

###########
# Utilities
###########

abstract type AbstractIndexer end

abstract type AbstractGetIndexer <: AbstractIndexer end
abstract type AbstractStateGetIndexer <: AbstractGetIndexer end
abstract type AbstractParameterGetIndexer <: AbstractGetIndexer end
abstract type AbstractSetIndexer <: AbstractIndexer end

(ai::AbstractStateGetIndexer)(prob) = ai(is_timeseries(prob), prob)
(ai::AbstractStateGetIndexer)(prob, i) = ai(is_timeseries(prob), prob, i)
(ai::AbstractParameterGetIndexer)(prob) = ai(is_parameter_timeseries(prob), prob)
(ai::AbstractParameterGetIndexer)(prob, i) = ai(is_parameter_timeseries(prob), prob, i)
function (ai::AbstractParameterGetIndexer)(buffer::AbstractArray, prob)
    ai(buffer, is_parameter_timeseries(prob), prob)
end
function (ai::AbstractParameterGetIndexer)(buffer::AbstractArray, prob, i)
    ai(buffer, is_parameter_timeseries(prob), prob, i)
end

abstract type IndexerTimeseriesType end

# Can only index parameter timeseries
struct IndexerOnlyTimeseries <: IndexerTimeseriesType end
# Can only index non-timeseries objects
struct IndexerMixedTimeseries <: IndexerTimeseriesType end
# Can index timeseres and non-timeseries objects, has the same value at all times
struct IndexerNotTimeseries <: IndexerTimeseriesType end
# Value changes over time, can index timeseries and non-timeseries objects
struct IndexerBoth <: IndexerTimeseriesType end

is_indexer_timeseries(x) = is_indexer_timeseries(typeof(x))
function indexer_timeseries_index end

const AtLeastTimeseriesIndexer = Union{IndexerOnlyTimeseries, IndexerBoth}
const AtLeastNotTimeseriesIndexer = Union{IndexerNotTimeseries, IndexerBoth}

as_timeseries_indexer(x) = as_timeseries_indexer(is_indexer_timeseries(x), x)
as_timeseries_indexer(::IndexerOnlyTimeseries, x) = x
as_timeseries_indexer(::IndexerNotTimeseries, x) = x
as_not_timeseries_indexer(x) = as_not_timeseries_indexer(is_indexer_timeseries(x), x)
as_not_timeseries_indexer(::IndexerNotTimeseries, x) = x

function _postprocess_tsidxs(ts_idxs)
    delete!(ts_idxs, ContinuousTimeseries())
    if isempty(ts_idxs)
        return nothing
    elseif length(ts_idxs) == 1
        return only(ts_idxs)
    else
        return collect(ts_idxs)
    end
end

struct CallWith{A}
    args::A

    CallWith(args...) = new{typeof(args)}(args)
end

function (cw::CallWith)(arg)
    arg(cw.args...)
end

function _call(f, args...)
    return f(args...)
end

struct Fix1Multiple{F, A}
    f::F
    arg::A
end

function (fn::Fix1Multiple)(args...)
    fn.f(fn.arg, args...)
end

struct OOPSetter{S, I, D}
    indp::I
    idxs::D
end

OOPSetter(indp, idxs, isstate) = OOPSetter{isstate, typeof(indp), typeof(idxs)}(indp, idxs)

function (os::OOPSetter{true})(valp, val)
    buffer = hasmethod(state_values, Tuple{typeof(valp)}) ? state_values(valp) : valp
    return remake_buffer(os.indp, buffer, (os.idxs,), (val,))
end

function (os::OOPSetter{false})(valp, val)
    buffer = hasmethod(parameter_values, Tuple{typeof(valp)}) ? parameter_values(valp) :
             valp
    return remake_buffer(os.indp, buffer, (os.idxs,), (val,))
end

function (os::OOPSetter{true})(valp, val::Union{Tuple, AbstractArray})
    buffer = hasmethod(state_values, Tuple{typeof(valp)}) ? state_values(valp) : valp
    if os.idxs isa Union{Tuple, AbstractArray}
        return remake_buffer(os.indp, buffer, os.idxs, val)
    else
        return remake_buffer(os.indp, buffer, (os.idxs,), (val,))
    end
end

function (os::OOPSetter{false})(valp, val::Union{Tuple, AbstractArray})
    buffer = hasmethod(parameter_values, Tuple{typeof(valp)}) ? parameter_values(valp) :
             valp
    if os.idxs isa Union{Tuple, AbstractArray}
        return remake_buffer(os.indp, buffer, os.idxs, val)
    else
        return remake_buffer(os.indp, buffer, (os.idxs,), (val,))
    end
end

function _root_indp(indp)
    if hasmethod(symbolic_container, Tuple{typeof(indp)}) &&
       (sc = symbolic_container(indp)) != indp
        return _root_indp(sc)
    else
        return indp
    end
end

###########
# Errors
###########

struct ParameterTimeseriesValueIndexMismatchError{P <: IsTimeseriesTrait} <: Exception
    valp::Any
    indexer::Any
    args::Any

    function ParameterTimeseriesValueIndexMismatchError{Timeseries}(valp, indexer, args)
        if is_parameter_timeseries(valp) != Timeseries()
            throw(ArgumentError("""
                This should never happen. Expected parameter timeseries value provider, \
                got $(valp). Open an issue in SymbolicIndexingInterface.jl with an MWE.
            """))
        end
        return new{Timeseries}(valp, indexer, args)
    end
    function ParameterTimeseriesValueIndexMismatchError{NotTimeseries}(valp, indexer)
        if is_parameter_timeseries(valp) != NotTimeseries()
            throw(ArgumentError("""
                This should never happen. Expected non-parameter timeseries value \
                provider, got $(valp). Open an issue in SymbolicIndexingInterface.jl \
                with an MWE.
            """))
        end
        if is_indexer_timeseries(indexer) != IndexerOnlyTimeseries()
            throw(ArgumentError("""
                This should never happen. Expected timeseries indexer, got $(indexer). \
                Open an issue in SymbolicIndexingInterface.jl with an MWE.
            """))
        end
        return new{NotTimeseries}(valp, indexer, nothing)
    end
end

function Base.showerror(io::IO, err::ParameterTimeseriesValueIndexMismatchError{Timeseries})
    print(io, """
        Invalid indexing operation: tried to access object of type $(typeof(err.valp)) \
        (which is a parameter timeseries object) with non-timeseries indexer \
        $(err.indexer) at index $(err.args) in the timeseries.
    """)
end

function Base.showerror(
        io::IO, err::ParameterTimeseriesValueIndexMismatchError{NotTimeseries})
    print(io, """
        Invalid indexing operation: tried to access object of type $(typeof(err.valp)) \
        (which is not a parameter timeseries object) using timeseries indexer \
        $(err.indexer).
    """)
end

struct MixedParameterTimeseriesIndexError <: Exception
    obj::Any
    ts_idxs::Any
end

function Base.showerror(io::IO, err::MixedParameterTimeseriesIndexError)
    print(io, """
        Invalid indexing operation: tried to access object of type $(typeof(err.obj)) \
        (which is a parameter timeseries object) with variables having mixed timeseries \
        indexes $(err.ts_idxs).
    """)
end

struct NotVariableOrParameter <: Exception
    fn::Any
    sym::Any
end

function Base.showerror(io::IO, err::NotVariableOrParameter)
    print(
        io, """
      `$(err.fn)` requires that the symbolic variable(s) passed to it satisfy `is_variable`
      or `is_parameter`. Got `$(err.sym)` which is neither.
  """)
end

function MustBeBothStateAndParameterProviderError(missing_state::Bool)
    ArgumentError("""
        A setter returned from `setsym_oop` must be called with a value provider that \
        contains both states and parameters. The given value provided does not \
        implement `$(missing_state ? "state_values" : "parameter_values")`.
        """)
end

function check_both_state_and_parameter_provider(valp)
    if !hasmethod(state_values, Tuple{typeof(valp)}) || state_values(valp) === valp
        throw(MustBeBothStateAndParameterProviderError(true))
    end
    if !hasmethod(parameter_values, Tuple{typeof(valp)}) || parameter_values(valp) === valp
        throw(MustBeBothStateAndParameterProviderError(false))
    end
end
