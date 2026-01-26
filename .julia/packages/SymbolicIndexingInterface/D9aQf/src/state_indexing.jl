function set_state!(sys, val, idx)
    state_values(sys)[idx] = val
end

"""
    getsym(indp, sym)

Return a function that takes a value provider and returns the value of the symbolic
variable `sym`. If `sym` is not an observed quantity, the returned function can also
directly be called with an array of values representing the state vector. `sym` can be an
index into the state vector, a symbolic variable, a symbolic expression involving symbolic
variables in the index provider `indp`, a parameter symbol, the independent variable
symbol, or an array/tuple of the aforementioned. If the returned function is called with
a timeseries object, it can also be given a second argument representing the index at
which to return the value of `sym`.

At minimum, this requires that the value provider implement [`state_values`](@ref). To
support symbolic expressions, the value provider must implement [`observed`](@ref),
[`parameter_values`](@ref) and [`current_time`](@ref).

This function typically does not need to be implemented, and has a default implementation
relying on the above functions.

If the value provider is a parameter timeseries object, the same rules apply as
[`getp`](@ref). The difference here is that `sym` may also contain non-parameter symbols,
and the values are always returned corresponding to the state timeseries.
"""
function getsym(sys, sym)
    symtype = symbolic_type(sym)
    elsymtype = symbolic_type(eltype(sym))
    _getsym(sys, symtype, elsymtype, sym)
end

struct GetStateIndex{I} <: AbstractStateGetIndexer
    idx::I
end
function (gsi::GetStateIndex)(::Timeseries, prob)
    getindex.(state_values(prob), (gsi.idx,))
end
function (gsi::GetStateIndex)(::Timeseries, prob, i::Union{Int, CartesianIndex})
    getindex(state_values(prob, i), gsi.idx)
end
function (gsi::GetStateIndex)(::Timeseries, prob, i)
    getindex.(state_values(prob, i), gsi.idx)
end
function (gsi::GetStateIndex)(::NotTimeseries, prob)
    state_values(prob, gsi.idx)
end

function _getsym(sys, ::NotSymbolic, ::NotSymbolic, sym)
    return GetStateIndex(sym)
end

struct GetIndepvar <: AbstractStateGetIndexer end

(::GetIndepvar)(::IsTimeseriesTrait, prob) = current_time(prob)
(::GetIndepvar)(::Timeseries, prob, i) = current_time(prob, i)

struct TimeDependentObservedFunction{I, F, H} <: AbstractStateGetIndexer
    ts_idxs::I
    obsfn::F
end

function TimeDependentObservedFunction{H}(ts_idxs, obsfn) where {H}
    return TimeDependentObservedFunction{typeof(ts_idxs), typeof(obsfn), H}(ts_idxs, obsfn)
end

const NonMarkovianObservedFunction = TimeDependentObservedFunction{I, F, false} where {I, F}

indexer_timeseries_index(t::TimeDependentObservedFunction) = t.ts_idxs
function is_indexer_timeseries(::Type{G}) where {G <:
                                                 TimeDependentObservedFunction{ContinuousTimeseries}}
    return IndexerBoth()
end
function is_indexer_timeseries(::Type{G}) where {G <:
                                                 TimeDependentObservedFunction{<:Vector}}
    return IndexerMixedTimeseries()
end
function (o::TimeDependentObservedFunction)(ts::IsTimeseriesTrait, prob, args...)
    return o(ts, is_indexer_timeseries(o), prob, args...)
end

function (o::TimeDependentObservedFunction)(::Timeseries, ::IndexerBoth, prob)
    return o.obsfn.(state_values(prob),
        (parameter_values(prob),),
        current_time(prob))
end
function (o::NonMarkovianObservedFunction)(::Timeseries, ::IndexerBoth, prob)
    return o.obsfn.(state_values(prob),
        (get_history_function(prob),),
        (parameter_values(prob),),
        current_time(prob))
end
function (o::TimeDependentObservedFunction)(
        ::Timeseries, ::IndexerBoth, prob, i::Union{Int, CartesianIndex})
    return o.obsfn(state_values(prob, i), parameter_values(prob), current_time(prob, i))
end
function (o::NonMarkovianObservedFunction)(
        ::Timeseries, ::IndexerBoth, prob, i::Union{Int, CartesianIndex})
    return o.obsfn(state_values(prob, i), get_history_function(prob),
        parameter_values(prob), current_time(prob, i))
end
function (o::TimeDependentObservedFunction)(ts::Timeseries, ::IndexerBoth, prob, ::Colon)
    return o(ts, prob)
end
function (o::TimeDependentObservedFunction)(
        ts::Timeseries, ::IndexerBoth, prob, i::AbstractArray{Bool})
    map(only(to_indices(current_time(prob), (i,)))) do idx
        o(ts, prob, idx)
    end
end
function (o::TimeDependentObservedFunction)(ts::Timeseries, ::IndexerBoth, prob, i)
    o.((ts,), (prob,), i)
end
function (o::TimeDependentObservedFunction)(::NotTimeseries, ::IndexerBoth, prob)
    return o.obsfn(state_values(prob), parameter_values(prob), current_time(prob))
end
function (o::NonMarkovianObservedFunction)(::NotTimeseries, ::IndexerBoth, prob)
    return o.obsfn(state_values(prob), get_history_function(prob),
        parameter_values(prob), current_time(prob))
end

function (o::TimeDependentObservedFunction)(
        ::Timeseries, ::IndexerMixedTimeseries, prob, args...)
    throw(MixedParameterTimeseriesIndexError(prob, indexer_timeseries_index(o)))
end
function (o::TimeDependentObservedFunction)(
        ::NotTimeseries, ::IndexerMixedTimeseries, prob, args...)
    return o.obsfn(state_values(prob), parameter_values(prob), current_time(prob))
end
function (o::NonMarkovianObservedFunction)(
        ::NotTimeseries, ::IndexerMixedTimeseries, prob, args...)
    return o.obsfn(state_values(prob), get_history_function(prob),
        parameter_values(prob), current_time(prob))
end

struct TimeIndependentObservedFunction{F} <: AbstractStateGetIndexer
    obsfn::F
end

function (o::TimeIndependentObservedFunction)(::IsTimeseriesTrait, prob)
    return o.obsfn(state_values(prob), parameter_values(prob))
end

function _getsym(sys, ::ScalarSymbolic, ::SymbolicTypeTrait, sym)
    if is_variable(sys, sym)
        idx = variable_index(sys, sym)
        return getsym(sys, idx)
    elseif is_parameter(sys, sym)
        return getp(sys, sym)
    elseif is_independent_variable(sys, sym)
        return GetIndepvar()
    elseif is_observed(sys, sym)
        if !is_time_dependent(sys)
            return TimeIndependentObservedFunction(observed(sys, sym))
        end

        ts_idxs = get_all_timeseries_indexes(sys, sym)
        if ContinuousTimeseries() in ts_idxs
            if length(ts_idxs) == 1
                ts_idxs = only(ts_idxs)
            else
                ts_idxs = collect(ts_idxs)
            end
            fn = observed(sys, sym)
            return TimeDependentObservedFunction{is_markovian(sys)}(ts_idxs, fn)
        else
            return getp(sys, sym)
        end
    end
    error("Invalid symbol $sym for `getsym`")
end

struct MultipleGetters{I, G} <: AbstractStateGetIndexer
    ts_idxs::I
    getters::G
end

indexer_timeseries_index(mg::MultipleGetters) = mg.ts_idxs
function is_indexer_timeseries(::Type{G}) where {G <: MultipleGetters{ContinuousTimeseries}}
    return IndexerBoth()
end
function is_indexer_timeseries(::Type{G}) where {G <: MultipleGetters{<:Vector}}
    return IndexerMixedTimeseries()
end
function is_indexer_timeseries(::Type{G}) where {G <: MultipleGetters{Nothing}}
    return IndexerNotTimeseries()
end

function (mg::MultipleGetters)(ts::IsTimeseriesTrait, prob, args...)
    return mg(ts, is_indexer_timeseries(mg), prob, args...)
end

function (mg::MultipleGetters)(ts::Timeseries, ::IndexerBoth, prob)
    return mg.((ts,), (prob,), eachindex(current_time(prob)))
end
function (mg::MultipleGetters)(
        ::Timeseries, ::IndexerBoth, prob, i::Union{Int, CartesianIndex})
    return map(CallWith(prob, i), mg.getters)
end
function (mg::MultipleGetters)(ts::Timeseries, ::IndexerBoth, prob, ::Colon)
    return mg(ts, prob)
end
function (mg::MultipleGetters)(ts::Timeseries, ::IndexerBoth, prob, i::AbstractArray{Bool})
    return map(only(to_indices(current_time(prob), (i,)))) do idx
        mg(ts, prob, idx)
    end
end
function (mg::MultipleGetters)(ts::Timeseries, ::IndexerBoth, prob, i)
    mg.((ts,), (prob,), i)
end
function (mg::MultipleGetters)(
        ::NotTimeseries, ::Union{IndexerBoth, IndexerNotTimeseries}, prob)
    return map(g -> g(prob), mg.getters)
end

function (mg::MultipleGetters)(::Timeseries, ::IndexerMixedTimeseries, prob, args...)
    throw(MixedParameterTimeseriesIndexError(prob, indexer_timeseries_index(mg)))
end
function (mg::MultipleGetters)(::NotTimeseries, ::IndexerMixedTimeseries, prob, args...)
    return map(g -> g(prob), mg.getters)
end

struct AsTupleWrapper{N, A, G} <: AbstractStateGetIndexer
    getter::G
end

AsTupleWrapper{N}(getter::G) where {N, G} = AsTupleWrapper{N, Nothing, G}(getter)
AsTupleWrapper{N, A}(getter::G) where {N, A, G} = AsTupleWrapper{N, A, G}(getter)

wrap_tuple(::AsTupleWrapper{N, Nothing}, val) where {N} = ntuple(i -> val[i], Val(N))
function wrap_tuple(::AsTupleWrapper{N, A}, val) where {N, A}
    NamedTuple{A}(ntuple(i -> val[i], Val(N)))
end

function (atw::AsTupleWrapper)(::Timeseries, prob)
    return wrap_tuple.((atw,), atw.getter(prob))
end
function (atw::AsTupleWrapper)(::Timeseries, prob, i::Union{Int, CartesianIndex})
    return wrap_tuple(atw, atw.getter(prob, i))
end
function (atw::AsTupleWrapper)(::Timeseries, prob, i)
    return wrap_tuple.((atw,), atw.getter(prob, i))
end
function (atw::AsTupleWrapper)(::NotTimeseries, prob)
    wrap_tuple(atw, atw.getter(prob))
end

for (t1, t2) in [
    (ScalarSymbolic, Any),
    (ArraySymbolic, Any),
    (NotSymbolic, Union{<:Tuple, <:NamedTuple, <:AbstractArray})
]
    @eval function _getsym(sys, ::NotSymbolic, elt::$t1, sym::$t2)
        if isempty(sym)
            return MultipleGetters(ContinuousTimeseries(), sym)
        end
        sym_arr = sym isa Union{Tuple, NamedTuple} ? collect(sym) : sym
        supports_tuple = supports_tuple_observed(sys)
        num_observed = 0
        for s in sym
            num_observed += is_observed(sys, s)
            num_observed > 1 && break # exit early, we only need to know whether 0, 1, or more
        end
        if !is_time_dependent(sys)
            if num_observed == 0 || num_observed == 1 && sym isa Tuple
                return MultipleGetters(nothing, getsym.((sys,), sym))
            else
                obs = supports_tuple ? observed(sys, sym) : observed(sys, sym_arr)
                getter = TimeIndependentObservedFunction(obs)
                if sym isa Tuple
                    getter = AsTupleWrapper{length(sym)}(getter)
                elseif sym isa NamedTuple
                    getter = AsTupleWrapper{length(sym), keys(sym)}(getter)
                end
                return getter
            end
        end
        ts_idxs = get_all_timeseries_indexes(sys, sym_arr)
        if !(ContinuousTimeseries() in ts_idxs)
            return getp(sys, sym)
        end
        if length(ts_idxs) == 1
            ts_idxs = only(ts_idxs)
        else
            ts_idxs = collect(ts_idxs)
        end

        if num_observed == 0 || num_observed == 1 && sym isa Union{Tuple, NamedTuple}
            _sym = sym isa NamedTuple ? Tuple(sym) : sym
            getters = getsym.((sys,), _sym)
            getter = MultipleGetters(ts_idxs, getters)
            if sym isa NamedTuple
                getter = AsTupleWrapper{length(sym), keys(sym)}(getter)
            end
            return getter
        else
            obs = supports_tuple ? observed(sys, sym) : observed(sys, sym_arr)
            getter = if is_time_dependent(sys)
                TimeDependentObservedFunction{is_markovian(sys)}(ts_idxs, obs)
            else
                TimeIndependentObservedFunction(obs)
            end
            if sym isa Tuple && !supports_tuple
                getter = AsTupleWrapper{length(sym)}(getter)
            elseif sym isa NamedTuple
                getter = AsTupleWrapper{length(sym), keys(sym)}(getter)
            end
            return getter
        end
    end
end

function _getsym(sys, ::ArraySymbolic, ::SymbolicTypeTrait, sym)
    if is_variable(sys, sym)
        idx = variable_index(sys, sym)
        return getsym(sys, idx)
    elseif is_parameter(sys, sym)
        return getp(sys, sym)
    end
    return getsym(sys, collect(sym))
end

# setsym doesn't need the same `let` blocks to be inferred for some reason

"""
    setsym(sys, sym)

Return a function that takes a value provider and a value, and sets the the state `sym` to
that value. Note that `sym` can be an index, a symbolic variable, or an array/tuple of the
aforementioned.

Requires that the value provider implement [`state_values`](@ref) and the returned
collection be a mutable reference to the state vector in the value provider. Alternatively,
if this is not possible or additional actions need to be performed when updating state,
[`set_state!`](@ref) can be defined. This function does not work on types for which
[`is_timeseries`](@ref) is [`Timeseries`](@ref).
"""
function setsym(sys, sym)
    symtype = symbolic_type(sym)
    elsymtype = symbolic_type(eltype(sym))
    _setsym(sys, symtype, elsymtype, sym)
end

struct SetStateIndex{I} <: AbstractSetIndexer
    idx::I
end

function (ssi::SetStateIndex)(prob, val)
    set_state!(prob, val, ssi.idx)
end

function _setsym(sys, ::NotSymbolic, ::NotSymbolic, sym)
    return SetStateIndex(sym)
end

function _setsym(sys, ::ScalarSymbolic, ::SymbolicTypeTrait, sym)
    if is_variable(sys, sym)
        idx = variable_index(sys, sym)
        return SetStateIndex(idx)
    elseif is_parameter(sys, sym)
        return setp(sys, sym)
    end
    error("Invalid symbol $sym for `setsym`")
end

struct NamedTupleSetter{S <: NamedTuple} <: AbstractSetIndexer
    setter::S
end

function (nts::NamedTupleSetter)(prob, val)
    _generated_setter(nts, prob, val)
end

@generated function _generated_setter(
        nts::NamedTupleSetter{<:NamedTuple{N1}}, prob, val::NamedTuple{N2}) where {N1, N2}
    expr = Expr(:block)
    for (i, name) in enumerate(N2)
        idx = findfirst(isequal(name), N1)
        if idx === nothing
            throw(ArgumentError("""
            Invalid name $(name) in value. Must be one of $(N1).
            """))
        end
        push!(expr.args, :(nts.setter[$idx](prob, val[$i])))
    end
    return expr
end

for (t1, t2) in [
    (ScalarSymbolic, Any),
    (ArraySymbolic, Any),
    (NotSymbolic, Union{<:Tuple, <:NamedTuple, <:AbstractArray})
]
    @eval function _setsym(sys, ::NotSymbolic, ::$t1, sym::$t2)
        if sym isa NamedTuple
            setters = NamedTuple{keys(sym)}(setsym.((sys,), values(sym)))
            return NamedTupleSetter(setters)
        end
        setters = setsym.((sys,), sym)
        return MultipleSetters(setters)
    end
end

function _setsym(sys, ::ArraySymbolic, ::SymbolicTypeTrait, sym)
    if is_variable(sys, sym)
        idx = variable_index(sys, sym)
        if idx isa AbstractArray
            return MultipleSetters(SetStateIndex.(idx))
        else
            return SetStateIndex(idx)
        end
    elseif is_parameter(sys, sym)
        return setp(sys, sym)
    end
    return setsym(sys, collect(sym))
end

const getu = getsym
const setu = setsym

"""
    setsym_oop(indp, sym)

Return a function which takes a value provider `valp` and a value `val`, and returns
`state_values(valp), parameter_values(valp)` with the states/parameters in `sym` set to the
corresponding values in `val`. This allows changing the types of values stored, and leverages
[`remake_buffer`](@ref). Note that `sym` can be an index, a symbolic variable, or an
array/tuple of the aforementioned. All entries `s` in `sym` must satisfy `is_variable(indp, s)`
or `is_parameter(indp, s)`.

Requires that the value provider implement `state_values`, `parameter_values` and `remake_buffer`.
"""
function setsym_oop(indp, sym)
    symtype = symbolic_type(sym)
    elsymtype = symbolic_type(eltype(sym))
    return _setsym_oop(indp, symtype, elsymtype, sym)
end

struct FullSetter{S, P, I, J}
    state_setter::S
    param_setter::P
    state_split::I
    param_split::J
end

FullSetter(ssetter, psetter) = FullSetter(ssetter, psetter, nothing, nothing)

function (fs::FullSetter)(valp, val)
    check_both_state_and_parameter_provider(valp)

    return fs.state_setter(valp, val[fs.state_split]),
    fs.param_setter(valp, val[fs.param_split])
end

function (fs::FullSetter{Nothing})(valp, val)
    check_both_state_and_parameter_provider(valp)

    return state_values(valp), fs.param_setter(valp, val)
end

function (fs::(FullSetter{S, Nothing} where {S}))(valp, val)
    check_both_state_and_parameter_provider(valp)

    return fs.state_setter(valp, val), parameter_values(valp)
end

function (fs::(FullSetter{Nothing, Nothing}))(valp, val)
    check_both_state_and_parameter_provider(valp)

    return state_values(valp), parameter_values(valp)
end

function _setsym_oop(indp, ::NotSymbolic, ::NotSymbolic, sym)
    return FullSetter(OOPSetter(_root_indp(indp), sym, true), nothing)
end

function _setsym_oop(indp, ::ScalarSymbolic, ::SymbolicTypeTrait, sym)
    if (idx = variable_index(indp, sym)) !== nothing
        return FullSetter(OOPSetter(_root_indp(indp), idx, true), nothing)
    elseif (idx = parameter_index(indp, sym)) !== nothing
        return FullSetter(nothing, OOPSetter(_root_indp(indp), idx, false))
    end
    throw(NotVariableOrParameter("setsym_oop", sym))
end

for (t1, t2) in [
    (ScalarSymbolic, Any),
    (NotSymbolic, Union{<:Tuple, <:AbstractArray})
]
    @eval function _setsym_oop(indp, ::NotSymbolic, ::$t1, sym::$t2)
        vars = []
        state_split = eltype(eachindex(sym))[]
        pars = []
        param_split = eltype(eachindex(sym))[]
        for (i, s) in enumerate(sym)
            if (idx = variable_index(indp, s)) !== nothing
                push!(vars, idx)
                push!(state_split, i)
            elseif (idx = parameter_index(indp, s)) !== nothing
                push!(pars, idx)
                push!(param_split, i)
            else
                throw(NotVariableOrParameter("setsym_oop", s))
            end
        end
        if sym isa Tuple
            vars = Tuple(vars)
            pars = Tuple(pars)
        end
        indp = _root_indp(indp)
        return FullSetter(isempty(vars) ? nothing : OOPSetter(indp, identity.(vars), true),
            isempty(pars) ? nothing : OOPSetter(indp, identity.(pars), false),
            state_split, param_split)
    end
end

function _setsym_oop(indp, ::ArraySymbolic, ::SymbolicTypeTrait, sym)
    if (idx = variable_index(indp, sym)) !== nothing
        return setsym_oop(indp, idx)
    elseif (idx = parameter_index(indp, sym)) !== nothing
        return FullSetter(
            nothing, OOPSetter(indp, idx, false))
    end
    return setsym_oop(indp, collect(sym))
end
