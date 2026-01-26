"""
    struct BatchedInterface{S <: AbstractVector, I}
    function BatchedInterface(indp_syms::Tuple...)

A struct which stores information for batched calls to [`getsym`](@ref) or [`setsym`](@ref).
Given `Tuple`s, where the first element of each tuple is an index provider and the second
an array of symbolic variables (either states or parameters) in the index provider,
`BatchedInterface` will compute the union of all symbols and associate each symbol with
the first index provider with which it occurs.

For example, given two index providers `s1 = SymbolCache([:x, :y, :z])` and
`s2 = SymbolCache([:y, :z, :w])`, `BatchedInterface((s1, [:x, :y]), (s2, [:y, :z]))` will
associate `:x` and `:y` with `s1` and `:z` with `s2`. The information that `s1` had
associated symbols `:x` and `:y` and `s2` had associated symbols `:y` and `:z` will also
be retained internally.

`BatchedInterface` implements [`variable_symbols`](@ref), [`is_variable`](@ref),
[`variable_index`](@ref) to query the order of symbols in the union.

See [`getsym`](@ref) and [`setsym`](@ref) for further details.

See also: [`associated_systems`](@ref).
"""
struct BatchedInterface{S <: AbstractVector, I, T, P}
    "Order of symbols in the union."
    symbol_order::S
    "Index of the index provider each symbol in the union is associated with."
    associated_systems::Vector{Int}
    "Index of symbol in the index provider it is associated with."
    associated_indexes::I
    "Whether the symbol is a state in the index provider it is associated with."
    isstate::BitVector
    "Map from index provider to indexes of its symbols in the union."
    system_to_symbol_subset::Vector{Vector{Int}}
    "Map from index provider to indexes of its symbols in the index provider."
    system_to_symbol_indexes::Vector{Vector{T}}
    "Map from index provider to whether each of its symbols is a state in the index provider."
    system_to_isstate::Vector{BitVector}
    "Index providers, in order"
    index_providers::Vector{P}
end

function BatchedInterface(syssyms::Tuple...)
    symbol_order = []
    associated_systems = Int[]
    associated_indexes = []
    isstate = BitVector()
    system_to_symbol_subset = Vector{Int}[]
    system_to_symbol_indexes = []
    system_to_isstate = BitVector[]
    index_providers = []
    for (i, (sys, syms)) in enumerate(syssyms)
        symbol_subset = Int[]
        symbol_indexes = []
        system_isstate = BitVector()
        allsyms = []
        root_indp = sys
        while applicable(symbolic_container, root_indp) &&
            (sc = symbolic_container(root_indp)) != root_indp
            root_indp = sc
        end
        push!(index_providers, root_indp)
        for sym in syms
            if symbolic_type(sym) === NotSymbolic()
                error("Only symbolic variables allowed in BatchedInterface.")
            end
            if symbolic_type(sym) === ArraySymbolic()
                append!(allsyms, collect(sym))
            else
                push!(allsyms, sym)
            end
        end
        for sym in allsyms
            if !is_variable(sys, sym) && !is_parameter(sys, sym)
                error("Only variables and parameters allowed in BatchedInterface.")
            end
            if !any(isequal(sym), symbol_order)
                push!(symbol_order, sym)
                push!(associated_systems, i)
                push!(isstate, is_variable(sys, sym))
                if isstate[end]
                    push!(associated_indexes, variable_index(sys, sym))
                else
                    push!(associated_indexes, parameter_index(sys, sym))
                end
            end
            push!(symbol_subset, findfirst(isequal(sym), symbol_order))
            push!(system_isstate, is_variable(sys, sym))
            push!(symbol_indexes,
                system_isstate[end] ? variable_index(sys, sym) : parameter_index(sys, sym))
        end
        push!(system_to_symbol_subset, symbol_subset)
        push!(system_to_symbol_indexes, identity.(symbol_indexes))
        push!(system_to_isstate, system_isstate)
    end
    symbol_order = identity.(symbol_order)
    associated_indexes = identity.(associated_indexes)
    system_to_symbol_indexes = identity.(system_to_symbol_indexes)

    return BatchedInterface{typeof(symbol_order), typeof(associated_indexes),
        eltype(eltype(system_to_symbol_indexes)), eltype(index_providers)}(
        symbol_order, associated_systems, associated_indexes, isstate,
        system_to_symbol_subset, system_to_symbol_indexes, system_to_isstate,
        identity.(index_providers))
end

variable_symbols(bi::BatchedInterface) = bi.symbol_order
variable_index(bi::BatchedInterface, sym) = findfirst(isequal(sym), bi.symbol_order)
is_variable(bi::BatchedInterface, sym) = variable_index(bi, sym) !== nothing

"""
    associated_systems(bi::BatchedInterface)

Return an array of integers of the same length as `variable_symbols(bi)` where each value
is the index of the index provider associated with the corresponding symbol in
`variable_symbols(bi)`.
"""
associated_systems(bi::BatchedInterface) = bi.associated_systems

"""
    getsym(bi::BatchedInterface)

Given a [`BatchedInterface`](@ref) composed from `n` index providers (and corresponding
symbols), return a function which takes `n` corresponding value providers and returns an
array of the values of the symbols in the union. The returned function can also be passed
an `AbstractArray` of the appropriate `eltype` and size as its first argument, in which
case the operation will populate the array in-place with the values of the symbols in the
union.

Note that all of the value providers passed to the function returned by `getsym` must satisfy
`is_timeseries(prob) === NotTimeseries()`.

The value of the `i`th symbol in the union (obtained through `variable_symbols(bi)[i]`) is
obtained from the problem corresponding to the associated index provider (i.e. the value
provider at index `associated_systems(bi)[i]`).

See also: [`variable_symbols`](@ref), [`associated_systems`](@ref), [`is_timeseries`](@ref),
[`NotTimeseries`](@ref).
"""
function getsym(bi::BatchedInterface)
    numprobs = length(bi.system_to_symbol_subset)
    probnames = [Symbol(:prob, i) for i in 1:numprobs]

    fnbody = quote end
    for (i, (prob, idx, isstate)) in enumerate(zip(
        bi.associated_systems, bi.associated_indexes, bi.isstate))
        symname = Symbol(:sym, i)
        getter = isstate ? state_values : parameter_values
        probname = probnames[prob]
        push!(fnbody.args, :($symname = $getter($probname, $idx)))
    end

    oop_expr = Expr(:vect)
    for i in 1:length(bi.symbol_order)
        push!(oop_expr.args, Symbol(:sym, i))
    end

    iip_expr = quote end
    for i in 1:length(bi.symbol_order)
        symname = Symbol(:sym, i)
        push!(iip_expr.args, :(out[$i] = $symname))
    end

    oopfn = Expr(
        :function,
        Expr(:tuple, probnames...),
        quote
            $fnbody
            $oop_expr
        end
    )
    iipfn = Expr(
        :function,
        Expr(:tuple, :out, probnames...),
        quote
            $fnbody
            $iip_expr
            out
        end
    )

    return let oop = @RuntimeGeneratedFunction(oopfn),
        iip = @RuntimeGeneratedFunction(iipfn)

        _getter(probs...) = oop(probs...)
        _getter(out::AbstractArray, probs...) = iip(out, probs...)
        _getter
    end
end

"""
    setsym(bi::BatchedInterface)

Given a [`BatchedInterface`](@ref) composed from `n` index providers (and corresponding
symbols), return a function which takes `n` corresponding problems and an array of the
values, and updates each of the problems with the values of the corresponding symbols.

Note that all of the value providers passed to the function returned by `setsym` must satisfy
`is_timeseries(prob) === NotTimeseries()`.

Note that if any subset of the `n` index providers share common symbols (among those passed
to `BatchedInterface`) then all of the corresponding value providers in the subset will be
updated with the values of the common symbols.

See also: [`is_timeseries`](@ref), [`NotTimeseries`](@ref).
"""
function setsym(bi::BatchedInterface)
    numprobs = length(bi.system_to_symbol_subset)
    probnames = [Symbol(:prob, i) for i in 1:numprobs]

    full_update_fnexpr = let fnbody = quote end
        for (sys_idx, subset) in enumerate(bi.system_to_symbol_subset)
            probname = probnames[sys_idx]
            for (idx_in_subset, idx_in_union) in enumerate(subset)
                idx = bi.system_to_symbol_indexes[sys_idx][idx_in_subset]
                isstate = bi.system_to_isstate[sys_idx][idx_in_subset]
                setter = isstate ? set_state! : set_parameter!
                push!(fnbody.args, :($setter($probname, vals[$idx_in_union], $idx)))
            end
            # also run hook
            if !all(bi.system_to_isstate[sys_idx])
                paramidxs = [bi.system_to_symbol_indexes[sys_idx][idx_in_subset]
                             for idx_in_subset in 1:length(subset)
                             if !bi.system_to_isstate[sys_idx][idx_in_subset]]
                push!(fnbody.args, :($finalize_parameters_hook!($probname, $paramidxs)))
            end
        end
        push!(fnbody.args, :(return vals))
        Expr(
            :function,
            Expr(:tuple, probnames..., :vals),
            fnbody
        )
    end

    partial_update_fnexpr = let fnbody = quote end
        curfnbody = fnbody
        for (sys_idx, subset) in enumerate(bi.system_to_symbol_subset)
            newcurfnbody = if sys_idx == 1
                Expr(:if, :(idx == $sys_idx))
            else
                Expr(:elseif, :(idx == $sys_idx))
            end
            push!(curfnbody.args, newcurfnbody)
            curfnbody = newcurfnbody

            ifbody = quote end
            push!(curfnbody.args, ifbody)

            probname = :prob
            for (idx_in_subset, idx_in_union) in enumerate(subset)
                idx = bi.system_to_symbol_indexes[sys_idx][idx_in_subset]
                isstate = bi.system_to_isstate[sys_idx][idx_in_subset]
                setter = isstate ? set_state! : set_parameter!
                push!(ifbody.args, :($setter($probname, vals[$idx_in_union], $idx)))
            end
            # also run hook
            if !all(bi.system_to_isstate[sys_idx])
                paramidxs = [bi.system_to_symbol_indexes[sys_idx][idx_in_subset]
                             for idx_in_subset in 1:length(subset)
                             if !bi.system_to_isstate[sys_idx][idx_in_subset]]
                push!(ifbody.args, :($finalize_parameters_hook!($probname, $paramidxs)))
            end
        end
        push!(curfnbody.args, :(error("Invalid problem index $idx")))
        push!(fnbody.args, :(return nothing))
        Expr(
            :function,
            Expr(:tuple, :prob, :idx, :vals),
            fnbody
        )
    end
    return let full_update = @RuntimeGeneratedFunction(full_update_fnexpr),
        partial_update = @RuntimeGeneratedFunction(partial_update_fnexpr)

        setter!(args...) = full_update(args...)
        setter!(prob, idx::Int, vals::AbstractVector) = partial_update(prob, idx, vals)
        setter!
    end
end

"""
    setsym_oop(bi::BatchedInterface)

Given a [`BatchedInterface`](@ref) composed from `n` index providers (and corresponding
symbols), return a function which takes `n` corresponding value providers and an array of
values, and returns an `n`-tuple where each element is a 2-tuple consisting of the updated
state values and parameter values of the corresponding value provider. Requires that the
value provider implement [`state_values`](@ref), [`parameter_values`](@ref). The updates are
performed out-of-place using [`remake_buffer`](@ref).

Note that all of the value providers passed to the returned function must satisfy
`is_timeseries(prob) === NotTimeseries()`.

Note that if any subset of the `n` index providers share common symbols (among those passed
to `BatchedInterface`) then all of the corresponding value providers in the subset will be
updated with the values of the common symbols.

See also: [`is_timeseries`](@ref), [`NotTimeseries`](@ref).
"""
function setsym_oop(bi::BatchedInterface)
    numprobs = length(bi.system_to_symbol_subset)
    probnames = [Symbol(:prob, i) for i in 1:numprobs]
    arg = :vals
    full_update = Expr(:block)

    function get_update_expr(prob::Symbol, sys_i::Int)
        union_idxs = bi.system_to_symbol_subset[sys_i]
        indp_idxs = bi.system_to_symbol_indexes[sys_i]
        isstate = bi.system_to_isstate[sys_i]
        indp = bi.index_providers[sys_i]
        curexpr = Expr(:block)

        statessym = Symbol(:states_, sys_i)
        if all(.!isstate)
            push!(curexpr.args, :($statessym = $state_values($prob)))
        else
            state_idxssym = Symbol(:state_idxs_, sys_i)
            state_idxs = indp_idxs[isstate]
            state_valssym = Symbol(:state_vals_, sys_i)
            vals_idxs = union_idxs[isstate]
            push!(curexpr.args, :($state_idxssym = $state_idxs))
            push!(curexpr.args, :($state_valssym = $view($arg, $vals_idxs)))
            push!(curexpr.args,
                :($statessym = $remake_buffer(
                    syss[$sys_i], $state_values($prob), $state_idxssym, $state_valssym)))
        end

        paramssym = Symbol(:params_, sys_i)
        if all(isstate)
            push!(curexpr.args, :($paramssym = $parameter_values($prob)))
        else
            param_idxssym = Symbol(:param_idxs_, sys_i)
            param_idxs = indp_idxs[.!isstate]
            param_valssym = Symbol(:param_vals, sys_i)
            vals_idxs = union_idxs[.!isstate]
            push!(curexpr.args, :($param_idxssym = $param_idxs))
            push!(curexpr.args, :($param_valssym = $view($arg, $vals_idxs)))
            push!(curexpr.args,
                :($paramssym = $remake_buffer(
                    syss[$sys_i], $parameter_values($prob), $param_idxssym, $param_valssym)))
        end

        return curexpr, statessym, paramssym
    end

    full_update_expr = Expr(:block)
    full_update_retval = Expr(:tuple)
    partial_update_expr = Expr(:block)
    cur_partial_update_expr = partial_update_expr
    for i in 1:numprobs
        update_expr, statesym, paramsym = get_update_expr(probnames[i], i)
        push!(full_update_expr.args, update_expr)
        push!(full_update_retval.args, Expr(:tuple, statesym, paramsym))

        cur_ifexpr = Expr(i == 1 ? :if : :elseif, :(idx == $i))
        update_expr, statesym, paramsym = get_update_expr(:prob, i)
        push!(update_expr.args, :(return ($statesym, $paramsym)))
        push!(cur_ifexpr.args, update_expr)
        push!(cur_partial_update_expr.args, cur_ifexpr)
        cur_partial_update_expr = cur_ifexpr
    end
    push!(full_update_expr.args, :(return $full_update_retval))
    push!(cur_partial_update_expr.args, :(error("Invalid problem index $idx")))

    full_update_fnexpr = Expr(
        :function, Expr(:tuple, :syss, probnames..., arg), full_update_expr)
    partial_update_fnexpr = Expr(
        :function, Expr(:tuple, :syss, :prob, :idx, arg), partial_update_expr)

    return let full_update = @RuntimeGeneratedFunction(full_update_fnexpr),
        partial_update = @RuntimeGeneratedFunction(partial_update_fnexpr),
        syss = Tuple(bi.index_providers)

        setter(args...) = full_update(syss, args...)
        setter(prob, idx::Int, vals::AbstractVector) = partial_update(syss, prob, idx, vals)
        setter
    end
end
