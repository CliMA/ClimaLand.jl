"""
    struct ParameterIndexingProxy

This struct wraps any struct implementing the value provider and index provider interfaces.
It allows `getindex` and `setindex!` operations to get/set parameter values. Requires that
the wrapped type support [`getp`](@ref) and [`setp`](@ref) for getting and setting
parameter values respectively.
"""
struct ParameterIndexingProxy{T}
    wrapped::T
end

function Base.getindex(p::ParameterIndexingProxy, idx, args...)
    getp(p.wrapped, idx)(p.wrapped, args...)
end

function Base.setindex!(p::ParameterIndexingProxy, val, idx)
    return setp(p.wrapped, idx)(p.wrapped, val)
end

function Base.show(io::IO, ::MIME"text/plain", pip::ParameterIndexingProxy)
    show_params(io, pip; num_rows = 20, show_all = false, scalarize = true)
end

"""
    show_params(io::IO, pip::ParameterIndexingProxy; num_rows = 20, show_all = false, scalarize = true, kwargs...)

Method for customizing the table output. Keyword args:
- num_rows
- show_all: whether to show all parameters. Overrides `num_rows`.
- scalarize: whether to scalarize array symbolics in the table output.
- kwargs... are passed to the pretty_table call (if PrettyTables is loaded).
"""
show_params(io, pip; kwargs...) = _show_params(io, pip; kwargs...)

# Fallback implementation when PrettyTables is not loaded
function _show_params(io::IO, pip::ParameterIndexingProxy; num_rows = 20,
        show_all = false, scalarize = true, kwargs...)
    params = Any[]
    vals = Any[]
    for p in parameter_symbols(pip.wrapped)
        if symbolic_type(p) === ArraySymbolic() && scalarize
            val = getp(pip.wrapped, p)(pip.wrapped)
            for (_p, _v) in zip(collect(p), val)
                push!(params, _p)
                push!(vals, _v)
            end
        else
            push!(params, p)
            val = getp(pip.wrapped, p)(pip.wrapped)
            push!(vals, val)
        end
    end

    num_shown = if show_all
        length(params)
    else
        if num_rows > length(params)
            length(params)
        else
            num_rows
        end
    end

    # Fallback implementation without PrettyTables
    println(io, "Parameter Indexing Proxy")
    println(io, "=" ^ 50)
    println(io, "Parameter                | Value")
    println(io, "-" ^ 50)
    for i in 1:num_shown
        println(io, rpad(string(params[i]), 24) * " | " * string(vals[i]))
    end

    if num_shown < length(params)
        println(io,
            "$num_shown of $(length(params)) params shown. To show all the parameters, call `show_params(io, ps, show_all = true)`. Adjust the number of rows with the num_rows kwarg. Consult `show_params` docstring for more options.")
    end
end
