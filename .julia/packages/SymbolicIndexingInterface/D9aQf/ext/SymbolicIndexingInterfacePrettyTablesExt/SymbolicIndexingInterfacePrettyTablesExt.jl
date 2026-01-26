module SymbolicIndexingInterfacePrettyTablesExt

using SymbolicIndexingInterface
using SymbolicIndexingInterface: ParameterIndexingProxy, parameter_symbols, symbolic_type,
                                 ArraySymbolic, getp
using PrettyTables

# Override the fallback implementation with the PrettyTables version
function SymbolicIndexingInterface.show_params(
        io::IO, pip::ParameterIndexingProxy; num_rows = 20,
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

    pretty_table(io, [params[1:num_shown] vals[1:num_shown]];
        column_labels = ["Parameter", "Value"],
        kwargs...)

    if num_shown < length(params)
        println(io,
            "$num_shown of $(length(params)) params shown. To show all the parameters, call `show_params(io, ps, show_all = true)`. Adjust the number of rows with the num_rows kwarg. Consult `show_params` docstring for more options.")
    end
end

end
