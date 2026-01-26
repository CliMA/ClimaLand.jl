function get_param_lim()
    config = CUDA.compiler_config(CUDA.device())
    (; ptx, cap) = config.params
    return cap >= v"7.0" && ptx >= v"8.1" ? 32764 : 4096
end
param_usage(arg) = sizeof(typeof(CUDA.cudaconvert(arg)))
param_usage_args(args) =
    sum(x -> param_usage(x), args) + param_usage(CUDA.KernelState)

function fused_multibroadcast_args(fmb::MBF.FusedMultiBroadcast)
    dest = first(fmb.pairs).first
    CI = CartesianIndices(axes(dest))
    return (fmb, CI)
end

"""
    Options(
        [types...];
        match_only::Bool = false
        print_types::Bool = false
        recursion_types = (UnionAll,DataType)
        recursion_depth = 1000
    )

Printing options for `@rprint_parameter_memory`:

 - `match_only`: only print properties that match the given types
 - `print_types`: print types (e.g., `prop::typeof(prop)`)
 - `recursion_types`: skip recursing through recursion types (e.g., `UnionAll` and `DataType`)
                      to avoid infinite recursion
 - `recursion_depth`: limit recursion depth (to avoid infinite recursion)
"""
struct Options{T}
    types::T
    match_only::Bool
    print_types::Bool
    recursion_types::Tuple
    recursion_depth::Int
    size_threshhold::Int
    max_type_depth::Int
    function Options(
        types...;
        match_only = false,
        print_types = true,
        recursion_types = (UnionAll, DataType),
        recursion_depth = 1000,
        size_threshhold = 10,
        max_type_depth = 1,
    )
        if (types isa AbstractArray || types isa Tuple) && length(types) > 0
            types = types[1]
        else
            types = (Union{},)
        end
        return new{typeof(types)}(
            types,
            match_only,
            print_types,
            recursion_types,
            recursion_depth,
            size_threshhold,
            max_type_depth,
        )
    end
end
Options(type::Type; kwargs...) = Options((type,); kwargs...)

Options() = Options(();)

function type_string(io, obj; maxdepth)
    sz = get(io, :displaysize, displaysize(io))::Tuple{Int, Int}
    S = max(sz[2], 120)
    slim = Base.type_depth_limit(string(typeof(obj)), S; maxdepth)
    return slim
end

function _rprint_parameter_memory(io, obj, pc; o::Options, name, counter = 0)
    counter > o.recursion_depth && return
    for pn in propertynames(obj)
        prop = getproperty(obj, pn)
        pc_full = (pc..., ".", pn)
        pc_string = name * string(join(pc_full))
        if any(map(type -> prop isa type, o.types))
            suffix =
                o.print_types ?
                "::$(type_string(io, prop; maxdepth=o.max_type_depth))" : ""
            s = sizeof(typeof(CUDA.cudaconvert(prop)))
            if s > o.size_threshhold
                println(io, "size: $s, $pc_string$suffix")
            end
            if !any(map(x -> prop isa x, o.recursion_types))
                _rprint_parameter_memory(
                    io,
                    prop,
                    pc_full;
                    o,
                    name,
                    counter = counter + 1,
                )
                counter > o.recursion_depth && return
            end
        else
            if !o.match_only
                suffix =
                    o.print_types ?
                    "::$(type_string(io, prop; maxdepth=o.max_type_depth))" : ""
                s = sizeof(typeof(CUDA.cudaconvert(prop)))
                if s > o.size_threshhold
                    println(io, "size: $s, $pc_string$suffix")
                end
            end
            if !any(map(x -> prop isa x, o.recursion_types))
                _rprint_parameter_memory(
                    io,
                    prop,
                    pc_full;
                    o,
                    name,
                    counter = counter + 1,
                )
            end
            counter > o.recursion_depth && return
        end
    end
end

print_name(io, name, o) = o.match_only || println(io, name)

function rprint_parameter_memory(io, obj, name, o::Options = Options())
    print_name(io, name, o)
    _rprint_parameter_memory(
        io,
        obj,
        (); # pc
        o,
        name,
    )
    println(io, "")
end

"""
    @rprint_parameter_memory obj options

Recursively print out propertynames and
parameter memory of `obj` given options
`options`. See [`Options`](@ref) for more
information on available options.
"""
macro rprint_parameter_memory(obj, o)
    return :(rprint_parameter_memory(
        stdout,
        $(esc(obj)),
        $(string(obj)),
        $(esc(o)),
    ))
end

"""
    @rprint_parameter_memory obj options

Recursively print out propertynames and
parameter memory of `obj` given options
`options`. See [`Options`](@ref) for more
information on available options.
"""
macro rprint_parameter_memory(obj)
    return :(rprint_parameter_memory(stdout, $(esc(obj)), $(string(obj))))
end
