module ClimaParams

using TOML
import Dates: DateTime

export ParamDict

export float_type,
    get_parameter_values,
    write_log_file,
    log_parameter_information,
    create_toml_dict,
    merge_toml_files,
    get_tagged_parameter_values,
    get_tagged_parameter_names,
    fuzzy_match

const NAMESTYPE =
    Union{AbstractVector{S}, NTuple{N, S} where {N}} where {S <: AbstractString}

"""
    ParamDict{FT}

A concrete parameter dictionary that stores parameter data from TOML files.

This struct holds the effective set of parameters (defaults merged with any
overrides) and tracks which override parameters have been used.

# Fields
- `data::Dict`: The main dictionary holding the complete, merged set of parameter values and their metadata.
- `override_dict::Union{Nothing, Dict}`: A dictionary containing only the parameters from an override file, used for tracking purposes. Is `nothing` if no override file was provided.
"""
struct ParamDict{FT}
    "The main dictionary holding the complete, merged set of parameter values and their metadata."
    data::Dict
    "A dictionary containing only the parameters from an override file, used for tracking purposes. Is `nothing` if no override file was provided."
    override_dict::Union{Nothing, Dict}
end

"""
    float_type(pd::ParamDict)

Returns the float type `FT` with which the parameter dictionary `pd` was initialized.
"""
float_type(::ParamDict{FT}) where {FT} = FT

Base.iterate(pd::ParamDict, state) = Base.iterate(pd.data, state)
Base.iterate(pd::ParamDict) = Base.iterate(pd.data)


"""
    Base.getindex(pd::ParamDict, key)

Retrieves a parameter by its `key`, converting it to the type specified in the TOML file.

This allows for direct, dictionary-like access to the typed value of a parameter.

# Arguments
- `pd::ParamDict`: The parameter dictionary.
- `key`: The name of the parameter to retrieve.

# Returns
- The parameter's value, cast to the type defined in its metadata (e.g., `Float64`, `Int`, `Bool`).

# Examples

    toml_dict = CP.create_toml_dict(Float64)
    param_value = toml_dict["planet_radius"]  # Returns the value, e.g., 6.371e6

"""
function Base.getindex(pd::ParamDict, i)
    param_data = getindex(pd.data, i)
    param_value = param_data["value"]
    param_type = get(param_data, "type", nothing)
    isnothing(param_type) && error("No type found for parameter `$i`")

    if param_value isa AbstractVector
        val = map(x -> _get_typed_value(pd, x, i, param_type), param_value)
    else
        val = _get_typed_value(pd, param_value, i, param_type)
    end
    log_component!(pd, (i,), "getindex")
    return val
end

"""
    log_component!(pd::ParamDict, names::NAMESTYPE, component::AbstractString)

Logs that a set of parameters are used by a specific model `component`.

This function modifies the parameter dictionary in-place by adding or appending
the `component` string to a `"used_in"` entry for each parameter specified in `names`.
This is crucial for tracking which parameters are active in a simulation.

# Arguments
- `pd::{ParamDict}`: The parameter dictionary to be modified.
- `names`: A vector or tuple of strings with the names of parameters to log.
- `component::{AbstractString}`: The name of the model component using the parameters.
"""
function log_component!(
    pd::ParamDict,
    names::NAMESTYPE,
    component::AbstractString,
)
    component_key = "used_in"
    data = pd.data
    for name in names
        for (key, val) in data
            name ≠ key && continue
            data[key][component_key] = if component_key in keys(data[key])
                unique([data[key][component_key]..., component])
            else
                [component]
            end
        end
    end
end

"""
    _get_typed_value(pd::ParamDict, val, valname, valtype)

An internal helper function that converts a raw parameter `val` to the
correct Julia type based on the `valtype` string from the TOML file.

# Arguments
- `pd::ParamDict`: The parameter dictionary, used to get the float type.
- `val`: The raw value of the parameter.
- `valname::{AbstractString}`: The name of the parameter (for error messages).
- `valtype::{AbstractString}`: The type string, e.g., "float", "integer", "string", "bool", "datetime".

# Returns
- The value `val` converted to the appropriate type. Throws an error for an unknown `valtype`.
"""
function _get_typed_value(
    pd::ParamDict{FT},
    val,
    valname::AbstractString,
    valtype,
) where {FT}
    valid_types = Dict(
        "float" => FT,
        "integer" => Int,
        "string" => String,
        "bool" => Bool,
        "datetime" => DateTime,
    )
    if valtype in keys(valid_types)
        return valid_types[valtype](val)
    else
        error(
            """ For parameter with identifier: $valname, the attribute: type = $valtype, is not recognised, 
            please select from: $(keys(valid_types))
            """,
        )
    end
end

"""
    get_parameter_values(pd, names, [component])
    get_parameter_values(pd, name_map, [component])

Retrieves parameter values from the dictionary `pd`, returning them in a `NamedTuple`.
This function has two main methods:

1.  Retrieve parameters by a list of `names`.
2.  Retrieve and rename parameters using a `name_map`.

If a `component` string is provided, it also logs the parameters as being used by that component.

# Arguments
- `pd::ParamDict`: The parameter dictionary.
- `names::Union{String,Vector{String}}`: A single name or vector of names to retrieve.
- `name_map`: A `Dict` or other iterable of `Pair`s mapping the parameter name in the TOML file to the desired variable name in the code (e.g., `"long_name_in_toml" => "short_name_in_code"`).
- `component::Union{AbstractString, Nothing}`: An optional string to log which model component uses these parameters.

# Returns
- A `NamedTuple` where keys are the parameter names (or the renamed variable names) and values are the corresponding typed parameter values.

# Examples

    # Method 1: Retrieve by name
    params = get_parameter_values(toml_dict, ["gravitational_acceleration", "planet_radius"])
    # params.gravitational_acceleration = 9.81

    # Method 2: Retrieve and rename
    name_map = Dict("gravitational_acceleration" => "g", "planet_radius" => "R_p")
    params_renamed = get_parameter_values(toml_dict, name_map)
    # params_renamed.g = 9.81

"""
function get_parameter_values(
    pd::ParamDict,
    names::AbstractString,
    component = nothing,
)
    return get_parameter_values(pd, [names], component)
end

function get_parameter_values(
    pd::ParamDict,
    names::NAMESTYPE,
    component::Union{AbstractString, Nothing} = nothing,
)
    if !isnothing(component)
        log_component!(pd, names, component)
    end
    return NamedTuple(map(x -> Symbol(x) => pd[x], names))
end



function get_parameter_values(
    pd::ParamDict,
    name_map::Union{AbstractVector{Pair{S, S}}, NTuple{N, Pair}},
    component = nothing,
) where {S, N}
    return get_parameter_values(pd, Dict(name_map), component)
end

function get_parameter_values(
    pd::ParamDict,
    name_map::Vararg{Pair};
    component = nothing,
)
    return get_parameter_values(
        pd,
        Dict(Symbol(key) => Symbol(value) for (key, value) in name_map),
        component,
    )
end

function get_parameter_values(
    pd::ParamDict,
    name_map::Dict{S, S},
    component = nothing,
) where {S <: AbstractString}

    return get_parameter_values(
        pd,
        Dict(Symbol(key) => Symbol(value) for (key, value) in name_map),
        component,
    )
end

function get_parameter_values(
    pd::ParamDict,
    name_map::NamedTuple,
    component = nothing,
)
    return get_parameter_values(pd, Dict(pairs(name_map)), component)
end

function get_parameter_values(
    pd::ParamDict,
    name_map::Dict{Symbol, Symbol},
    component = nothing,
)
    params = get_parameter_values(pd, string.(keys(name_map)), component)
    return (;
        [
            short_name => getfield(params, long_name) for
            (long_name, short_name) in name_map
        ]...
    )
end

"""
    create_parameter_struct(param_struct_type, toml_dict, name_map, [nested_structs])

Constructs an instance of a parameter struct from a TOML dictionary.

This function retrieves all necessary parameter values using a `name_map` and
instantiates the `param_struct_type`, including any `nested_structs`.

This function makes several assumptions about the parameter struct:
- It has a constructor that accepts keyword arguments for its fields.
- Its first type parameter is the floating-point type (e.g., `MyParams{FT}`).
- All nested parameter structs required by the constructor are passed via `nested_structs`.

# Arguments
- `param_struct_type`: The type of the parameter struct to be created (e.g., `MyParams`).
- `toml_dict::ParamDict`: The TOML dictionary containing the parameter values.
- `name_map`: A `Dict` or other iterable of `Pair`s to map TOML names to struct field names.
- `nested_structs`: A `NamedTuple` of already-constructed nested parameter structs, if any.
"""
function create_parameter_struct(
    param_struct_type,
    toml_dict,
    name_map,
    nested_structs = (;),
)
    params = get_parameter_values(toml_dict, name_map)
    FT = float_type(toml_dict)
    return param_struct_type{FT, typeof.(values(nested_structs))...}(;
        params...,
        nested_structs...,
    )
end

"""
    merge_toml_files(filepaths; override::Bool=false)

Parses and merges multiple TOML files into a single dictionary.

# Arguments
- `filepaths`: An iterable of strings, where each string is a path to a TOML file.
- `override::Bool`: If `false` (the default), an error is thrown for duplicate parameter entries across files. If `true`, a warning is issued and later files in the `filepaths` list will overwrite earlier entries.

# Returns
- `Dict{String, Any}`: A dictionary containing the merged data from all TOML files.
"""
function merge_toml_files(filepaths; override = false)
    merged_dict = Dict{String, Any}()
    for filepath in filepaths
        toml_data = TOML.parsefile(filepath)
        for (table_name, table_data) in toml_data
            if haskey(merged_dict, table_name)
                override || error("Duplicate TOML entry: $table_name")
                @warn """
'$table_name' is being overwritten by '$filepath'
Current entry: $(merged_dict[table_name]["type"])($(merged_dict[table_name]["value"]))
New entry: $(table_data["type"])($(table_data["value"]))"""
            end
        end
        merge!(merged_dict, toml_data)
    end
    return merged_dict
end

"""
    check_override_parameter_usage(pd::ParamDict, strict::Bool)

Verifies that all parameters supplied in an override file were actually used
during the simulation by checking for the `"used_in"` log entry.

# Arguments
- `pd::{ParamDict}`: The parameter dictionary to check.
- `strict::Bool`: If `true`, throws an error if any override parameter is unused. If `false`, only a warning is issued.
"""
check_override_parameter_usage(pd::ParamDict, strict::Bool) =
    check_override_parameter_usage(pd, strict, pd.override_dict)

check_override_parameter_usage(pd::ParamDict, strict::Bool, ::Nothing) = nothing

function check_override_parameter_usage(
    pd::ParamDict,
    strict::Bool,
    override_dict,
)
    unused_override = Dict()
    for (key, _) in override_dict
        logged_val = pd.data[key]
        unused_override[key] = !("used_in" in keys(logged_val))
    end
    if any(values(unused_override))
        unused_override_keys = collect(keys(unused_override))
        filter!(key -> unused_override[key], unused_override_keys)
        @warn(
            string(
                "Keys are present in parameter file but not used ",
                "in the simulation. \n Typically this is due to ",
                "a mismatch in parameter name in toml and in source. ",
                "Offending keys: $(unused_override_keys)",
            )
        )
        if strict
            @error(
                "At least one override parameter set and not used in simulation"
            )
            error(
                "Halting simulation due to unused parameters." *
                "\n Typically this is due to a typo in the parameter name." *
                "\n change `strict` flag to `false` to prevent this causing an exception",
            )
        end
    end
    return nothing
end

"""
    write_log_file(pd::ParamDict, filepath::AbstractString)

Saves all *used* parameters to a TOML file at the specified `filepath`.

This function filters the dictionary to include only parameters that have been
logged with [`log_component!`](@ref), creating a file that can be used to
reproduce an experiment with the exact same parameter set.

# Arguments
- `pd::ParamDict`: The parameter dictionary containing usage logs.
- `filepath::{AbstractString}`: The path where the log file will be saved.
"""
function write_log_file(pd::ParamDict, filepath::AbstractString)
    used_parameters = Dict()
    for (key, val) in pd.data
        if "used_in" in keys(val)
            used_parameters[key] = val
        end
    end
    open(filepath, "w") do io
        TOML.print(io, used_parameters)
    end
end


"""
    log_parameter_information(pd::ParamDict, filepath; strict::Bool = false)

A convenience function that performs end-of-run parameter handling.

It calls [`write_log_file`](@ref) to save used parameters and then
[`check_override_parameter_usage`](@ref) to validate that all override
parameters were used.

# Arguments
- `pd::ParamDict`: The parameter dictionary.
- `filepath::{AbstractString}`: The path for the output log file.
- `strict::Bool`: If `true`, errors if override parameters are unused.
"""
function log_parameter_information(
    pd::ParamDict,
    filepath::AbstractString;
    strict::Bool = false,
)
    #[1.] write the parameters to log file
    write_log_file(pd, filepath)
    #[2.] send warnings or errors if parameters were not used
    check_override_parameter_usage(pd, strict)
end

"""
    merge_override_default_values(override_toml_dict, default_toml_dict)

An internal helper that merges two `ParamDict` objects, with values from the
`override_toml_dict` taking precedence over the `default_toml_dict`.
"""
function merge_override_default_values(
    override_toml_dict::ParamDict{FT},
    default_toml_dict::ParamDict{FT},
) where {FT}
    data = default_toml_dict.data
    override_dict = override_toml_dict.override_dict
    for (key, val) in override_toml_dict.data
        if !(key in keys(data))
            data[key] = val
        else
            for (kkey, vval) in val # as val is a Dict too
                data[key][kkey] = vval
            end
        end
    end
    return ParamDict{FT}(data, override_dict)
end

"""
    create_toml_dict(
        FT;
        override_file::Union{String, Dict, Nothing}=nothing,
        default_file::Union{String, Dict}="parameters.toml",
    )

Creates a `ParamDict{FT}` by reading and merging default and override parameter sources.

This is the main entry point for constructing a parameter dictionary. It reads a
`default_file` and optionally an `override_file`, with parameters from the
override file taking precedence. The sources can be file paths or already-parsed
Julia `Dict`s.

# Arguments
- `FT::{Type{<:AbstractFloat}}`: The floating-point type to be used for all "float" parameters.

# Keywords
- `override_file`: Path to a TOML file or a `Dict` containing override parameters.
- `default_file`: Path to the default TOML file or a `Dict` containing default parameters. Defaults to the `parameters.toml` file in the package directory.

# Returns
- A `ParamDict{FT}` containing the merged and typed parameters.
"""
function create_toml_dict(
    ::Type{FT};
    override_file::Union{Nothing, String, Dict} = nothing,
    default_file::Union{String, Dict} = joinpath(@__DIR__, "parameters.toml"),
) where {FT <: AbstractFloat}

    default_dict =
        default_file isa String ? TOML.parsefile(default_file) : default_file
    default_toml_dict = ParamDict{FT}(default_dict, nothing)
    isnothing(override_file) && return default_toml_dict

    override_dict =
        override_file isa String ? TOML.parsefile(override_file) : override_file
    override_toml_dict = ParamDict{FT}(override_dict, override_dict)

    return merge_override_default_values(override_toml_dict, default_toml_dict)
end

Base.print(td::ParamDict, io = stdout) = TOML.print(io, td.data)

function Base.show(io::IO, d::ClimaParams.ParamDict{FT}) where {FT}
    n = length(d.data)
    print(io, "ParamDict{$FT} with $n parameters")
end

function Base.:(==)(pd1::ParamDict{FT1}, pd2::ParamDict{FT2}) where {FT1, FT2}
    return FT1 == FT2 &&
           pd1.data == pd2.data &&
           pd1.override_dict == pd2.override_dict
end

"""
    get_tagged_parameter_names(pd::ParamDict, tag)

Retrieves the names of all parameters associated with a given `tag` or list of `tags`.

Tag matching is case-insensitive and ignores punctuation and whitespace.

# Arguments
- `pd::ParamDict`: The parameter dictionary.
- `tag::Union{AbstractString, Vector{<:AbstractString}}`: The tag or vector of tags to search for.

# Returns
- `Vector{String}`: A list of parameter names that have the specified tag(s).
"""
function get_tagged_parameter_names(pd::ParamDict, tag::AbstractString)
    data = pd.data
    ret_values = String[]
    for (key, val) in data
        if any(fuzzy_match.(tag, get(val, "tag", [])))
            push!(ret_values, key)
        end
    end
    return ret_values
end

get_tagged_parameter_names(
    pd::ParamDict,
    tags::Vector{S},
) where {S <: AbstractString} =
    vcat(map(x -> get_tagged_parameter_names(pd, x), tags)...)

"""
    fuzzy_match(s1::AbstractString, s2::AbstractString)

Compares two strings for equality, ignoring case and select punctuation.

The characters `[' ', '_', '*', '.', ',', '-', '(', ')']` are stripped from both strings before comparison.
"""
function fuzzy_match(s1::AbstractString, s2::AbstractString)
    strip_chars(x) = replace(x, [' ', '_', '*', '.', ',', '-', '(', ')'] => "")
    return lowercase(strip_chars(s1)) == lowercase(strip_chars(s2))
end

"""
    get_tagged_parameter_values(pd::ParamDict, tag)

Retrieves the values of all parameters associated with a given `tag` or list of `tags`.

# Arguments
- `pd::ParamDict`: The parameter dictionary.
- `tag::Union{AbstractString, Vector{<:AbstractString}}`: The tag or vector of tags to search for.

# Returns
- A `NamedTuple` of the tagged parameters, where keys are parameter names and values are their typed values.
"""
get_tagged_parameter_values(pd::ParamDict, tag::AbstractString) =
    get_parameter_values(pd, get_tagged_parameter_names(pd, tag))

get_tagged_parameter_values(
    pd::ParamDict,
    tags::Vector{S},
) where {S <: AbstractString} =
    merge(map(x -> get_tagged_parameter_values(pd, x), tags)...)

end # module
