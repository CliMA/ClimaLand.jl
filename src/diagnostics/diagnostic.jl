# ClimaLand diagnostics contains # a dictionary `ALL_DIAGNOSTICS` with all the
# diagnostics we know how to compute, keyed over their short name. If you want
# to add more diagnostics, look at the included files. You can add your own file
# if you want to define several new diagnostics that are conceptually related.
# The dictionary `ALL_DIAGNOSTICS` should be considered an implementation
# detail, use the getters/setters.

const ALL_DIAGNOSTICS = Dict{String, DiagnosticVariable}()

"""

    add_diagnostic_variable!(; short_name,
                               long_name,
                               standard_name,
                               units,
                               description,
                               compute!)


Add a new variable to the `ALL_DIAGNOSTICS` dictionary (this function mutates the state of
`ClimaLand.ALL_DIAGNOSTICS`).

If possible, please follow the naming scheme outline in
https://airtable.com/appYNLuWqAgzLbhSq/shrKcLEdssxb8Yvcp/tblL7dJkC3vl5zQLb

Keyword arguments
=================

- `short_name`: Name used to identify the variable in the output files and in the file
                names. Short but descriptive. `ClimaLand` diagnostics are identified by the
                short name. We follow the Coupled Model Intercomparison Project conventions.

- `long_name`: Name used to identify the variable in the output files.

- `standard_name`: Standard name, as in
                   http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html

- `units`: Physical units of the variable.

- `comments`: More verbose explanation of what the variable is, or comments related to how
              it is defined or computed.

- `compute!`: Function that compute the diagnostic variable from the state.
                             It has to take two arguments: the `integrator`, and a
                             pre-allocated area of memory where to write the result of the
                             computation. If no pre-allocated area is available, a new
                             one will be allocated. To avoid extra allocations, this
                             function should perform the calculation in-place (i.e., using
                             `.=`).

"""
function add_diagnostic_variable!(;
    short_name,
    long_name,
    standard_name = "",
    units,
    comments = "",
    compute!,
)
    ALL_DIAGNOSTICS[short_name] = DiagnosticVariable(;
        short_name,
        long_name,
        standard_name,
        units,
        comments,
        compute!,
    )
end

"""

    conditional_add_diagnostic_variable!(
                               possible_diags;
                               short_name,
                               long_name,
                               standard_name,
                               units,
                               description,
                               compute!)


Calls `add_diagnostic_variable!` if `short_name` is in `possible_diags`, otherwise does nothing.
If `possible_diags` is `Val{:all}()`, then the diagnostic variable is added regardless of its short name.

This is useful to avoid adding diagnostics that are not relevant for a particular model,
which can reduce compilation and inference time.
"""
function conditional_add_diagnostic_variable!(
    possible_diags::Union{AbstractArray, Tuple};
    short_name,
    kwargs...,
)
    if short_name in possible_diags
        add_diagnostic_variable!(; short_name, kwargs...)
    end
    return
end

conditional_add_diagnostic_variable!(::Val{:all}; kwargs...) =
    add_diagnostic_variable!(; kwargs...)
conditional_add_diagnostic_variable!(possible_diags; kwargs...) = error(
    "Invalid argument for possible_diags. It should be either `Val{:all}()` or a collection of short names of diagnostics to add.",
)


"""

    get_diagnostic_variable!(short_name)

Return a `DiagnosticVariable` from its `short_name`, if it exists.
"""
function get_diagnostic_variable(short_name)
    haskey(ALL_DIAGNOSTICS, short_name) ||
        error("diagnostic $short_name does not exist")

    return ALL_DIAGNOSTICS[short_name]
end

# General helper functions for undefined diagnostics for a particular model
error_diagnostic_variable(variable, land_model::T) where {T} =
    error("Cannot compute $variable with model = $(nameof(T))")

# with_error is a helper macro that generates the error message
# when the user tries calling something that is incompatible with the model.
# It should be called when defining compute functions
macro with_error(compute_function_expr)
    compute_function_expr.head == :function ||
        error("Cannot parse this function, head is not a :function")

    # Two firsts:
    # 1st: extract the function signature
    # 2nd: extract the name
    function_name = first(first(compute_function_expr.args).args)
    function_name isa Symbol || error("Cannot parse this function!")

    # qualified_name ensures that this macro can be used outside of this module while
    # still defining compute functions in this module
    qualified_name = GlobalRef(Diagnostics, function_name)
    # Assuming the convention that functions are called "compute_variable!",
    # otherwise the error might look a little less informative
    variable_name = replace(string(function_name), "compute_" => "", "!" => "")
    return esc(
        quote
            # Paste back the definition of the function
            $compute_function_expr
            # And add the error method, unless it was added by another macro call
            if !(hasmethod($qualified_name, (Any, Any, Any, Any, Any)))
                function $qualified_name(_, _, _, _, land_model)
                    error_diagnostic_variable($variable_name, land_model)
                end
            end
        end,
    )
end

"""
    nlayers(field::Fields.Field)

Returns the number of layers in the vertical for the `field`;
this differs from the ClimaCore.Spaces `nlevels` function in
that it considers non-extruded spaces (spectral element or point
spaces) to have a single layer.
"""
nlayers(field::Fields.Field) = nlayers(axes(field))

"""
    nlayers(space::Spaces.AbstractSpace)

Returns the default number of layers in the vertical for a ClimaCore
space. For extruded spaces, this is equal to the number of ClimaCore levels.
"""
nlayers(space::Spaces.AbstractSpace) = Spaces.nlevels(space)

"""
    nlayers(::Spaces.PointSpace)
    nlayers(::Spaces.MultiPointSpace)

Returns the number of layers in the vertical; this is one
for the PointSpace and MultiPointSpace.
"""
nlayers(::Spaces.PointSpace) = 1
nlayers(::Spaces.MultiPointSpace) = 1

"""
    depth_interpolation_weights(z, depth)

Return `(i_lo, i_hi, w_hi)` such that the value at `depth` [m] below the
surface of a quantity known at the layer center heights `z` [m] (increasing
upward, negative below the surface) is `(1 - w_hi) v[i_lo] + w_hi v[i_hi]`, by
linear interpolation between the bracketing centers. Depths above the top
center or below the bottom center take the value of that center.
"""
function depth_interpolation_weights(z::AbstractVector, depth)
    z_target = -depth
    z_target >= z[end] && return (lastindex(z), lastindex(z), zero(z_target))
    z_target <= z[1] && return (firstindex(z), firstindex(z), zero(z_target))
    i_hi = findfirst(>(z_target), z)
    i_lo = i_hi - 1
    w_hi = (z_target - z[i_lo]) / (z[i_hi] - z[i_lo])
    return (i_lo, i_hi, w_hi)
end

"""
    layer_center_heights(field::Fields.Field)

Return the heights [m] of the layer centers of a column `field`, from the
bottom up, as a host `Array` (the field may live on a GPU).
"""
layer_center_heights(field::Fields.Field) =
    Array(vec(parent(Fields.coordinate_field(axes(field)).z)[:, 1]))

"""
    at_depth(field::Fields.Field, depth)

Return the value of the column `field` at `depth` [m] below the surface,
interpolated linearly between layer centers (see `depth_interpolation_weights`).
"""
function at_depth(field::Fields.Field, depth)
    (i_lo, i_hi, w_hi) =
        depth_interpolation_weights(layer_center_heights(field), depth)
    values = Array(parent(field))
    return (1 - w_hi) * values[i_lo] + w_hi * values[i_hi]
end

"""
    diagnostic_as_vectors(writer::ClimaDiagnostics.DictWriter, diagnostic;
                          layer = nothing, depth = nothing)

Extract `diagnostic` from given `writer` as tuple of vectors (time and value).
By default, if `layer` and `depth` are nothing, it gets the surface value;
otherwise it returns the layer requested, or the value interpolated linearly
between layer centers to `depth` [m] below the surface (see
[`depth_interpolation_weights`](@ref)). If both are given, `depth` is used and
`layer` is ignored. Variables without vertical resolution ignore `depth`.

For variables resolved in depth, the bottom layer is indicated by `1` and the
top layer by the number of layers.

`diagnostic` is typically a string with the short name of the diagnostic.
"""
function diagnostic_as_vectors(
    writer::DictWriter,
    diagnostic;
    layer = nothing,
    depth = nothing,
)

    # writer[diagnostic] is a dictionary with keys the times and with values Fields. We need
    # to be a little careful because dictionaries are not ordered, so we have to sort them
    # by time.
    times = collect(keys(writer[diagnostic]))
    sort_indices = sortperm(times)
    values_all = parent.(values(writer[diagnostic]))[sort_indices]
    field = first(values(writer[diagnostic]))
    # Variables without vertical resolution ignore `depth`
    if isnothing(depth) || nlayers(field) == 1
        layer_id = layer isa Nothing ? nlayers(field) : layer
        vector_layer =
            vcat([values_all[i][layer_id, :] for i in eachindex(values_all)]...)
    else
        (i_lo, i_hi, w_hi) =
            depth_interpolation_weights(layer_center_heights(field), depth)
        vector_layer = vcat(
            [
                (1 - w_hi) .* values_all[i][i_lo, :] .+
                w_hi .* values_all[i][i_hi, :] for
                i in eachindex(values_all)
            ]...,
        )
    end

    return times, vector_layer
end

"""
    close_output_writers(diagnostics)

Close the output writers in the `diagnostics`, an iterable of
`ClimaDiagnostics.ScheduledDiagnostic` or `nothing`.

This function should be called at the end of every simulation.
"""
function close_output_writers(diagnostics)
    isnothing(diagnostics) && return nothing
    for diagnostic in diagnostics
        close(diagnostic.output_writer)
    end
    return nothing
end


function default_diagnostic_num_points(
    domain::Union{SphericalShell, SphericalSurface},
)
    num_horizontal_elements, num_vertical_elements = domain.nelements

    # 4 panels cover the cubed sphere from -180 to 180 in long
    num_long = 4num_horizontal_elements
    # 2 panels cover the cubed sphere from 0 to 180 in lat
    num_lat = 2num_horizontal_elements

    num_z = num_vertical_elements

    return num_long, num_lat, num_z
end

function default_diagnostic_num_points(domain::Union{HybridBox, Plane})
    num_lat, num_long, num_z = domain.nelements
    return num_long, num_lat, num_z
end

function default_diagnostic_num_points(domain::Union{Point, Column})
    return domain.nelements
end

# Do you want to define more diagnostics? Add them here
include("land_compute_methods.jl")

# define_diagnostics.jl contains the list of all the diagnostics
include("define_diagnostics.jl")

# Default diagnostics and higher level interfaces
include("default_diagnostics.jl")
