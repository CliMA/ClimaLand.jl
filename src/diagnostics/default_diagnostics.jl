using Dates
export default_diagnostics

# This file is included by Diagnostics.jl and defines all the defaults for
# various models (e.g., Bucket, SoilCanopyModel). A model here is either a
# standalone (e.g., Bucket) or integrated (e.g., SoilCanopy) model.
#
# If you are developing new models, add your defaults here. If you want to add
# more high level interfaces, add them here. Feel free to include extra files.

"""
    function common_diagnostics(reduction_period,
                                reduction_type,
                                output_writer,
                                start_date,
                                short_names...;
                                pre_output_hook! = nothing,
                               )

For each variable specified in `short_names`, create a `ClimaDiagnostics.ScheduledDiagnostic`.
Diagnostics are computed at every step, and output at the specified period with the provided reduction.

The type of `output_writer` determines where the output will be saved to.
This must be a `ClimaDiagnostics.AbstractWriter` object: concrete options include
`NetCDFWriter`, `HDF5Writer`, and `DictWriter`.

The provided `pre_output_hook!` is a function that is called after applying the reduction,
but before outputting the diagnostic. This is typically not needed (`nothing`), but to
compute an average we use `ClimaDiagnostics.average_pre_output_hook!`.
"""
function common_diagnostics(
    reduction_period,
    reduction_type,
    output_writer,
    start_date,
    short_names...;
    dt = nothing,
    pre_output_hook! = nothing,
)
    # Convert the provided period and reduction to the appropriate types
    period = get_period(reduction_period, dt)
    reduction = get_reduction(reduction_type)
    reduction_type == Val(:average) &&
        (pre_output_hook! = average_pre_output_hook!)
    # Make a list of ScheduledDiagnostic objects
    return vcat(
        map(short_names) do short_name
            output_schedule_func =
                EveryCalendarDtSchedule(period; start_date)
            return ScheduledDiagnostic(
                variable = get_diagnostic_variable(short_name),
                compute_schedule_func = EveryStepSchedule(),
                output_schedule_func = output_schedule_func,
                reduction_time_func = reduction,
                output_writer = output_writer,
                pre_output_hook! = pre_output_hook!,
            )
        end...,
    )
end

"""
    get_period(val, dt)

Helper function to convert from a user-provided Symbol
to the corresponding Dates.Period.

Currently, the following periods are supported:
- :every_dt (requires `dt` to be provided separately)
- :halfhourly
- :hourly
- :daily
- :tendaily
- :monthly
"""

get_period(::Val{:every_dt}, dt) = Second(dt)
get_period(::Val{:halfhourly}, dt) = Minute(30)
get_period(::Val{:hourly}, dt) = Hour(1)
get_period(::Val{:daily}, dt) = Day(1)
get_period(::Val{:tendaily}, dt) = Day(10)
get_period(::Val{:monthly}, dt) = Month(1)
get_period(val, dt) = @error("Diagnostic reduction period $val not supported.")

"""
    get_reduction(val)

Helper function to convert from a user-provided Symbol
to the corresponding diagnostic reduction Function.

Currently, the following reduction types are supported:
- :instantaneous (should be used with `reduction_period = :every_dt`)
- :average (requires `pre_output_hook! = average_pre_output_hook!` when creating the `ScheduledDiagnostic`)
- :max
- :min

New methods of this function can be added to compute different types of reductions.
For any new methods, the returned function should be one that can act on a list of Float inputs.
"""
get_reduction(::Val{:instantaneous}) = nothing
get_reduction(::Val{:average}) = (+)
get_reduction(::Val{:max}) = max
get_reduction(::Val{:min}) = min
get_reduction(val) = @error("Diagnostic reduction type $val not supported.")

"""
    default_output_writer(domain::Union{SphericalShell, SphericalSurface}, start_date, outdir)

Creates a NetCDF diagnostics output writer using the start_date and outdir provided;
the default output coordinates are spaced evenly at locations which depend on the resolution
of the `domain`.
"""
function default_output_writer(
    domain::Union{SphericalShell, SphericalSurface},
    start_date,
    outdir,
)
    #If num_long = 360, num_lat = 180, the output points will align with the ERA5 1 degree grid
    num_long, num_lat, _ =
        ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
    Δ_long = 360.0 / num_long
    Δ_lat = 180.0 / num_lat
    longs = collect(range(-180.0; length = num_long, step = Δ_long))
    lats = collect(range(-90.0; length = num_lat, step = Δ_lat))
    space =
        haskey(domain.space, :subsurface) ? domain.space.subsurface :
        domain.space.surface
    output_writer =
        NetCDFWriter(space, outdir; start_date, horizontal_pts = (longs, lats))
    return output_writer
end

"""
    default_output_writer(domain::Union{Plane,HybridBox}, start_date, outdir)

Creates a NetCDF diagnostics output writer using the start_date and outdir provided;
the default output coordinates are spaced evenly at locations which depend on the resolution
of the `domain`.
"""
function default_output_writer(
    domain::Union{Plane, HybridBox},
    start_date,
    outdir,
)
    space =
        haskey(domain.space, :subsurface) ? domain.space.subsurface :
        domain.space.surface
    if domain.longlat isa Nothing
        num_x, num_y, num_z =
            ClimaLand.Diagnostics.default_diagnostic_num_points(domain)
        output_writer = NetCDFWriter(
            space,
            outdir;
            start_date,
            num_points = (num_x, num_y, num_z),
        )
    else
        coords = ClimaLand.Domains.coordinates(domain).surface
        longs = unique(Array(parent(coords.long))[:])
        lats = unique(Array(parent(coords.lat))[:])
        output_writer = NetCDFWriter(
            space,
            outdir;
            start_date,
            horizontal_pts = (longs, lats),
        )
    end

    return output_writer
end

"""
    default_output_writer(domain::Union{Column, Point}}, start_date, outdir)

Creates an in memory diagnostics output writer for Column and Point domains.
"""
function default_output_writer(domain::Union{Column, Point}, start_date, outdir)
    output_writer = DictWriter()
    return output_writer
end

default_diagnostics(
    model::ClimaLand.AbstractModel,
    start_date::ITime{<:Any, <:Any, <:DateTime},
    outdir;
    kwargs...,
) = default_diagnostics(model, date(start_date), outdir; kwargs...)

# The default diagnostics currently require a start date because they use Dates.Period.
function default_diagnostics(
    model::ClimaLand.AbstractModel,
    start_date::Union{Nothing, ITime{<:Any, <:Any, Nothing}},
    outdir;
    kwargs...,
)
    @warn "Default diagnostics not available when running without a start date."
    return []
end

"""
    default_diagnostics(model::Union{
                            CanopyModel{FT},
                            SoilCanopyModel{FT},
                            BucketModel{FT},
                        },
                        start_date::DateTime,
                        outdir;
                        output_writer = default_output_writer(get_domain(model), start_date, outdir),
                        output_vars = :short,
                        reduction_period = :monthly,
                        reduction_type = :average,
                        dt = nothing)

Construct a list of `ScheduledDiagnostics` that outputs the given variables at the specified average period.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.
Please see the method `get_possible_diagnostics` for the list of available diagnostics for each model.

`reduction_period` specifies the frequency at which to average the diagnostics, and
`reduction_type` specifies the type of reduction to apply.
Please see the docstring of `get_period` for the list of available periods,
and the docstring of `get_reduction` for the list of available reduction types.

This method can be extended for any model that extends `get_possible_diagnostics` and `get_short_diagnostics`.
Note that `EnergyHydrology` and `LandModel` have a specialized method that handles conservation diagnostics.
"""
function default_diagnostics(
    model::Union{CanopyModel{FT}, SoilCanopyModel{FT}, BucketModel{FT}},
    start_date::DateTime,
    outdir;
    output_writer = default_output_writer(
        get_domain(model),
        start_date,
        outdir,
    ),
    output_vars = :short,
    reduction_period = :monthly,
    reduction_type = :average,
    dt = nothing,
) where {FT}
    define_diagnostics!(model)

    possible_diags = get_possible_diagnostics(model)
    if output_vars == :long
        diagnostics = possible_diags
    elseif output_vars == :short
        diagnostics = get_short_diagnostics(model)
    else
        @assert typeof(output_vars) <: Vector{String}
        @assert all([var in possible_diags for var in output_vars])
        diagnostics = output_vars
    end

    default_outputs = common_diagnostics(
        Val(reduction_period),
        Val(reduction_type),
        output_writer,
        start_date,
        diagnostics...;
        dt,
        pre_output_hook! = nothing,
    )

    return [default_outputs...]
end

"""
    default_diagnostics(
        land_model::Union{EnergyHydrology{FT},LandModel{FT}},
        start_date::DateTime,
        outdir;
        output_writer = default_output_writer(get_domain(model), start_date, outdir),
        output_vars = :short,
        reduction_period = :monthly,
        reduction_type = :average,
        conservation = false,
        conservation_period = Day(10),
        dt = nothing,
    ) where {FT}

Define a method specific to the EnergyHydrology and LandModel models so that we can
handle conservation diagnostics specially.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.

`reduction_period` specifies the frequency at which to average the diagnostics, and
`reduction_type` specifies the type of reduction to apply.
Please see the docstring of `get_period` for the list of available periods,
and the docstring of `get_reduction` for the list of available reduction types.

Conservation diagnostics should not be provided as part of the `output_vars` argument,
but rather included by providing `conservation = true`.
Please see the method `get_possible_diagnostics` for the list of available diagnostics.
"""
function default_diagnostics(
    land_model::Union{EnergyHydrology{FT}, LandModel{FT}},
    start_date::DateTime,
    outdir;
    output_writer = default_output_writer(
        get_domain(land_model),
        start_date,
        outdir,
    ),
    output_vars = :short,
    reduction_period = :monthly,
    reduction_type = :average,
    conservation = false,
    conservation_period = Day(10),
    dt = nothing,
) where {FT}

    define_diagnostics!(land_model)

    possible_diags = get_possible_diagnostics(land_model)
    if output_vars == :long
        diagnostics = possible_diags
    elseif output_vars == :short
        diagnostics = get_short_diagnostics(land_model)
    else
        @assert typeof(output_vars) <: Vector{String}
        @assert all([var in possible_diags for var in output_vars])
        diagnostics = output_vars
    end

    default_outputs = common_diagnostics(
        Val(reduction_period),
        Val(reduction_type),
        output_writer,
        start_date,
        diagnostics...;
        dt,
        pre_output_hook! = nothing,
    )

    if conservation
        additional_diags = ["epa", "epac", "wvpa", "wvpac"]
        additional_outputs = vcat(
            map(additional_diags) do short_name
                output_schedule_func =
                    conservation_period isa Period ?
                    EveryCalendarDtSchedule(conservation_period; start_date) : EveryDtSchedule(conservation_period)
                return ScheduledDiagnostic(
                    variable = get_diagnostic_variable(short_name),
                    compute_schedule_func = EveryStepSchedule(),
                    output_schedule_func = output_schedule_func,
                    output_writer = output_writer,
                )
            end...,
        )
    else
        additional_outputs = []
    end

    return [default_outputs..., additional_outputs...]
end

function default_diagnostics(
    model::AbstractModel,
    start_date = nothing,
    outdir = nothing;
    output_writer = nothing,
    output_vars = nothing,
    reduction_period = nothing,
)
    @warn(
        "No default diagnostics defined for model type $(nameof(typeof(model))); consider extending `default_diagnostics` for this model type."
    )
    return []
end

"""
    add_diagnostics!(diagnostics, model, subcomponent)

A function to add diagnostics for a specific subcomponent to the list of diagnostics
for the provided model. This should be extended for any model with diagnostics
that depend on the types of subcomponents, and then called from the corresponding
`get_possible_diagnostics` method.

The fallback method does nothing.
"""
add_diagnostics!(_, _::AbstractModel, _) = nothing

function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.AtmosDrivenFluxBC,
)
    append!(diagnostics, ["soilrn", "soillhf", "soilshf", "salb"])
    return nothing
end
function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.Soil.Runoff.SurfaceRunoff,
)
    append!(diagnostics, ["sr"])
    return nothing
end
function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.Soil.Runoff.TOPMODELRunoff,
)
    append!(diagnostics, ["sr", "ssr", "sfsat", "sath", "infc"])
    return nothing
end
function add_diagnostics!(
    diagnostics,
    model::Union{SoilCanopyModel, LandModel},
    subcomponent::ClimaLand.PrescribedRadiativeFluxes,
)
    append!(diagnostics, ["rn"])
    return nothing
end
function add_diagnostics!(
    diagnostics,
    model::CanopyModel,
    subcomponent::ClimaLand.PrescribedAtmosphere,
)
    append!(diagnostics, ["ws"])
    return nothing
end

## Possible diagnostics for standalone models
"""
    get_possible_diagnostics(model)

Return a list containing all possible diagnostics for the given model.
See the file `src/diagnostics/land_compute_methods.jl` to see which model
variable(s) each diagnostic comes from.
"""
function get_possible_diagnostics(model::EnergyHydrology)
    diagnostics = [
        "swc",
        "si",
        "sie",
        "tsoil",
        "et",
        "shc",
        "stc",
        "swp",
        "infil",
        "iwc",
        "precip",
        "sdr",
    ]

    # Add diagnostics based on the top boundary condition type and runoff model
    add_diagnostics!(diagnostics, model, model.boundary_conditions.top)
    add_diagnostics!(diagnostics, model, model.boundary_conditions.top.runoff)
    return unique!(diagnostics)
end

function get_possible_diagnostics(model::SoilCO2Model)
    return ["sco2", "hr", "scd", "scms", "so2", "soc"]
end
function get_possible_diagnostics(model::CanopyModel)
    diagnostics = [
        "sif",
        "ra",
        "gs",
        "trans",
        "clhf",
        "cshf",
        "lwp",
        # "fa", # return a Tuple
        "far",
        "lai",
        "msf",
        "rai",
        "sai",
        "gpp",
        "an",
        "rd",
        "vcmax25",
        "nir",
        "anir",
        "rnir",
        "tnir",
        "par",
        "apar",
        "rpar",
        "tpar",
        "lwn",
        "swn",
        "ct",
        "shf",
        "lhf",
        "er",
        "et",
        "airp", # start of driver diagnostics
        "rain",
        "snow",
        "lwd",
        "swd",
    ]

    # Add conditional diagnostics based on atmosphere type
    add_diagnostics!(diagnostics, model, model.boundary_conditions.atmos)
    return diagnostics
end
function get_possible_diagnostics(model::SnowModel)
    return ["swe", "snd", "snowc"]
end
function get_possible_diagnostics(model::BucketModel)
    return [
        "swa",
        "rn",
        "tsfc",
        "lhf",
        "shf",
        "swu",
        "lwu",
        "vflux",
        "tsoil",
        "wsoil",
        "wsfc",
        "ssfc",
    ]
end

## Possible diagnostics for integrated models
"""
    get_component_diagnostics(model::AbstractLandModel, diagnostics_function::Function)

Helper function for integrated models that combines the sets of diagnostics
for all of its component models.

Based on the input `diagnostics_function`, this function can be used to return either
all possible diagnostics or the short diagnostics of the component models.
`diagnostics_function` should be a Function that takes in a model as its single argument.
This may be `get_possible_diagnostics` or `get_short_diagnostics`.
"""
function get_component_diagnostics(
    model::AbstractLandModel,
    diagnostics_function::Function,
)
    diagnostics = String[]
    @assert diagnostics_function in
            (get_possible_diagnostics, get_short_diagnostics) "Invalid diagnostics_function $(diagnostics_function)"
    # Add diagnostics for each component of the integrated model
    append!(
        diagnostics,
        map(
            component_name ->
                diagnostics_function(getproperty(model, component_name)),
            ClimaLand.land_components(model),
        )...,
    )
    return diagnostics
end
function get_possible_diagnostics(model::SoilCanopyModel)
    component_diagnostics =
        get_component_diagnostics(model, get_possible_diagnostics)

    # Add conditional diagnostics based on radiative forcing type
    add_diagnostics!(
        component_diagnostics,
        model,
        model.canopy.boundary_conditions.radiation,
    )

    additional_diagnostics = ["swa", "swu", "lwu", "infil"]

    return unique!(append!(component_diagnostics, additional_diagnostics))
end
function get_possible_diagnostics(model::LandModel)
    component_diagnostics =
        get_component_diagnostics(model, get_possible_diagnostics)

    # Add conditional diagnostics based on radiative forcing type
    add_diagnostics!(
        component_diagnostics,
        model,
        model.canopy.boundary_conditions.radiation,
    )

    additional_diagnostics = ["swa", "swu", "lwu", "tair", "precip"]

    return unique!(append!(component_diagnostics, additional_diagnostics))
end

"""
    get_short_diagnostics(model)

Return a shortened list of highlighted diagnostics for the given model.
"""
function get_short_diagnostics(model::EnergyHydrology)
    return ["swc", "si", "sie", "tsoil", "et"]
end
function get_short_diagnostics(model::SoilCO2Model)
    return ["sco2"]
end
function get_short_diagnostics(model::CanopyModel)
    return ["gpp", "ct", "lai", "trans", "er", "sif"]
end
function get_short_diagnostics(model::SnowModel)
    return get_possible_diagnostics(model)
end
function get_short_diagnostics(model::SoilCanopyModel)
    component_diagnostics =
        get_component_diagnostics(model, get_short_diagnostics)

    # Add conditional diagnostics based on soil runoff type, since this
    # isn't done for the soil short diagnostics.
    add_diagnostics!(
        component_diagnostics,
        model.soil,
        model.soil.boundary_conditions.top.runoff,
    )
    additional_diagnostics = ["swa", "lwu", "swu"]
    return unique!(append!(component_diagnostics, additional_diagnostics))
end
function get_short_diagnostics(model::LandModel)
    component_diagnostics =
        get_component_diagnostics(model, get_short_diagnostics)
    additional_diagnostics = [
        "swa",
        "lwu",
        "swu",
        "shf",
        "swd",
        "lwd",
        "tair",
        "precip",
        "lhf",
        "msf",
        "lwp",
        "iwc",
        "sdr",
    ]

    # Add conditional diagnostics based on soil runoff type, since this
    # isn't done for the soil short diagnostics.
    add_diagnostics!(
        component_diagnostics,
        model.soil,
        model.soil.boundary_conditions.top.runoff,
    )
    return unique!(append!(component_diagnostics, additional_diagnostics))
end
function get_short_diagnostics(model::BucketModel)
    return get_possible_diagnostics(model)
end
