export default_diagnostics

# This file is included by Diagnostics.jl and defines all the defaults for
# various models (e.g., Bucket, SoilCanopyModel). A model here is either a
# standalone (e.g., Bucket) or integrated (e.g., SoilCanopy) model.
#
# If you are developing new models, add your defaults here. If you want to add
# more high level interfaces, add them here. Feel free to include extra files.

"""
    function common_diagnostics(
                                period,
                                reduction,
                                output_writer,
                                start_date,
                                short_names...;
                                pre_output_hook! = nothing,
                               )

For each variable specified in `short_names`, create a `ClimaDiagnostics.ScheduledDiagnostic`.
Diagnostics are computed at every step, and output at the
specified period with the provided reduction.
The type of `output_writer` determines where the output will be saved to.
This must be a `ClimaDiagnostics.AbstractWriter` object: concrete options include
`NetCDFWriter`, `HDF5Writer`, and `DictWriter`.
"""
function common_diagnostics(
    period,
    reduction,
    output_writer,
    start_date,
    short_names...;
    pre_output_hook! = nothing,
)
    return vcat(
        map(short_names) do short_name
            output_schedule_func =
                period isa Period ?
                EveryCalendarDtSchedule(period; start_date) :
                EveryDtSchedule(period)
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

include("standard_diagnostic_frequencies.jl")

default_diagnostics(
    model::ClimaLand.AbstractModel,
    start_date::ITime{<:Any, <:Any, <:DateTime},
    outdir,
) = default_diagnostics(model, date(start_date), outdir)

# The default diagnostics currently require a start date because they use Dates.Period.
function default_diagnostics(
    model::ClimaLand.AbstractModel,
    start_date::Union{Nothing, ITime{<:Any, <:Any, Nothing}},
    outdir,
)
    @warn "Default diagnostics not available when running without a start date."
    return []
end

function default_diagnostics(model::ClimaLand.AbstractModel, start_date, outdir)
    # a start date is required for default diagnostics
    domain = ClimaLand.get_domain(model)
    default_diagnostic_domain =
        haskey(domain.space, :subsurface) ? domain.space.subsurface :
        domain.space.surface
    output_writer = NetCDFWriter(default_diagnostic_domain, outdir; start_date)
    return default_diagnostics(model, start_date; output_writer)
end

"""
    default_diagnostics(model::AbstractModel{FT}, start_date; output_writer, output_vars = :short, average_period = :monthly, dt = nothing)

Construct a list of `ScheduledDiagnostics` that outputs the given variables at the specified average period.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.
Please see the method `get_possible_diagnostics` for the list of available diagnostics for each model.

`average_period` specifies the frequency at which to average the diagnostics. The following options are currently supported:
    - `:instantaneous` (note: `dt` must be specified in this case)
    - `:halfhourly`
    - `:hourly`
    - `:daily`
    - `:monthly`

This method can be extended for any model that extends `get_possible_diagnostics` and `get_short_diagnostics`.
Note that `EnergyHydrology` has a specialized method that handles conservation diagnostics.
"""
function default_diagnostics(
    model::Union{
        CanopyModel{FT},
        SoilCanopyModel{FT},
        LandModel{FT},
        BucketModel{FT},
    },
    start_date;
    output_writer,
    output_vars = :short,
    average_period = :monthly,
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
        @assert all([var ∈ possible_diags for var in output_vars])
        diagnostics = output_vars
    end


    if average_period == :halfhourly
        default_outputs =
            halfhourly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :hourly
        default_outputs =
            hourly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :daily
        default_outputs =
            daily_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :monthly
        default_outputs =
            monthly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :instantaneous
        @assert !isnothing(dt) "dt must be specified when `average_period = :instantaneous`"
        default_outputs =
            every_dt_inst(FT, dt, diagnostics...; output_writer, start_date)
    else
        @error("Invalid diagnostics average period $(average_period)")
    end
    return [default_outputs...]
end

"""
    default_diagnostics(
        land_model::EnergyHydrology{FT},
        start_date;
        output_writer,
        output_vars = :short,
        average_period = :monthly,
        conservation = false,
        conservation_period = Day(10),
        dt = nothing,
    ) where {FT}

Define a method specific to the EnergyHydrology model so that we can
handle conservation diagnostics specially.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.

The following options for `average_period` are currently supported:
    - `:instantaneous` (note: `dt` must be specified in this case)
    - `:halfhourly`
    - `:hourly`
    - `:daily`
    - `:monthly`

Conservation diagnostics should not be provided as part of the `output_vars` argument,
but rather included by providing `conservation = true`.
Please see the method `get_possible_diagnostics` for the list of available diagnostics.
"""
function default_diagnostics(
    land_model::EnergyHydrology{FT},
    start_date;
    output_writer,
    output_vars = :short,
    average_period = :monthly,
    conservation = false,
    conservation_period = Day(10),
    dt = nothing,
) where {FT}

    define_diagnostics!(land_model)

    possible_diags = get_possible_diagnostics(land_model)
    if output_vars in (:long, :short)
        diagnostics = possible_diags
    else
        @assert typeof(output_vars) <: Vector{String}
        @assert all([var ∈ possible_diags for var in output_vars])
        diagnostics = output_vars
    end

    if average_period == :halfhourly
        default_outputs =
            halfhourly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :hourly
        default_outputs =
            hourly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :daily
        default_outputs =
            daily_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :monthly
        default_outputs =
            monthly_averages(FT, diagnostics...; output_writer, start_date)
    elseif average_period == :instantaneous
        @assert !isnothing(dt) "dt must be specified when `average_period = :instantaneous`"
        default_outputs =
            every_dt_inst(FT, dt, diagnostics...; output_writer, start_date)
    else
        @error("Invalid diagnostics average period $(average_period)")
    end

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
    start_date = nothing;
    output_writer = nothing,
    output_vars = nothing,
    average_period = nothing,
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

## Possible diagnostics for standalone models
"""
    get_possible_diagnostics(model)

Return a list containing all possible diagnostics for the given model.
See the file `src/diagnostics/land_compute_methods.jl` to see which model
variable(s) each diagnostic comes from.
"""
function get_possible_diagnostics(model::EnergyHydrology)
    diagnostics = ["swc", "si", "sie", "tsoil", "et", "shc", "stc", "swp"]

    # Add diagnostics based on the top boundary condition type and runoff model
    add_diagnostics!(diagnostics, model, model.boundary_conditions.top)
    add_diagnostics!(diagnostics, model, model.boundary_conditions.top.runoff)
    return unique!(diagnostics)
end
function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.AtmosDrivenFluxBC,
)
    append!(diagnostics, ["soilrn", "soillhf", "soilshf", "salb"])
end
function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.Soil.Runoff.SurfaceRunoff,
)
    append!(diagnostics, ["sr"])
end
function add_diagnostics!(
    diagnostics,
    model::EnergyHydrology,
    subcomponent::ClimaLand.Soil.Runoff.TOPMODELRunoff,
)
    append!(diagnostics, ["sr", "ssr"])
end

function get_possible_diagnostics(model::SoilCO2Model)
    return ["sco2", "hr", "scd", "scms", "soc"]
end
function get_possible_diagnostics(model::CanopyModel)
    return [
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
        "qsfc",
        "ws",
    ]
end
function get_possible_diagnostics(model::SnowModel)
    return ["swe", "snd", "snowc"]
end
function get_possible_diagnostics(model::BucketModel)
    return [
        "swa",
        "rn",
        "tsfc",
        "qsfc",
        "lhf",
        "rae",
        "shf",
        "vflux",
        "rhosfc",
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
    component_diags = get_component_diagnostics(model, get_possible_diagnostics)
    additional_diagnostics = ["swa", "swu", "lwu", "rn", "infil"]

    return unique!(append!(component_diags, additional_diagnostics))
end
function get_possible_diagnostics(model::LandModel)
    component_diags = get_component_diagnostics(model, get_possible_diagnostics)
    additional_diagnostics =
        ["swa", "swu", "lwu", "rn", "infil", "iwc", "tair", "precip"]

    return unique!(append!(component_diags, additional_diagnostics))
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
