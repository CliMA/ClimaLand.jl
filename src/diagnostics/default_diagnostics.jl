export default_diagnostics

using Dates

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

Helper function to define functions like `daily_max`.
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
    output_writer = NetCDFWriter(
        default_diagnostic_domain,
        outdir;
        start_date,
        sync_schedule = EveryCalendarDtSchedule(Dates.Month(6); start_date),
    )
    return default_diagnostics(model, start_date; output_writer)
end

"""
    default_diagnostics(model::AbstractModel{FT}, start_date; output_writer, output_vars = :short, average_period = :monthly)

For a general AbstractModel, we need a specification of output_vars to determine which diagnostics to output.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.

This method can be extended for any model that extends `get_possible_diagnostics` and `get_short_diagnostics`.
Note that `EnergyHydrology` has a specialized method that handles conservation diagnostics.

Please see the method `get_possible_diagnostics` for the list of available diagnostics for each model.
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
    ) where {FT}

Define a method specific to the EnergyHydrology model so that we can
handle conservation diagnostics specially.

The input `output_vars` can have 3 values:
- `:long` - all diagnostics are output
- `:short` - a short list of diagnostics is output
- `_::Vector{String}` - a user-defined list of diagnostics is output

If a user-defined list is provided for `output_vars`, it must be a vector of strings that are
valid short names of diagnostics for the model.

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
    get_possible_diagnostics(model)

Return a list containing all possible diagnostics for the given model.
See the file `src/diagnostics/land_compute_methods.jl` to see which model
variable(s) each diagnostic comes from.
"""
function get_possible_diagnostics(model::EnergyHydrology{FT}) where {FT}
    return ["swc", "si", "sie", "tsoil", "et"]
end
function get_possible_diagnostics(model::CanopyModel{FT}) where {FT}
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
function get_possible_diagnostics(model::SoilCanopyModel{FT}) where {FT}
    return [
        "swa",
        "sif", # start of canopy diagnostics
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
        "soc", # start of driver diagnostics
        "airp",
        "rain",
        "lwd",
        "swd",
        "snow",
        "qsfc",
        "ws",
        "infil", # start of soil diagnostics
        "shc",
        "stc",
        "swp",
        "soilrn",
        "tsoil",
        "soillhf",
        "soilshf",
        "hr",
        "scd",
        "scms",
        "sco2",
        "swc",
        # "pwc", # return a Tuple
        "si",
        "sie",
        "swu",
        "lwu",
        "er",
        "et",
        "sr",
        "rn",
        "lhf",
        "shf",
        "salb",
    ]
end
function get_possible_diagnostics(model::LandModel{FT}) where {FT}
    return [
        "swa",
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
        "soc",
        "airp",
        "rain",
        "lwd",
        "swd",
        "snow",
        "qsfc",
        "infil",
        "shc",
        "stc",
        "swp",
        "soilrn",
        "tsoil",
        "soillhf",
        "soilshf",
        "hr",
        "scd",
        "scms",
        "ct",
        "sco2",
        "swc",
        # "pwc", # return a Tuple
        "si",
        "sie",
        "swu",
        "lwu",
        "er",
        "et",
        "sr",
        "swe",
        "snd",
        "rn",
        "lhf",
        "shf",
        "iwc",
        "snowc",
        "tair",
        "precip",
    ]
end
function get_possible_diagnostics(model::BucketModel{FT}) where {FT}
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

"""
    get_short_diagnostics(model)

Return a shortened list of highlighted diagnostics for the given model.
"""
function get_short_diagnostics(model::EnergyHydrology{FT}) where {FT}
    return get_possible_diagnostics(model)
end
function get_short_diagnostics(model::CanopyModel{FT}) where {FT}
    return ["gpp", "ct", "lai", "trans", "er", "sif"]
end
function get_short_diagnostics(model::SoilCanopyModel{FT}) where {FT}
    return [
        "gpp",
        "ct",
        "lai",
        "sco2",
        "swc",
        "si",
        "swa",
        "lwu",
        "et",
        "er",
        "sr",
        "sif",
    ]
end
function get_short_diagnostics(model::LandModel{FT}) where {FT}
    return [
        "gpp",
        "swc",
        "si",
        "sie",
        "swu",
        "lwu",
        "et",
        "sr",
        "ssr",
        "swe",
        "shf",
        "lhf",
        "trans",
        "msf",
        "lwp",
        "tsoil",
        "lai",
        "iwc",
        "swd",
        "lwd",
        "snowc",
        "tair",
        "precip",
    ]
end
function get_short_diagnostics(model::BucketModel{FT}) where {FT}
    return get_possible_diagnostics(model)
end
