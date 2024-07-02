export default_diagnostics

# This file is included by Diagnostics.jl and defines all the defaults for
# various models (e.g., Bucket, SoilCanopyModel). A model here is either a
# standalone (e.g., Bucket) or integrated (e.g., SoilCanopy) model.
#
# If you are developing new models, add your defaults here. If you want to add
# more high level interfaces, add them here. Feel free to include extra files.

# Bucket model

"""
    function common_diagnostics(
                                period,
                                reduction,
                                output_writer,
                                t_start,
                                short_names...;
                                pre_output_hook! = nothing,
                               )

Helper function to define functions like `daily_max`.
"""
function common_diagnostics(
    period,
    reduction,
    output_writer,
    t_start,
    short_names...;
    pre_output_hook! = nothing,
)
    return [
        ScheduledDiagnostic(
            variable = get_diagnostic_variable(short_name),
            compute_schedule_func = EveryStepSchedule(),
            output_schedule_func = EveryDtSchedule(period; t_start),
            reduction_time_func = reduction,
            output_writer = output_writer,
            pre_output_hook! = pre_output_hook!,
        ) for short_name in short_names
    ]
end

include("standard_diagnostic_frequencies.jl")

# Bucket
function default_diagnostics(land_model::BucketModel, t_start; output_writer)

    define_diagnostics!(land_model)

    bucket_diagnostics = [
        "alpha",
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
    ] # TO DO: would it be helpful to return this list?

    default_outputs =
        hourly_averages(bucket_diagnostics...; output_writer, t_start)
    return [default_outputs...]
end

# SoilCanopyModel
# TO DO: short and long list of variables
function default_diagnostics(
    land_model::SoilCanopyModel,
    t_start;
    output_writer,
    output_vars = "long",
)

    define_diagnostics!(land_model)
    
    if output_vars == "short"
        soilcanopy_diagnostics = [
            "rn",
            "lhf",
            "rae",
            "shf",
            "vflux",
            "tsoil",
            "slw",
            "infil",
            "scd",
            "scms",
            "gs",
            "mt",
            "trans",
            "rain", # do we want?
            "an",
            "gpp",
            "rd",
            "vcmax25",
            "par",
            "apar",
            "rpar",
            "tpar",
            "nir",
            "anir",
            "rnir",
            "tnir",
            "swn",
            "lwn",
            "ra",
            "soilco2",
        ]
    elseif output_vars == "short"
        soilcanopy_diagnostics = [
            #"lhf",
            #"shf",
            "gpp",
            #"tsoil",
            #"rn",
            #"trans",
        ]
    end

    default_outputs =
        hourly_averages(soilcanopy_diagnostics...; output_writer, t_start)
    return [default_outputs...]
end
