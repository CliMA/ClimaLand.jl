# ClimaLand Diagnostics: why and how

ClimaLand simulations generates variables in the integrator state (Y) and cache (p) at each time step.
A user will need to use these variables in some form, i.e., access them from a file that contains variables at a given temporal and spatial resolution.
The user will also want to retrieve metadata about those variables, such as name and units.
This is where ClimaLand diagnostics comes in, it writes simulations variables (in a file, such as NetCDF or HDF5, or in Julia Dict), at a specified spatio-temporal reduction
(e.g., hourly averages, monthly max, instantaneous, integrated through soil depth...), along with metadata (e.g., soil temperature short name is t\_soil, expressed in "K" units).
We want to provide users with default options, but also the possibility to define their own variables and reductions.

Internally, this is done by using the [`ClimaDiagnostics.jl`](https://github.com/CliMA/ClimaDiagnostics.jl) package, that provides the functionality to produce a
[`ClimaLand.Diagnostics`](https://github.com/CliMA/ClimaLand.jl/tree/main/src/Diagnostics/Diagnostics.jl) module in the src/Diagnostics.jl folder. In this folder,
 - `Diagnostics.jl` defines the module,
 - `diagnostic.jl` defines `ALL_DIAGNOSTICS`, a Dict containing all diagnostics variables defined in `define_diagnostics.jl`,
along with a boolean indicating if each diagnostic variable is a default. It also defines the function
 `add_diagnostic_variable!` which defines a method to add diagnostic variables to ALL\_DIAGNOSTICS, finally it contains a function `get_diagnostic_variable` which returns a
 `DiagnosticVariable` from its `short_name`, if it exists.
 - `define_diagnostics.jl`, mentioned above, contains a function `define_diagnostics!(land_model)` which contains all default diagnostic variables by calling.
 `add_diagnostic_variable!`, and dispatch off the type of land\_model to define how to compute a diagnostic (for example, surface temperature is computed in `p.bucket.T_sfc` in the bucket model).
 - compute methods are defined in a separate file, for example, `bucket_compute_methods.jl`.
 - `standard_diagnostic_frequencies.jl` defines standard functions to schedule diagnostics, for example, hourly average or monthly max, these functions are called on a list of diagnostic variables. As developers, we can add more standard functions that users may want to have access to easily in this file.
 - `default_diagnostics.jl` defines default diagnostics functions to use on a model simulation. For example, `default_diagnostics(land_model::BucketModel, output_writer)`.
 will return a `ScheduledDiagnostics` that computes hourly averages for all Bucket variables, along with their metadata, ready to be written on a NetCDF file when running a Bucket simulation.

The following section give more details on these functions, along with examples. As developers, we want to extand these functionality as ClimaLand progresses.

# Compute methods

Each model defines all its compute methods in a file (bucket\_compute\_methods.jl for the bucket model, for example).
The structure of a diagnostic variable compute method is, for example:
```Julia
@with_error function compute_albedo!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.α_sfc)
    else
        out .= p.bucket.α_sfc
    end
end
```

It defines how to access your diagnostic (here, p.bucket.α\_sfc) with the land\_model `BucketModel`.
Note that you can also use the @diagnostic\_compute macro to do the same thing:

```Julia
@diagnostic_compute "albedo" BucketModel p.bucket.α\_sfc
```

The `@with_error` macro define helper functions returning error messages if a user tries to compute a diagnostic variable that doesn't exist in their model type.

# Define diagnostics

Once the compute functions have been defined, they are added to `define_diagnostics!(land_model)`, which adds diagnostics variables to ALL\_DIAGNOSTICS dict,
defined in diagnostic.jl. In these functions, you also define a `short_name`, `long_name`, `standard_name`, `units` and `comment`. For example:

```Julia
add_diagnostic_variable!(
        short_name = "alpha",
        long_name = "Albedo",
        standard_name = "albedo",
        units = "",
        compute! = (out, Y, p, t) -> compute_albedo!(out, Y, p, t, land_model),
    )
```

# Default diagnostics

For each model, we define a function `default_diagnostics` which will define what diagnostic variables to compute by default for a specific model, and
on what schedule (for example, hourly average). For example,

```Julia
function default_diagnostics(land_model::BucketModel{FT}; output_writer) where {FT}

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
        "t",
        "w",
        "ws",
        "sigmas",
    ]

    default_outputs =
        hourly_averages(FT, bucket_diagnostics...; output_writer)
    return [default_outputs...]
end
```

is the default for the BucketModel, it will return hourly averages for the variables listed in `bucket_diagnostics` (which are all variables in the BucketModel).

For SoilCanopyModel, LandModel, and SoilModel, we added two keyword arguments: `output_vars` (can be :long or :short) and `average_period` (can be :hourly, :daily, or :monthly).
If `output_vars = :long` (the default), then `soilcanopy_diagnostics` is an Array of all short\_name, if `output_vars = :short`, then we have a subset of all short\_name: `soilcanopy_diagnostics = ["gpp", "ct", "lai", "sco2", "swc", "si", "swa", "lwu", "et", "er", "sr", "sif"]`.
If `average_period = :hourly`, `default_outputs` calls `hourly_averages`, et cetera.

# Standard diagnostic frequencies

We defined some functions of diagnostic schedule that may often be used in `standard_diagnostic_frequencies.jl`, for example

```Julia
hourly_averages(FT, short_names...; output_writer) = common_diagnostics(
    60 * 60 * one(FT),
    (+),
    output_writer,
    nothing, # start_date
    short_names...;
    pre_output_hook! = average_pre_output_hook!,
)
```

will return a list of `ScheduledDiagnostics` that compute the hourly average for the given variables listed in `short_names`.
We also, so far, provide functions for mins, maxs and averages aggregated monthly, over ten days, daily, and hourly.
As a developer, you may want to add more standard diagnostics here.
