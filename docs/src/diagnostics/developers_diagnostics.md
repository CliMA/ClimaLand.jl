# ClimaLand Diagnostics: why and how

ClimaLand simulations generates variables in the integrator state and cache at each time step.
A user will need to use these variables in some form, i.e., access them from a file that contains variables at a given temporal and spatial resolution.
The user will also want to retrieve metadata about those variables, such as name and units.
This is where ClimaLand diagnostics comes in, it writes simulations variables (in a file, such as NetCDF or HDF5, or in Julia Dict), at a specified spatio-temporal reduction
(e.g., hourly averages, monthly max, instantaneous, integrated through soil depth...), along with metadata (e.g., soil temperature short name is t_soil, expressed in "K" units).
We want to provide users with default options, but also the possibility to define their own variables and reductions.

Internally, this is done by using the [`ClimaDiagnostics.jl`](https://github.com/CliMA/ClimaDiagnostics.jl) package, that provides the functionality to produce a
[`ClimaLand.Diagnostics`](https://github.com/CliMA/ClimaLand.jl/tree/main/src/Diagnostics/Diagnostics.jl) module in the src/Diagnostics.jl folder. In this folder,
 - `Diagnostics.jl` defines the module,
 - `diagnostic.jl` defines `ALL_DIAGNOSTICS`, a Dict containing all diagnostics variables defined in `define_diagnostics.jl`, it also defines the function 
 `add_diagnostic_variable!` which defines a method to add diagnostic variables to ALL_DIAGNOSTICS, finally it contains a function `get_diagnostic_variable` which returns a
 `DiagnosticVariable` from its `short_name`, if it exists. 
 - `define_diagnostics.jl`, mentioned above, contains a function `define_diagnostics!(land_model)` which contains all default diagnostic variables by calling.
 `add_diagnostic_variable!`, and dispatch off the type of land_model to define how to compute a diagnostic (for example, surface temperature is computed in `p.bucket.T_sfc` in the bucket model).
 - compute methods are defined in a separate file, for example, `bucket_compute_methods.jl`. 
 - `standard_diagnostic_frequencies.jl` defines standard functions to schedule diagnostics, for example, hourly average or monthly max, these functions are called on a list of diagnostic variables. As developers, we can add more standard functions that users may want to have access to easily in this file. 
 - `default_diagnostics.jl` defines default diagnostics functions to use on a model simulation. For example, `default_diagnostics(land_model::BucketModel, t_start; output_writer)`.
 will return a `ScheduledDiagnostics` that computes hourly averages for all Bucket variables, along with their metadata, ready to be written on a NetCDF file when running a Bucket simulation.

The following section give more details on these functions, along with examples. As developers, we want to extand these functionality as ClimaLand progresses.

# Compute methods

Each model defines all its compute methods in a file (bucket_compute_methods.jl for the bucket model, for example). 
The structure of a diagnostic variable compute method is, for example:
```
function compute_albedo!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.α_sfc)
    else
        out .= p.bucket.α_sfc
    end
end
```

It defines how to access your diagnostic (here, p.bucket.α_sfc), in your model type (here, ::BucketModel).
Note that, as explained in the [ClimaDiagnostics.jl documentation](https://clima.github.io/ClimaDiagnostics.jl/dev/user_guide/), `out` will probably not be needed in the future.

We also define helper functions returning error messages if a user tries to compute a diagnostic variable that doesn't exist in their model type. 

```
error_diagnostic_variable(variable, land_model::T) where {T} =
    error("Cannot compute $variable with model = $T")

compute_albedo!(_, _, _, _, land_model) =
    error_diagnostic_variable("albedo", land_model)
```

# Define diagnostics

Once the compute functions have been defined, they are added to `define_diagnostics!(land_model)`, which adds diagnostics variables to ALL_DIAGNOSTICS dict,
defined in diagnostic.jl. In these functions, you also define a `short_name`, `long_name`, `standard_name`, `units` and `comment`. For example:

```
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
```
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
        "t",
        "w",
        "ws",
        "sigmas",
    ]

    default_outputs =
        hourly_averages(bucket_diagnostics...; output_writer, t_start)
    return [default_outputs...]
end
```

is the default for the BucketModel, it will return hourly averages for the variables listed in `bucket_diagnostics` (which are all variables in the BucketModel). 

# Standard diagnostic frequencies

We defined some functions of diagnostic schedule that may often be used in `standard_diagnostic_frequencies.jl`, for example

```
hourly_averages(short_names...; output_writer, t_start) = common_diagnostics(
    60 * 60 * one(t_start),
    (+),
    output_writer,
    t_start,
    short_names...;
    pre_output_hook! = average_pre_output_hook!,
)
```

will return a list of `ScheduledDiagnostics` that compute the hourly average for the given variables listed in `short_names`.
We also, so far, provide functions for mins, maxs and averages aggregated monthly, over ten days, daily, and hourly.
As a developer, you may want to add more standard diagnostics here.
