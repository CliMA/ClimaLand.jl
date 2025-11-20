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
 - `diagnostic.jl` defines `ALL_DIAGNOSTICS`, a Dict containing all diagnostics variables defined in `construct_diagnostics.jl`, it also defines the function
 `add_diagnostic_variable!` which defines a method to add diagnostic variables to ALL\_DIAGNOSTICS, finally it contains a function `get_diagnostic_variable` which returns a
 `DiagnosticVariable` from its `short_name`, if it exists.
 - `construct_diagnostics.jl`, mentioned above, contains a function `construct_diagnostics(land_model)` which contains all default diagnostic variables by calling.
 `add_diagnostic_variable!`, and dispatch off the type of land\_model to define how to compute a diagnostic (for example, surface temperature is computed in `p.bucket.T_sfc` in the bucket model).
 - `land_compute_methods.jl` defines how to compute diagnostics for various ClimaLand models.

The following section give more details on these functions, along with examples.

# Default diagnostics

For each model, we define a function `default_diagnostics` which retrieves the set of diagnostic
variables to compute for this model, and over which time period to average the values.

## "Long" diagnostics (`output_vars = :long`)
Each model has a method of the function `get_possible_diagnostics` which contains a
list of all available diagnostics for the model. This is exactly the set of diagnostics
output when `output_vars = :long`.
For standalone models, the list is mostly hardcoded, though it may also depend on
which parameterizations are available in the model. For integrated models, the list
contains all possible diagnostics for each component model, as well as any additional
diagnostics specific to the integrated model.

## "Short" diagnostics (`output_vars = :short`)
Similarly, each model has a method of `get_short_diagnostics`, which is a selected subset
of the overall available diagnostics for the model that gets used when `output_vars = :short`.
As is the case with `get_possible_diagnostics`, the set of variables for integrated models is
based on its component models, and any additional integrated model-specific variables.

# Standard diagnostic frequencies

We have defined functions which compute statistical metrics of the instantaneous diagnostic variables at different standard frequencies. For example,

```Julia
hourly_averages(FT, short_names...; output_writer) = common_diagnostics(
    :hourly,
    :average,
    output_writer,
    nothing, # start_date
    short_names...;
)
```

will return a list of `ScheduledDiagnostics` that compute the hourly average for the given variables listed in `short_names`.
We also, so far, provide functions for mins, maxs and averages aggregated monthly, over ten days, daily, hourly, and half-hourly.
As a developer, you may want to add more standard diagnostics here.

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

Once the compute functions have been defined, they are added to `construct_diagnostics(land_model)`, which adds diagnostics variables to ALL\_DIAGNOSTICS dict,
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
