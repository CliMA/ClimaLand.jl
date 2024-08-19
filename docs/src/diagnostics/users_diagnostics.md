# Using ClimaLand Diagnostics when running a simulation

When running a ClimaLand simulations, you have multiple options on how to write the outputs of that simulation.
You may want all variables, or just a selected few.
You may want instantaneous values, at the highest temporal and spatial resolution, or you may want to get averages at hourly or monthly time scale, and integrate in space
(for example soil moisture from 0 to 1 meter depth).
You may want to get more specific reductions, such as 10 days maximums, or compute a new variables that is a function of others.
You may want to get your outputs in memory in a Julia Dict, or write them in a NetCDF file.

This is where ClimaLand Diagnostics comes in for users.

In this documentation page, we first explain how to use default diagnostics and what are the defaults, and then explain how to define your own diagnostics for more advanced users.

# Default Diagnostics

Once you have defined your model and are ready to run a simulation, and after adding ClimaDiagnostics (using ClimaDiagnostics),
 you can add default diagnostics to it by doing the following steps:

1. define an output folder

```
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path("base_output_dir/")
```

2. define a space

Your diagnostics will be written in time and space. These may be defined in your model, but usually land model space is a sphere with no vertical dimension.
You may have variables varying with soil depth, and so you will need:

```
space = bucket_domain.space.subsurface
```

3. define your writter

Your diagnostics will be written in a Julia Dict or a netcdf file, for example. This is up to you. For a netcdf file, you define your writter like this:

```
nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(space, output_dir)
```

providing the space and output_dir defined in steps 1. and 2.

4. make your diagnostics on your model, using your writter, and define a callback

Now that you defined your model and your writter, you can create a callback function to be called when solving your model. For example:

```
diags = ClimaLand.default_diagnostics(model, 1.0, reference_date; output_writer = nc_writer)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

sol = SciMLBase.solve(prob, ode_algo; dt = Δt, callback = diag_cb)
```

Your diagnostics have now been written in netcdf files in your output folder.

Note that by default, `default_diagnostics` assign two optional kwargs: `output_vars = :long` and `average_period` = :daily.
`output_vars = :long` will write all available diagnostics, whereas `output_vars = :short` will only write essentials diagnostics.
`average_period` defines the period over which diagnostics are averaged, it can be set to `:hourly`, `:daily` and `:monthly`.

# Custom Diagnostics

When defining a custom diagnostic, follow these steps:
 - Define how to compute your diagnostic variable from your model state and cache.
 For example, let's say you want the bowen ratio (ratio between sensible heat and latent heat) in the Bucket model.
 ```
 function compute_bowen_ratio!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.shf / p.bucket.turbulent_fluxes.lhf)
    else
        out .= p.bucket.turbulent_fluxes.shf / p.bucket.turbulent_fluxes.lhf
    end
end
 ```
 - Add that diagnostic variable to your list of variables
 ```
 add_diagnostic_variable!(
        short_name = "bor",
        long_name = "Bowen ratio",
        standard_name = "bowen_ratio",
        units = "",
        compute! = (out, Y, p, t) -> compute_bowen_ratio!(out, Y, p, t, land_model),
    )
 ```
 - Define how to schedule your variables. For example, you want the seasonal maximum of your variables, where season is defined as 90 days.
 ```
 seasonal_maxs(short_names...; output_writer, t_start) = common_diagnostics(
    90 * 24 * 60 * 60 * one(t_start),
    max,
    output_writer,
    t_start,
    short_names...,
)
 ```
  - Define a function to return your `ScheduledDiagnostics`
  ```
  function default_diagnostics(land_model::BucketModel, t_start; output_writer)

    define_diagnostics!(land_model)

	add_diagnostic_variable!(
        short_name = "bor",
        long_name = "Bowen ratio",
        standard_name = "bowen_ratio",
        units = "",
        compute! = (out, Y, p, t) -> compute_bowen_ratio!(out, Y, p, t, land_model),
    )

    my_custom_diagnostics = [
        "lhf",
        "shf",
        "bor",
    ]

    my_custom_outputs =
        seasonal_maxs(my_custom_diagnostics...; output_writer, t_start)
    return [my_custom_outputs...]
end
  ```
