# Using ClimaLand Diagnostics when running a simulation

When running a ClimaLand simulations, you have multiple options on how to write the outputs of that simulation.
You may want all variables, or just a selected few.
You may want instantaneous values, at the highest temporal and spatial resolution, or you may want to get averages at hourly or monthly time scale, and integrate in space
(for example soil moisture from 0 to 1 meter depth).
You may want to get more specific reductions, such as 10 days maximums, or compute a new variables that is a function of others.
You may want to get your outputs in memory in a Julia Dict, or write them in a NetCDF file.

This is where ClimaLand Diagnostics comes in for users.

In this documentation page, we first explain how to use default diagnostics and what are the defaults, and then explain how to define your own diagnostics for more advanced users.

## Default Diagnostics

Once you have defined your model and are ready to run a simulation, and after adding ClimaDiagnostics (using ClimaDiagnostics),
 you can add default diagnostics to it by doing the following steps:

### define an output folder

```Julia
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path("base_output_dir/")
```

### define a space

Your diagnostics will be written in time and space. These may be defined in your model, but usually land model space is a sphere with no vertical dimension.
You may have variables varying with soil depth, and so you will need:

```Julia
space = bucket_domain.space.subsurface
```

### define your writter

Your diagnostics will be written in a Julia Dict or a netcdf file, for example. This is up to you. For a netcdf file, you define your writter like this:

```Julia
nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(space, output_dir; start_date)
```

providing the space and output_dir defined in steps 1. and 2.

### make your diagnostics on your model, using your writter, and define a callback

Now that you defined your model and your writter, you can create a callback function to be called when solving your model. For example:

```Julia
t0 = 0 # the start date of your simulation

start_date = DateTime(2024) # start_date is the DateTime of your start date

diags = ClimaLand.default_diagnostics(model, start_date; output_writer = nc_writer)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

sol = SciMLBase.solve(prob, ode_algo; dt = Δt, callback = diag_cb)
```

Your diagnostics have now been written in netcdf files in your output folder.

Note that by default, `default_diagnostics` assign two optional kwargs: `output_vars = :long` and `average_period` = :daily.
`output_vars = :long` will write all available diagnostics, whereas `output_vars = :short` will only write essentials diagnostics.
`average_period` defines the period over which diagnostics are averaged, it can be set to `:hourly`, `:daily` and `:monthly`.

## Custom Diagnostics

When defining a custom diagnostic, follow these steps:

### Define how to compute your diagnostic variable from your model state and cache.

For example, let's say you want the bowen ratio (ratio between sensible heat and latent heat) in the Bucket model.

```Julia
function compute_bowen_ratio!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.shf / p.bucket.turbulent_fluxes.lhf)
    else
        out .= p.bucket.turbulent_fluxes.shf / p.bucket.turbulent_fluxes.lhf
    end
end
```

Or, for convenience, you can use the `@diagnostic_compute` macro which generates the same function.
However, it is better to use that macro only if you are getting a defined variable, such as latent heat flux.
(without an operation like the bowen ratio above). For example,

```Julia
@diagnostic_compute "latent_heat_flux" BucketModel p.bucket.turbulent_fluxes.lhf
```

### Add that diagnostic(s) variable to your list of variables

```Julia
 add_diagnostic_variable!(
    short_name = "bor",
    long_name = "Bowen ratio",
    standard_name = "bowen_ratio",
    units = "",
    comments = "Ratio of sensible to latent heat flux.",
    compute! = (out, Y, p, t) -> compute_bowen_ratio!(out, Y, p, t, land_model),
)

add_diagnostic_variable!(
    short_name = "lhf",
    long_name = "Latent Heat Flux",
    standard_name = "latent_heat_flux",
    units = "W m^-2",
    comments = "Exchange of energy at the land-atmosphere interface due to water evaporation or sublimation.",
    compute! = (out, Y, p, t) ->
    compute_latent_heat_flux!(out, Y, p, t, land_model),
)
```

### Define how to schedule your variables. For example, you want the seasonal maximum of your variables, where season is defined as 90 days.

```Julia
seasonal_maxs(FT, short_names...; output_writer) = common_diagnostics(
    90 * 24 * 60 * 60 * one(FT),
    max,
    output_writer,
    short_names...,
)
```

### Define a function to return your `ScheduledDiagnostics`

Now, you can call your schedule with your variables.

```Julia
my_custom_diagnostics = ["lhf", "bor"]

diags = seasonal_maxs(FT, my_custom_diagnostics...; output_writer)
```

### Analyze your simulation output

Once you've run your simulation and created an output folder (e.g., output\_dir) with diagnostics, you can use [ClimaAnalysis](https://github.com/CliMA/ClimaAnalysis.jl)
to access and analyze your data. For in depth documentation about ClimaAnalysis, see its [documentation](https://clima.github.io/ClimaAnalysis.jl/stable/).

Here is an example of how to plot a variable:

```Julia
import ClimaAnalysis

import ClimaAnalysis.Visualize as viz

import CairoMakie # the plotting package used by ClimaAnalysis

simdir = ClimaAnalysis.SimDir(output_dir) # where output_dir is where you saved your diagnostics.

var = get(simdir; "lhf") # assuming lhf, latent_heat_flux used as an example above, is one of your diagnostics variables.

fig = CairoMakie.Figure() # creates an empty figure object

viz.plot!(fig, var) # creates an axis inside fig, and plot your var in it.

CairoMakie.save(fig) # saves the figure in current working directory
```

