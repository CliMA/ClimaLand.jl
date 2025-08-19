# Using ClimaLand Diagnostics to save simulation output

When running a ClimaLand simulations, you have multiple options on how to write the outputs of that simulation.
You may want all variables, or just a selected few.
You may want instantaneous values, at the highest temporal and spatial resolution, or you may want to get averages at hourly or monthly time scale, and integrate in space
(for example soil moisture from 0 to 1 meter depth).
You may want to get more specific reductions, such as 10 days maximums, or compute a new variables as a function of others.
You may want to get your outputs in memory in a Julia dictionary, or write them in a NetCDF file.

This is where ClimaLand Diagnostics comes in for users.

In this documentation page, we first explain how to use default diagnostics and what are the defaults, and then explain how to define your own diagnostics for more advanced users.

## Default Diagnostics

Diagnostics refer to output saved during the simulation, which may be prognostic or diagnostic variables.
Note that this is different from checkpointing the simulation, which has specific requirements.
For information about checkpointing and restarting simulations, please see the page titled
[Restarting Simulations](@ref).

The main user-facing function in the ClimaLand.Diagnostics module is `default_diagnostics`. This function defines
what diagnostic variables to compute by default for a specific model, and
on what schedule (for example, hourly average).

`default_diagnostics` takes in the following arguments:
- `model`: The ClimaLand model to generate diagnostics for. Currently the following models support diagnostics: `CanopyModel`, `EnergyHydrologyModel`, `SoilCanopyModel`, `LandModel`, `BucketModel`
- `start_date`: The start date of the simulation.
- `output_writer`: An object of type `ClimaDiagnostics.AbstractWriter`. Specifically this may be a `NetCDFWriter`, `HDF5Writer`, or `DictWriter`, which save output to a NetCDF file, HDF5 file, or dictionary in Julia memory, respectively. For more details about the diagnostics writers, please see the ClimaDiagnostics.jl documentation.
- `output_vars`: This argument may be `:short` to output a default list of variables defined for each model, `:long` to output all
available variables for the model, or a user-defined list of variable "short names".
- `average_period`: How often to compute and save the average of the diagnostics.
- `dt`: Simulation timestep; only required for instantaneous diagnostics.

## Example: diagnostics for a `CanopyModel`

The following code sets up default short diagnostics to be averaged hourly and written in memory as a Julia dictionary:
```julia
diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    canopy,
    start_date;
    output_vars = :short,
    output_writer = diag_writer,
    average_period = :hourly,
);
```

To instead output a list of specific diagnostics, you can change the value of `output_vars`.
For example, to output gross primary productivity (GPP) and transpiration you would do the following:
```julia
diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    canopy,
    start_date;
    output_vars = ["gpp", "trans"],
    output_writer = diag_writer,
    average_period = :hourly,
);
```
A description of available diagnostics and their short names can be found on the [Available diagnostic variables](@ref) documentation page.

To write the diagnostics to a NetCDF file instead of saving it in memory, the `diag_writer` should be constructed as a `NetCDFWriter` and then passed to `default_diagnostics` as before:
```julia
outdir = "my_output_dir"
nc_writer =
    ClimaDiagnostics.Writers.NetCDFWriter(space, outdir; start_date)
```

!!! note

    The `NetCDFWriter` currently writes to file for each diagnostic output, which can be quite slow when saving variables at every step.
    On the other hand, the `DictWriter` saves output to memory which may be too large in global runs, so
    `DictWriter` usually should not be used for global runs.
    In general, we recommend using `DictWriter` for column simulations, and `NetCDFWriter` for global simulations.

## Diagnostics output naming and format

Diagnostics are typically named using the format `$(short_name)_$(period)_$(reduction)`.
For example, with the NetCDFWriter, hourly-averaged GPP would be saved in an output file titled `gpp_1h_average.nc`.

The specific output format depends on which output writer is being used; for more details,
please see the [ClimaDiagnostics documentation](https://clima.github.io/ClimaDiagnostics.jl/stable/writers/).

## Adding new diagnostics

To define a new, custom diagnostic, you must follow these steps:
- specify how to compute the diagnostic
- manually define the diagnostic via `add_diagnostic_variable!`
- add the diagnostic to the list of possible diagnostics for the relevant model(s)

### Define how to compute your diagnostic variable from your model state and cache.

These functions are defined in `src/diagnostics/land_compute_methods.jl` and must be named
in the format `compute_[standard_name]!`. Be sure to write method(s) for each model you want
this diagnostic to be available for.
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
However, it is better to use that macro only if you are retrieving an existing variable, such as latent heat flux,
rather than computing one, like the Bowen ratio above. For example,

```Julia
@diagnostic_compute "latent_heat_flux" BucketModel p.bucket.turbulent_fluxes.lhf
```

### Add that diagnostic(s) variable to your list of variables

These functions are defined in `src/diagnostics/define_diagnostics.jl`.

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

### Define how often to output your variables, and how to combine them.

We have many common output periods and reductions defined in `src/diagnostics/standard_diagnostic_frequencies.jl`,
including averages, maximums, and minimums over different periods of time, and instantaneous outputs.

You can also define your own output frequencies and reductions. For example, if you want the seasonal maximum of your
variables, where season is defined as 90 days, you could add the following function.

```Julia
seasonal_maxs(FT, short_names...; output_writer) = common_diagnostics(
    FT(90 * 24 * 60 * 60), # 90 days in seconds
    max,
    output_writer,
    nothing, # start_date
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
