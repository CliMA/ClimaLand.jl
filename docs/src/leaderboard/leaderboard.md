# Leaderboard

## Long run

### Add a new variable to compare against observations

The infrastructure to compute errors against observations is in the
`leaderboard` folder. This folder contains two files: `data_sources.jl`,
responsible for loading and preprocessing variables of interest, and
`leaderboard.jl`, which computes error and draws plots. To add a new variable to
the comparison, you modify `data_sources.jl`.

The file `data_sources.jl` is organized around two data loader structs,
`ILAMBDataLoader` and `ERA5DataLoader`, both subtypes of `AbstractDataLoader`.
Each loader exposes a `Base.get` method that returns a preprocessed
[`OutputVar`](https://clima.github.io/ClimaAnalysis.jl/stable/var/) and a
`preprocess` method that dispatches on the variable short name via
`Val{:shortname}`.

### Computation

As of now, the leaderboard produces bias plots with the global bias and global
root mean squared error (RMSE). These quantities are computed for each month
with the first year of the simulation not considered as that is the spinup time.
The start date of the simulation is 2008 which means that only the year 2009 is
used to compare against observational data.

### Add a new variable to the bias plots

To add a new variable you need to touch four places in `data_sources.jl`:
`_preprocess_sim_var`, the appropriate data loader constructor and `preprocess`
method, `get_mask_dict`, and `get_compare_vars_biases_plot_extrema`.

**1. Preprocess the simulation variable**

Add a `_preprocess_sim_var` method dispatching on `Val{:shortname}` to convert
units or perform any other preprocessing on the simulation output.

```julia
function _preprocess_sim_var(var, ::Val{:et})
    (ClimaAnalysis.units(var) == "kg m^-2 s^-1") && (
        var = ClimaAnalysis.convert_units(
            var,
            "mm / day",
            conversion_function = units -> units * 86400.0,
        )
    )
    return var
end
```

`preprocess_sim_var` is called automatically by `leaderboard.jl` on the
`OutputVar` loaded from the simulation diagnostics.

**2. Add the observational variable to a data loader**

Observational data is provided by `ILAMBDataLoader` (ILAMB datasets) or
`ERA5DataLoader` (ERA5 monthly averages). To add a new variable to
`ILAMBDataLoader`, register the NetCDF file in the constructor and add a
`preprocess` method for the new short name.

In the `ILAMBDataLoader()` constructor:

```julia
new_var_filepath = ClimaLand.Artifacts.ilamb_dataset_path("new_var_dataset.nc")
ClimaAnalysis.add_file!(catalog, new_var_filepath, "nc_varname" => "new_var")
```

Then add a preprocessing method:

```julia
function preprocess(::ILAMBDataLoader, var, ::Val{:new_var})
    replace!(var, missing => NaN)
    return _preprocess_var(var)
end
```

`Base.get(loader, short_name)` will call `ClimaAnalysis.transform_dates!`
which shifts the times to the first day of the month before dispatching to
`preprocess`, so date alignment is handled automatically.

!!! tip "Preprocessing"
    Observational and simulation data should use the same units and time
    conventions. We follow the convention that a monthly average is associated
    with the first day of that month. The function `transform_dates!` is applied
    automatically in `Base.get` for both loaders, so you only need to handle
    unit conversion and `missing`/NaN cleanup in `preprocess`.

**3. Add a mask**

Add an entry to `get_mask_dict` for the loader that provides the new variable.
The value is a function that takes `sim_var` and `obs_var` and returns a masking
function. The masking function is used to correctly normalize the global bias
and global RMSE.

```julia
mask_dict["new_var"] =
    (sim_var, obs_var) -> begin
        return ClimaAnalysis.make_lonlat_mask(
            # We do this to get a `OutputVar` with only two dimensions:
            # longitude and latitude
            ClimaAnalysis.slice(
                obs_var,
                time = ClimaAnalysis.times(obs_var) |> first,
            );
            # Any values that are NaN should be 0.0
            set_to_val = isnan,
            true_val = 0.0
        )
    end
```

**4. Set bias plot limits**

Add a key-value pair to `get_compare_vars_biases_plot_extrema` whose value is a
tuple `(lower, upper)` setting the color scale range for the bias plots.

```julia
compare_vars_biases_plot_extrema = Dict(
    "et" => (-2.0, 2.0) .* factor,
    "gpp" => (-6.0, 6.0) .* factor,
    "new_var" => (-10.0, 10.0) .* factor,
    ...
)
```
