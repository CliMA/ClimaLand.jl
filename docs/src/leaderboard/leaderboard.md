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

Note that `get_mask_dict(::ERA5DataLoader)` gives every ERA5 variable the ocean
mask, so step 3 is only needed for loaders whose observations have gaps over
land.

## Annual leaderboard columns

`compute_seasonal_leaderboard` draws one row per variable and the columns below.
The names are the short headings drawn on the figure:

- `SIM` (simulation), the simulated annual mean,
- `ANN` (annual bias), its bias against observations,
- `LAT` (zonal mean), the mean along each latitude band of both, banded by the
  standard deviation along that band,
- `MON` (seasonal cycle), their global monthly means, banded by the spread
  across years,
- `IAV` (interannual variability), the area-weighted global annual mean of one
  against the other, one point per year, with the 1:1 line, the least squares
  fit, its slope, r² and `σ_sim / σ_obs`.

Only years covering all twelve months enter `IAV`, since the annual mean of a
partial year is aliased by whichever part of the seasonal cycle it sampled, and
the column is drawn only when some variable has more than `_MIN_IAV_YEARS` of
them.

## Partitioning leaderboard

`partitioning.jl` produces a separate leaderboard of the dimensionless fractions
describing how the model splits energy and water: the evaporative fraction
`lhf / (lhf + shf)`, the surface runoff fraction `sr / (sr + ssr)`, the
evaporative index `et / precip`, and the transpiration fraction `trans / et`.
Its layout matches the annual leaderboard, with the same `SIM`, `ANN`, `LAT`,
`MON` and `IAV` columns.

Each row is a `PartitionFraction`, so adding one means appending an entry to
`PARTITION_FRACTIONS` rather than touching the plotting code:

```julia
PartitionFraction(
    "ef",                # short name used in file names
    "Evaporative Fraction",
    "LHF/(LHF+SHF)",     # label shown on rows and axes
    ["lhf"],             # diagnostics summed for the numerator
    ["lhf", "shf"],      # diagnostics summed for the denominator
    String[],            # components taken from the simulation on the obs side
    (-0.2, 0.2),         # bias panel color scale
    nothing,             # published constraint, when there are no gridded obs
)
```

Two things to keep in mind:

**Aggregate, then divide.** Every panel sums the numerator and the denominator
over its domain or period and divides only at the end. The mean of a ratio is
not the ratio of the means, and a per-cell `et / precip` diverges wherever
precipitation approaches zero, so the order matters both physically and
numerically. Cells whose aggregated denominator falls below
`_MIN_DENOMINATOR_FRACTION` of its global mean are masked out for the same
reason.

**Observations are optional.** A row is compared against ERA5 only when every
component that is not listed as prescribed is available from the loader, which
is built from the `era5_to_clima_names` keyword argument
(`ERA5_PARTITION_TO_CLIMA_NAMES` by default). Components that *are* listed as prescribed
(such as precipitation, which is forcing) are taken from the simulation on both
sides. A row with no gridded observations, like the transpiration fraction,
shows the simulation alone next to the published global constraint given in its
`reference` field.
