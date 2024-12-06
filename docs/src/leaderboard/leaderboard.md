# Leaderboard

## Long run

### Add a new variable to compare against observations
The infrastructure to compute errors against observations is in the `leaderboard` folder.
This folder contains two files: `data_sources.jl`, responsible for loading and preprocessing
variables of interest, and `leaderboard.jl`, which computes error and draw plots. To add a
new variable to the comparison, you modify the `data_sources.jl`.

### Computation
As of now, the leaderboard produces bias plots with the global bias and global root mean
squared error (RMSE). These quantities are computed for each month with the first year of
the simulation not considered as that is the spinup time. The start date of the simulation
is 2008 which means that only the year 2009 is used to compare against observational data.

### Add a new variable to the bias plots
There are four functions that you need to modify to add a new variable which are
`get_sim_var_dict`, `get_obs_var_dict`, `get_mask_dict`, and
`get_compare_vars_biases_plot_extrema`. Each function returns a dictionary that must be
modified to add a new variable to the leaderboard. The dictionaries are `sim_var_dict`,
`obs_var_dict`, `mask_dict`, and `compare_vars_biases_plot_extrema`.

To add a variable for the leaderboard, add a key-value pair to the dictionary `sim_var_dict`
whose key is the short name of the variable and the value is a function that returns a
[`OutputVar`](https://clima.github.io/ClimaAnalysis.jl/dev/var/). Any preprocessing is done
in the function which includes unit conversion and shifting the dates.

```julia
sim_var_dict["et"] =
        () -> begin
            # Load in variable
            sim_var = get(
                ClimaAnalysis.SimDir(diagnostics_folder_path),
                short_name = "et",
            )
            # Shift to the first day and subtract one month as preprocessing
            sim_var =
                ClimaAnalysis.shift_to_start_of_previous_month(sim_var)
            return sim_var
        end
```

Then, add a key-value pair to the dictionary `obs_var_dict` whose key is the same short name
as before and the value is a function that takes in a start date and returns a `OutputVar`.
Any preprocessing is done in the function.

```julia
obs_var_dict["et"] =
    (start_date) -> begin
    # We use ClimaArtifacts to use a dataset from ILAMB
        obs_var = ClimaAnalysis.OutputVar(
            ClimaLand.Artifacts.ilamb_dataset_path(;
                context = "evspsbl_MODIS_et_0.5x0.5.nc",
            ),
            "et",
            # start_date is used to align the dates in the observational data
            # with the simulation data
            new_start_date = start_date,
            # Shift dates to the first day of the month before aligning the dates
            shift_by = Dates.firstdayofmonth,
        )
        # More preprocessing to match the units with the simulation data
        ClimaAnalysis.units(obs_var) == "kg/m2/s" &&
            (obs_var = ClimaAnalysis.set_units(obs_var, "kg m^-2 s^-1"))
        # ClimaAnalysis cannot handle `missing` values, but does support handling NaNs
        obs_var = ClimaAnalysis.replace(obs_var, missing => NaN)
        return obs_var
    end
```

!!! tip "Preprocessing"
    Observational and simulational data should be preprocessed for dates and units. When
    using ClimaDiagnostics to report monthly averages from a simulation, monthly averages
    are output on the first day following the month when the average was computed. For
    instance, the monthly average corresponding to January 2010 is on the date 1 Feb 2010.
    Preprocessing is done to shift this date to 1 Jan 2010. When preprocessing data, we
    follow the convention that the first day corresponds to the monthly average for that
    month. For observational data, you should check the convention being followed and
    preprocess the dates if necessary.

    For `obs_var_dict`, the anonymous function must take in a start date. The start date is
    used in `leaderboard.jl` to adjust the seconds in the `OutputVar` to match between start
    date in the simulation data.

    Units should be the same between the simulation and observational data.

Next, add a key-value pair to the dictionary `mask_dict` whose key is the same short name
as before and the value is a function that takes in a `OutputVar` representing simulation
data and a `OutputVar` representing observational data and returns a masking function or
`nothing` if no masking function is needed. The masking function is used to correctly
normalize the global bias and global RMSE. See the example below where a mask is made using
the observational data.

```julia
mask_dict["et"] =
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

Finally, add a key-value pair to the dictionary `compare_vars_biases_plot_extrema` whose
key is the same short name as before and the value is a tuple of floats which determine
the range of the bias plots.

```julia
compare_vars_biases_plot_extrema = Dict(
    "et" => (-0.00001, 0.00001),
    "gpp" => (-8.0, 8.0),
    "lwu" => (-40.0, 40.0),
)
```
