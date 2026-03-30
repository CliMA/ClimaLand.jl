import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaComms
import ClimaLand
import ClimaCalibrate
import Dates
import ClimaCore
import EnsembleKalmanProcesses as EKP
import JLD2

# Access CalibrateConfig
include(
    joinpath(pkgdir(ClimaLand), "experiments/calibration/run_calibration.jl"),
)

include("observation_utils.jl")

# For now, we will reuse `data_sources.jl` that is used for making the
# leaderboard, since it is the easiest option.
include(
    joinpath(pkgdir(ClimaLand), "ext/land_sim_vis/leaderboard/data_sources.jl"),
)

"""
    make_era5_observation_vector(
        covar_estimator,
        short_names,
        sample_date_ranges,
        nelements,
    )

Return a vector of `EKP.Observation` consisting of variables from ERA5 dataset
with `short_names`.

The covariance matrix for each observation is determined by `covar_estimator`
and the date ranges for each observation is determined by `sample_date_ranges`.

The ocean mask is determined by `nelements`.

!!! note "Add new variable"
    To add a new variable from the ERA5 dataset, you must add the variable to
    `get_era5_obs_var_dict` in `data_sources.jl`. In addition, if the varaible
    requires any further or different preprocessing than the preprocessing done
    in `preprocess_single_era5_var`, you should also make any necessary changes
    to that function and in `process_member_data` in `observation_map.jl`.
"""
function make_era5_observation_vector(
    covar_estimator,
    short_names,
    sample_date_ranges,
    nelements,
)
    # TODO: The start_date argument can be removed by making data_sources.jl not
    # dependent on a start date.

    # The start date doesn't matter since we never resample along the
    # time dimension, so we grab the first date in sample_date_ranges
    start_date = first(first(sample_date_ranges))
    era5_vars = preprocess_era5_vars(short_names, start_date, nelements)
    observation_vector = map(sample_date_ranges) do (start_date, stop_date)
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            era5_vars,
            start_date,
            stop_date,
        )
    end
    return observation_vector
end

# TODO: Add a note on how to add a new variable and what you need to do
# Modify dict, maybe modify how they should be processed

"""
    preprocess_era5_vars(short_names, start_date, nelements)

Preprocess each variable from the ERA5 dataset with `short_names`.
"""
function preprocess_era5_vars(short_names, start_date, nelements)
    obs_vars = get_calibration_obs_var_dict()
    for short_name in short_names
        short_name ∉ keys(obs_vars) && error(
            "There is no variable with the short name $short_name. Add this variable to get_calibration_obs_var_dict",
        )
    end
    vars = map(short_names) do short_name
        var = obs_vars[short_name](start_date)
        preprocess_single_era5_var(var, short_name, nelements)
    end
    return vars
end

"""
    preprocess_single_era5_var(var::OutputVar, short_name, nelements)

Specifies how each individual `OutputVar` from the ERA5 dataset should be
processed.

The currently supported variables are `lhf`, `shf`, `lwu`, `swu`, and `gpp`. The
preprocessing is
- computing seasonal averages,
- resampling to fit the model grid,
- applying an ocean mask,
- removing the last longitude point to avoid double counting,
- and excluding the poles.

Note that an `EKP.Observation` is not generated in this function.
"""
function preprocess_single_era5_var(var::OutputVar, short_name, nelements)
    lats, lons = get_lat_lon_from_resolution(nelements)

    # Resampling is not `NaN` aware, so replace `NaN`s with zeros before
    # resampling. This is safe because the ocean mask is applied afterward,
    # removing ocean points where NaNs typically occur (e.g., FLUXCOM GPP).
    if any(isnan, var.data)
        @warn "Replacing NaNs with zeros in $short_name before resampling"
        var = ClimaAnalysis.replace(var, NaN => 0.0)
    end

    # Window to ensure that each season contains all three months.
    # The ERA5 bounds (1979-03 to 2024-08) are only applied when the data's
    # time range covers them; non-ERA5 sources (e.g., FLUXCOM GPP) may have a
    # narrower range and are used as-is.
    start_date = Dates.DateTime(var.attributes["start_date"])
    time_vals = ClimaAnalysis.times(var)
    t_first = start_date + Dates.Second(round(Int, time_vals[begin]))
    t_last = start_date + Dates.Second(round(Int, time_vals[end]))
    window_left = max(Dates.DateTime(1979, 3), t_first)
    window_right = min(Dates.DateTime(2024, 8), t_last)
    var = ClimaAnalysis.window(
        var,
        "time",
        left = window_left,
        right = window_right,
        by = ClimaAnalysis.MatchValue(),
    )

    # Take seasonal average, resample, and apply mask Resampling is an expensive
    # operation, so it is good to do as many reductions as we can.
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)

    var = ClimaAnalysis.resampled_as(var, lon = lons, lat = lats)

    # Cannot apply ClimaLand.apply_oceanmask because of the small
    # differences between the ClimaLand mask and ClimaAnalysis.apply_ocean_mask
    # For now, it is better to manually create a mask from ClimaCore and
    # generate a mask from it using ClimaAnalysis
    ocean_mask = make_ocean_mask(nelements)
    var = ocean_mask(var)

    # To prevent double counting along the longitudes since -180 and 180 degrees
    # are the same point
    var = ClimaAnalysis.window(
        var,
        "longitude",
        right = length(lons) - 1,
        by = ClimaAnalysis.Index(),
    )

    # Exclude the poles
    var = ClimaAnalysis.window(
        var,
        "latitude",
        left = 2,
        right = length(lats) - 1,
        by = ClimaAnalysis.Index(),
    )

    var = ClimaCalibrate.ObservationRecipe.change_data_type(var, Float32)
    return var
end

if abspath(PROGRAM_FILE) == @__FILE__
    (; obs_vec_filepath) = CALIBRATE_CONFIG
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 3.0, # for GPP, 3 g C m-2 day-1
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )
    nelements = CALIBRATE_CONFIG.nelements
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    @info "The number of samples is $(length(sample_date_ranges))"

    isfile(obs_vec_filepath) &&
        @warn "Overwriting the file $obs_vec_filepath to generate the vector of observations"

    observation_vector = make_era5_observation_vector(
        covar_estimator,
        short_names,
        sample_date_ranges,
        nelements,
    )
    JLD2.save_object(obs_vec_filepath, observation_vector)
end
