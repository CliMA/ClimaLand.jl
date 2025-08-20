import ClimaAnalysis
import ClimaAnalysis: OutputVar
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
    era5_obs_vars = get_era5_obs_var_dict()
    for short_name in short_names
        short_name ∉ keys(era5_obs_vars) && error(
            "There is no variable with the short name $short_name. Add this variable to get_era5_obs_var_dict",
        )
    end
    vars = map(short_names) do short_name
        var = era5_obs_vars[short_name](start_date)
        preprocess_single_era5_var(var, short_name, nelements)
    end
    return vars
end

"""
    preprocess_single_era5_var(var::OutputVar, short_name, nelements)

Specifies how each individual `OutputVar` from the ERA5 dataset should be
processed.

The currently supported variables are `lhf`, `shf`, `lwu`, and `swu`. The
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

    # If there are `NaN`s, then resampling cannot be performed as resampling is
    # not `NaN` aware
    any(isnan, var.data) && error(
        "Cannot process OutputVar with name $short_name because `NaN`s are present in the data",
    )

    # Window to ensure that each season contains all three months
    # The dates are found by inspecting the ERA5 data and choosing
    # the earliest and latest dates that contain full seasons
    var = ClimaAnalysis.window(
        var,
        "time",
        left = Dates.DateTime(1979, 3),
        right = Dates.DateTime(2024, 8),
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
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 25.0,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )
    nelements = CALIBRATE_CONFIG.nelements
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    @info "The number of samples is $(length(sample_date_ranges))"

    observation_vector = make_era5_observation_vector(
        covar_estimator,
        short_names,
        sample_date_ranges,
        nelements,
    )
    JLD2.save_object(
        "experiments/calibration/land_observation_vector.jld2",
        observation_vector,
    )
end
