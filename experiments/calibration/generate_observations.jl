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
    make_observation_vector(
        covar_estimator,
        short_names,
        sample_date_ranges,
        nelements,
    )

Return a vector of `EKP.Observation` consisting of observational variables
with `short_names`.

The covariance matrix for each observation is determined by `covar_estimator`
and the date ranges for each observation is determined by `sample_date_ranges`.

The ocean mask is determined by `nelements`.

Supports both ERA5 variables (lhf, shf, lwu, swu) and ILAMB variables (gpp, et)
via `get_calibration_obs_var_dict` in `data_sources.jl`.
"""
function make_observation_vector(
    covar_estimator,
    short_names,
    sample_date_ranges,
    nelements,
)
    # The start date doesn't matter since we never resample along the
    # time dimension, so we grab the first date in sample_date_ranges
    start_date = first(first(sample_date_ranges))
    obs_vars = preprocess_obs_vars(short_names, start_date, nelements)
    observation_vector = map(sample_date_ranges) do (start_date, stop_date)
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            obs_vars,
            start_date,
            stop_date,
        )
    end
    return observation_vector
end

"""
    preprocess_obs_vars(short_names, start_date, nelements)

Preprocess each observational variable with `short_names`.
"""
function preprocess_obs_vars(short_names, start_date, nelements)
    obs_var_dict = get_calibration_obs_var_dict()
    for short_name in short_names
        short_name ∉ keys(obs_var_dict) && error(
            "There is no variable with the short name $short_name. Add this variable to get_calibration_obs_var_dict in data_sources.jl",
        )
    end
    vars = map(short_names) do short_name
        var = obs_var_dict[short_name](start_date)
        preprocess_single_obs_var(var, short_name, nelements)
    end
    return vars
end

"""
    preprocess_single_obs_var(var::OutputVar, short_name, nelements)

Specifies how each individual `OutputVar` should be processed for calibration.

The preprocessing is:
- replacing NaNs with zeros (required for resampling),
- windowing to full seasons within the data's time range,
- computing seasonal averages,
- resampling to fit the model grid,
- applying an ocean mask,
- removing the last longitude point to avoid double counting,
- and excluding the poles.
"""
function preprocess_single_obs_var(var::OutputVar, short_name, nelements)
    lats, lons = get_lat_lon_from_resolution(nelements)

    # Replace NaNs with zeros so that resampling can be performed.
    # This is safe because the ocean mask will remove ocean points later.
    if any(isnan, var.data)
        @info "Replacing NaNs with zeros in $short_name for resampling"
        var = ClimaAnalysis.replace(var, NaN => 0.0)
    end

    # Window to ensure that each season contains all three months.
    # Use the data's own time range, clamped to full seasons.
    times = ClimaAnalysis.times(var)
    t_min = Dates.DateTime(Dates.year(first(times)), 3)
    t_max = Dates.DateTime(Dates.year(last(times)), 8)
    var = ClimaAnalysis.window(
        var,
        "time",
        left = t_min,
        right = t_max,
        by = ClimaAnalysis.MatchValue(),
    )

    # Take seasonal average, resample, and apply mask. Resampling is an
    # expensive operation, so it is good to do as many reductions as we can.
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)

    var = ClimaAnalysis.resampled_as(var, lon = lons, lat = lats)

    # Cannot apply ClimaLand.apply_oceanmask because of the small
    # differences between the ClimaLand mask and ClimaAnalysis.apply_ocean_mask
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
        scalar = 3.0,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )
    nelements = CALIBRATE_CONFIG.nelements
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    @info "The number of samples is $(length(sample_date_ranges))"

    isfile(obs_vec_filepath) &&
        @warn "Overwriting the file $obs_vec_filepath to generate the vector of observations"

    observation_vector = make_observation_vector(
        covar_estimator,
        short_names,
        sample_date_ranges,
        nelements,
    )
    JLD2.save_object(obs_vec_filepath, observation_vector)
end
