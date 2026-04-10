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
        noise_scalars,
        short_names,
        sample_date_ranges,
        nelements,
    )

Return a vector of `EKP.Observation` consisting of observational variables
with `short_names`.

Each variable gets its own covariance scaling from `noise_scalars`, a
`Dict{String, Float64}` mapping short names to scalar values. Observations
are created per-variable and combined via `EKP.combine_observations`.

The date ranges for each observation are determined by `sample_date_ranges`.

The ocean mask is determined by `nelements`.

Supports both ERA5 variables (lhf, shf, lwu, swu) and ILAMB variables (gpp, et)
via `get_calibration_obs_var_dict` in `data_sources.jl`.
"""
function make_observation_vector(
    noise_scalars,
    short_names,
    sample_date_ranges,
    nelements,
)
    # The start date doesn't matter since we never resample along the
    # time dimension, so we grab the first date in sample_date_ranges
    start_date = first(first(sample_date_ranges))
    obs_vars = preprocess_obs_vars(short_names, start_date, nelements)

    # Build per-variable covariance estimators
    covar_estimators = Dict(
        name => ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
            scalar = noise_scalars[name],
            use_latitude_weights = true,
            min_cosd_lat = 0.1,
        ) for name in short_names
    )

    observation_vector = map(sample_date_ranges) do (start_date, stop_date)
        # Create a separate observation for each variable with its own scalar
        per_var_obs = map(zip(short_names, obs_vars)) do (name, var)
            ClimaCalibrate.ObservationRecipe.observation(
                covar_estimators[name],
                [var],
                start_date,
                stop_date,
            )
        end
        # Combine into a single observation with block-diagonal covariance
        EKP.combine_observations(per_var_obs)
    end
    return observation_vector
end

"""
    preprocess_obs_vars(short_names, start_date, nelements)

Preprocess each observational variable with `short_names`.
"""
function preprocess_obs_vars(short_names, start_date, nelements)
    obs_var_dict = get_calibration_obs_var_dict(; short_names)
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
- windowing to full seasons within the data's time range,
- computing seasonal averages,
- resampling to fit the model grid,
- applying an ocean mask,
- removing the last longitude point to avoid double counting,
- and excluding the poles.
"""
function preprocess_single_obs_var(var::OutputVar, short_name, nelements)
    lats, lons = get_lat_lon_from_resolution(nelements)

    # NaNs are kept so that resampling propagates them rather than
    # interpolating good observations with zeros. Some valid points near
    # NaN regions may be lost, but this is preferred over corrupting them.

    # Window to ensure that each season contains all three months.
    # Use the data's own time range, clamped to full seasons.
    times = ClimaAnalysis.times(var)
    first_time = first(times)
    last_time = last(times)
    t_min = Dates.DateTime(Dates.year(first_time), 3)
    t_max = Dates.DateTime(Dates.year(last_time), 8)
    # Ensure bounds are within the data range
    if t_max > last_time
        t_max = Dates.DateTime(Dates.year(last_time) - 1, 8)
    end
    if t_min < first_time
        t_min = Dates.DateTime(Dates.year(first_time) + 1, 3)
    end
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
    (; obs_vec_filepath, nelements, sample_date_ranges, short_names) =
        CALIBRATE_CONFIG
    @info "The number of samples is $(length(sample_date_ranges))"
    @info "Noise scalars: $NOISE_SCALARS"

    isfile(obs_vec_filepath) &&
        @warn "Overwriting the file $obs_vec_filepath to generate the vector of observations"

    observation_vector = make_observation_vector(
        NOISE_SCALARS,
        short_names,
        sample_date_ranges,
        nelements,
    )
    JLD2.save_object(obs_vec_filepath, observation_vector)
end
