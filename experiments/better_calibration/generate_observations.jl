import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaLand
import ClimaCalibrate
import Dates
import ClimaCore
import EnsembleKalmanProcesses as EKP
import JLD2

include(
    joinpath(pkgdir(ClimaLand), "experiments/better_calibration/config.jl"),
)

include("observation_utils.jl")

# For now, we will use `data_sources.jl` for the leaderboard, since it is the
# easiest option, but it would be better to make your own `data_source.jl` and
# preprocess the observational data to match the simulation data as opposed to
# processing both the simulation and observational data (e.g. with ILAMB data).
include(
    joinpath(pkgdir(ClimaLand), "ext/land_sim_vis/leaderboard/data_sources.jl"),
)

"""
    make_observation(short_name, start_date, end_date, nelements)

Make an observation from the OutputVar with the name `short_name` from
`start_date` to `end_date`.
"""
function make_era5_observation_vector(
    covar_estimator,
    short_names,
    sample_date_ranges,
    nelements,
)
    # TODO: The start_date argument can be removed by making data_sources.jl not
    # depend on a start date, but it is nice to be able to reuse the code for
    # the leaderboard and calibration

    # The start date doesn't matter since we never resample along the
    # time dimension, so we grab the first date in sample_date_ranges
    start_date = first(first(sample_date_ranges))
    era5_vars = preprocess_era5_vars(short_names, start_date, nelements)
    observation_vector = map(sample_date_ranges) do (start_date, end_date)
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            era5_vars,
            start_date,
            end_date,
        )
    end
    return observation_vector
end

"""
    preprocess_era5_vars(short_names, start_date)

Preprocess each ERA5 variable.
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
    preprocess_single_era5_var(var::OutputVar, short_name)

Specify how each individual `OutputVar` from the ERA5 dataset should be
processed.

Note that the individual observation is not generated in this function.
"""
function preprocess_single_era5_var(var::OutputVar, short_name, nelements)
    lats, lons = diagnostics_lat_lon(nelements)

    # If there are `NaN`s, then resampling can't be done as resampling is
    # not `NaN` aware
    any(isnan, var.data) && error(
        "Cannot process OutputVar with name $short_name because `NaN`s are in the data",
    )

    # Window to ensure that each season contain all three months
    # The dates are found by inspecting the ERA5 data and choosing
    # the earliest and latest dates that would contain full seasons
    var = ClimaAnalysis.window(
        var,
        "time",
        left = Dates.DateTime(1979, 3),
        right = Dates.DateTime(2024, 8),
        by = ClimaAnalysis.MatchValue(),
    )

    # Take seasonal average, resample, and apply mask
    # Resampling is an expensive operation, so it is good to do as much
    # reductions that we can.
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)

    var = ClimaAnalysis.resampled_as(var, lon = lons, lat = lats)

    # Cannot apply ClimaLand.apply_oceanmask because of the small
    # differences between ClimaLand mask and ClimaAnalysis.apply_ocean_mask
    # For now, it is better to manually create a mask from ClimaCore and
    # generate a mask from it using ClimaAnalysis

    # TODO: This can be improved by saving the mask as a diagnostic
    ocean_mask = make_ocean_mask(nelements)
    var = ocean_mask(var)

    # To prevent double count near the poles, we also apply another window
    # operation
    var = ClimaAnalysis.window(
        var,
        "longitude",
        right = length(lons) - 1,
        by = ClimaAnalysis.Index(),
    )

    # Throw away the poles for latitudes
    var = ClimaAnalysis.window(
        var,
        "latitude",
        left = 2,
        right = length(lats) - 1,
        by = ClimaAnalysis.Index(),
    )

    # TODO: Remove this for generic calibration
    var = ClimaAnalysis.window(var, "latitude", left = 0.0, right = 45.0)

    var = ClimaAnalysis.window(var, "longitude", left = -60.0, right = 60.0)

    # TODO: Change the comment? not sure if covariance matrix is correct
    # Want the covariance matrix to be Float32
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
        "experiments/better_calibration/land_observation_vector.jld2",
        observation_vector,
    )
end
