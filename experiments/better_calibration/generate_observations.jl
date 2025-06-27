import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaLand
import ClimaCalibrate
import Dates
import ClimaCore

# Include because we need access to the constants
include("calibrate_land.jl")

"""
    make_observation(short_name, start_date, end_date)

Make an observation from the OutputVar with the name `short_name` from
start_date` to `end_date`.
"""
function make_era5_observation(covar_estimator, short_names, start_date, end_date)
    era5_vars = preprocess_era5_vars(short_names, START_DATE)
    return ClimaCalibrate.ObservationRecipe.observation(
        covar_estimator,
        era5_vars,
        start_date,
        end_date,
    )
end

function preprocess_era5_vars(short_names, args...; kwargs...)
    era5_obs_vars = get_era5_obs_var_dict()
    for short_name in short_names
        short_name ∉ keys(era5_obs_vars) && error(
            "There is not variable with the short name $short_name. Add this variable to get_era5_obs_var_dict",
        )
    end
    vars = map(short_names) do short_name
        var = era5_obs_vars[short_name](args...; kwargs...)
        preprocess_era5_vars(var)
    end
    return vars
end

"""
    preprocess_era5_vars(var::OutputVar, short_name)

Specify how each `OutputVar` from the ERA5 dataset should be processed.
"""
# TODO: It would be nice to just use the short name from `var`, but it need to
# be changed in the data processing
function preprocess_era5_vars(var::OutputVar, short_name)
    if short_name in ("lhf", "shf", "swu", "lwu")

        lats, lons = diagnostics_lat_lon(NELEMENTS)

        # Check for `NaN`s
        # If there are `NaN`s, then resampling can't be done as resampling is
        # not `NaN` aware
        any(isnan, var.data) && error(
            "Cannot process OutputVar with name $short_name because `NaN`s are in the data",
        )

        # Take seasonal average, resample, and apply mask
        # Resampling is an expensive operation, so it is good to do as much
        # reductions that we can.
        var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)

        var = ClimaAnalysis.resampled_as(var, lon = lons, lat = lats)

        # Cannot apply ClimaLand.apply_oceanmask because of the small
        # differences between ClimaLand mask and ClimaAnalysis mask (TODO: Need
        # to verify this) For now, it is better to manually make a mask using
        # ClimaCore and make a mask from that
        ocean_mask = make_ocean_mask(NELEMENTS)
        var = ocean_mask(var)

        # To prevent double count near the poles, we also apply a window operation here
        var = window(var, "longitude", right = length(lons) - 1, by = ClimaAnalysis.Index())

        # Want the covariance matrix to be Float32
        var = ObservationRecipe.change_data_type(var, Float32)
        return var
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    covar_estimator = ClimaCalibrate.ObservationRecipe.SeasonalDiagonalCovariance(
        model_error_scale = Float32(0.05),
        regularization = Float32(25.0),
        ignore_nan = true,
    )
    observation_vector = map(SAMPLE_DATE_RANGES) do start_date, end_date
        make_era5_observation(covar_estimator, SHORT_NAMES, start_date, end_date)
    end

    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, end_date) in SAMPLE_DATE_RANGES
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                MINIBATCH_SIZE,
            ),
        ),
    )
end
