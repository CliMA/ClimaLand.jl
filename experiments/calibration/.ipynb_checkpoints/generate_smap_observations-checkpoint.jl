import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaLand
import ClimaCalibrate
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2

# Grab CALIBRATE_CONFIG from run_calibration.jl (safe: guarded by PROGRAM_FILE)
include(joinpath(pkgdir(ClimaLand), "experiments/calibration/run_calibration.jl"))

# Shared helpers (ocean mask, grid lat/lon, etc.)
include(joinpath(pkgdir(ClimaLand), "experiments/calibration/observation_utils.jl"))

# NEW: SMAP data source registry you’ll create (see section 2)
include(joinpath(pkgdir(ClimaLand), "ext/land_sim_vis/leaderboard/data_sources_smap.jl"))

"""
    make_smap_observation_vector(covar_estimator, short_names, sample_date_ranges, nelements)

Build a vector of EKP.Observation from SMAP soil moisture, resampled to `nelements`.
"""
function make_smap_observation_vector(
    covar_estimator,
    short_names,
    sample_date_ranges,
    nelements,
)
    start_date = first(first(sample_date_ranges))  # not used for resampling; keeps API parity
    smap_vars = preprocess_smap_vars(short_names, start_date, nelements)
    observation_vector = map(sample_date_ranges) do (start_date, stop_date)
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            smap_vars,
            start_date,
            stop_date,
        )
    end
    return observation_vector
end

"""
    preprocess_smap_vars(short_names, start_date, nelements)

Load & preprocess each SMAP variable listed in `short_names`.
"""
function preprocess_smap_vars(short_names, start_date, nelements)
    smap_obs_vars = get_smap_obs_var_dict()  # defined in data_sources_smap.jl
    for short_name in short_names
        short_name ∉ keys(smap_obs_vars) && error(
            "There is no SMAP variable with short name $short_name. Add it to get_smap_obs_var_dict().",
        )
    end
    return map(short_names) do short_name
        var = smap_obs_vars[short_name](start_date)  # returns ClimaAnalysis.OutputVar
        preprocess_single_smap_var(var, short_name, nelements)
    end
end

"""
    preprocess_single_smap_var(var::OutputVar, short_name, nelements)

Preprocess SMAP OutputVar:
- window time to a safe span,
- fill NaNs over land per-time-slice,
- (optionally) aggregate to seasonal average (to match your ERA5-style pipeline),
- resample to model grid,
- ocean mask, trim last longitude, drop poles,
- cast to Float32.
"""
function preprocess_single_smap_var(var::OutputVar, short_name, nelements)
    lats, lons = get_lat_lon_from_resolution(nelements)

    # Window to a broad era covering your date ranges (tune as needed)
    var = ClimaAnalysis.window(
        var, "time",
        left  = Dates.DateTime(2015, 4),  # SMAP start
        right = Dates.DateTime(2025, 1),  # or latest you have
        by = ClimaAnalysis.MatchValue(),
    )

    # Fill NaNs by land-median per time slice (resampling is not NaN-aware)
    function fill_time_slice(v::OutputVar)
        # flatten over space for this time slice
        med = ClimaAnalysis.nanmedian(v.data)
        return ClimaAnalysis.replace(x -> (isnan(x) ? med : x), v)
    end
    # Apply per time index (ClimaAnalysis.apply_along can vary; use loop if needed)
    var = ClimaAnalysis.map_along(var, "time") do v_t
        fill_time_slice(v_t)
    end

    # Choose aggregation to match model diagnostics:
    # If your model uses seasonal means in observation_map, do the same here:
    var = ClimaAnalysis.average_season_across_time(var, ignore_nan = true)

    # Resample to model grid, then mask/trim same as ERA5 path
    var = ClimaAnalysis.resampled_as(var, lon = lons, lat = lats)
    ocean_mask = make_ocean_mask(nelements)
    var = ocean_mask(var)

    var = ClimaAnalysis.window(var, "longitude",
        right = length(lons) - 1, by = ClimaAnalysis.Index())

    var = ClimaAnalysis.window(var, "latitude",
        left = 2, right = length(lats) - 1, by = ClimaAnalysis.Index())

    var = ClimaCalibrate.ObservationRecipe.change_data_type(var, Float32)
    return var
end

# ------------------ Script entry point ------------------
if abspath(PROGRAM_FILE) == @__FILE__
    (; obs_vec_filepath, nelements, sample_date_ranges, short_names) = CALIBRATE_CONFIG

    # Covariance for soil moisture; tune scalar. Keep latitude weighting.
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 0.02,                # e.g., 0.02 m³/m³ (adjust for SMAP uncertainty)
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )

    @info "Generating SMAP observation vector for $(length(sample_date_ranges)) windows"

    isfile(obs_vec_filepath) &&
        @warn "Overwriting $obs_vec_filepath"

    observation_vector = make_smap_observation_vector(
        covar_estimator, short_names, sample_date_ranges, nelements)

    JLD2.save_object(obs_vec_filepath, observation_vector)
    @info "Saved SMAP observation vector to $obs_vec_filepath"
end
