"""
Generate SMAP observation vector from pre-regridded lat-lon SMAP data.

This simplified version loads SMAP that's already been regridded to the ClimaAnalysis
output grid (404×202 lat-lon), eliminating the need for coordinate matching.
"""

using Dates
using JLD2
using LinearAlgebra: I
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP

include(joinpath(@__DIR__, "api.jl"))
include(joinpath(@__DIR__, "run_uspac_calibration.jl"))

"""
    load_regridded_smap(regridded_filepath)

Load pre-regridded SMAP data from lat-lon file.
Returns observation values and valid indices (non-NaN points).
"""
function load_regridded_smap(regridded_filepath::String)
    if !isfile(regridded_filepath)
        error("Regridded SMAP file not found: $regridded_filepath. Run regrid_smap_to_model_grid.jl first.")
    end
    
    @info "Loading regridded SMAP data" filepath=regridded_filepath
    
    data = JLD2.load(regridded_filepath)
    sm_latlon = data["sm_latlon"]  # 2D array (n_lon, n_lat)
    
    # Flatten to 1D
    sm_flat = vec(sm_latlon)
    
    # Find valid (non-NaN) points
    valid_mask = .!isnan.(sm_flat)
    valid_indices = findall(valid_mask)
    valid_values = sm_flat[valid_mask]
    
    @info "Loaded regridded SMAP" grid_size=size(sm_latlon) total_points=length(sm_flat) valid_points=length(valid_values) coverage_pct=round(100*length(valid_values)/length(sm_flat), digits=2)
    
    return (
        observation_values = valid_values,
        valid_indices = valid_indices,
        grid_size = size(sm_latlon),
        total_points = length(sm_flat),
        metadata = data
    )
end

"""
    make_smap_observation_vector_from_regridded(regridded_filepath, covar_estimator)

Create EKP.Observation from pre-regridded SMAP data.
Much simpler than the original version - no spatial aggregation needed!
"""
function make_smap_observation_vector_from_regridded(
    regridded_filepath::String,
    covar_estimator::ClimaCalibrate.ObservationRecipe.ScalarCovariance
)
    # Load regridded SMAP
    smap_data = load_regridded_smap(regridded_filepath)
    
    obs_values = Float32.(smap_data.observation_values)
    n_obs = length(obs_values)
    
    @info "Creating observation vector" n_observations=n_obs
    
    # Create covariance matrix
    cov_matrix = covar_estimator.scalar * I(n_obs)
    
    # Create EKP.Observation
    obs = EKP.Observation(
        Dict(
            "samples" => [obs_values],
            "covariances" => [cov_matrix],
            "names" => ["sm_surface"]
        )
    )
    
    # Return as vector (for consistency with multi-period format)
    return [obs], smap_data
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # Input: regridded SMAP file
    regridded_smap_file = joinpath(@__DIR__, "smap_regridded_latlon.jld2")
    
    # Output: observation vector file  
    obs_output_file = joinpath(@__DIR__, "land_observation_vector.jld2")
    
    # Covariance estimator (same as before)
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 0.01,  # (m³/m³)²
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )
    
    @info "Generating SMAP observation vector from regridded data" input=regridded_smap_file output=obs_output_file
    
    if !isfile(regridded_smap_file)
        error("Regridded SMAP file not found: $regridded_smap_file\nPlease run: julia regrid_smap_to_model_grid.jl")
    end
    
    # Create observation vector
    observation_vector, smap_data = make_smap_observation_vector_from_regridded(
        regridded_smap_file,
        covar_estimator
    )
    
    # Save observation vector with metadata
    isfile(obs_output_file) && @warn "Overwriting $obs_output_file"
    
    JLD2.jldsave(
        obs_output_file;
        observation_vector,
        valid_indices = smap_data.valid_indices,  # Which grid points have observations
        grid_size = smap_data.grid_size,  # (404, 202)
        total_points = smap_data.total_points,  # 81,608
        sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges,
        short_names = ["sm_surface"],
        created = now(),
        description = "SMAP observations from pre-regridded 404×202 lat-lon data (ClimaAnalysis output grid)"
    )
    
    @info "SMAP observation vector created" filepath=obs_output_file n_observations=length(smap_data.observation_values)
    @info "Next steps:" 
    @info "  1. observation_map.jl will use valid_indices to extract matching simulation data"
    @info "  2. Run calibration with: qsub PBS_calibration.pbs"
end
