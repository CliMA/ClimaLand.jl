"""
Generate SMAP observation vector from pre-regridded lat-lon SMAP data.

This simplified version loads SMAP that's already been regridded to the ClimaAnalysis
output grid (404×202 lat-lon), eliminating the need for coordinate matching.
"""

using Dates
using JLD2
using LinearAlgebra: I, Diagonal
using Statistics: mean
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP

include(joinpath(@__DIR__, "api.jl"))
include(joinpath(@__DIR__, "run_uspac_calibration.jl"))

"""
    load_regridded_smap(regridded_filepath)

Load pre-regridded SMAP data from lat-lon file.
Returns observation values, valid indices (non-NaN points), and optional quality weights.
"""
function load_regridded_smap(regridded_filepath::String)
    if !isfile(regridded_filepath)
        error("Regridded SMAP file not found: $regridded_filepath. Run regrid_smap_to_model_grid.jl first.")
    end
    
    @info "Loading regridded SMAP data" filepath=regridded_filepath
    
    data = JLD2.load(regridded_filepath)
    sm_latlon = data["sm_latlon"]  # 2D array (n_lon, n_lat)
    
    # Check if quality weights are available (from new weighted approach)
    has_weights = haskey(data, "weights_latlon") && haskey(data, "use_quality_weights") && data["use_quality_weights"]
    
    if has_weights
        weights_latlon = data["weights_latlon"]
        sm_std_latlon = data["sm_std_latlon"]
        n_obs_latlon = data["n_obs_latlon"]
        @info "Quality-weighted data detected" weighting_scheme=get(data, "weighting_scheme", :unknown) threshold=get(data, "quality_flag_threshold", :unknown)
    end
    
    # Flatten to 1D
    sm_flat = vec(sm_latlon)
    
    # Find valid (non-NaN) points
    valid_mask = .!isnan.(sm_flat)
    valid_indices = findall(valid_mask)
    valid_values = sm_flat[valid_mask]
    
    # Also extract quality metrics for valid points if available
    weights_flat = has_weights ? vec(weights_latlon)[valid_indices] : nothing
    sm_std_flat = has_weights ? vec(sm_std_latlon)[valid_indices] : nothing
    n_obs_flat = has_weights ? vec(n_obs_latlon)[valid_indices] : nothing
    
    @info "Loaded regridded SMAP" grid_size=size(sm_latlon) total_points=length(sm_flat) valid_points=length(valid_values) coverage_pct=round(100*length(valid_values)/length(sm_flat), digits=2) has_quality_weights=has_weights
    
    return (
        observation_values = valid_values,
        valid_indices = valid_indices,
        grid_size = size(sm_latlon),
        total_points = length(sm_flat),
        weights = weights_flat,
        sm_std = sm_std_flat,
        n_obs = n_obs_flat,
        has_weights = has_weights,
        metadata = data
    )
end

"""
    make_smap_observation_vector_from_regridded(regridded_filepath, covar_estimator; use_quality_weights)

Create EKP.Observation from pre-regridded SMAP data.
If quality weights are available and `use_quality_weights=true`, incorporates them into the covariance matrix
so that higher quality observations have lower variance (more influence on calibration).

# Arguments
- `regridded_filepath`: Path to regridded SMAP file
- `covar_estimator`: ScalarCovariance object defining base observation error
- `use_quality_weights::Bool`: If true, use quality weights to modify covariance (default: true)
"""
function make_smap_observation_vector_from_regridded(
    regridded_filepath::String,
    covar_estimator::ClimaCalibrate.ObservationRecipe.ScalarCovariance;
    use_quality_weights::Bool = true
)
    # Load regridded SMAP
    smap_data = load_regridded_smap(regridded_filepath)
    
    obs_values = Float32.(smap_data.observation_values)
    n_obs = length(obs_values)
    
    @info "Creating observation vector" n_observations=n_obs has_quality_weights=smap_data.has_weights
    
    # Create covariance matrix
    if smap_data.has_weights && use_quality_weights
        # Incorporate quality weights into observation error covariance
        # Higher quality (higher weight) -> lower variance -> more influence
        weights = Float64.(smap_data.weights)
        
        # Normalize weights to have mean = 1.0 for interpretability
        # This keeps the base_variance scale meaningful
        normalized_weights = weights ./ mean(weights)
        
        # Variance inversely proportional to quality weight
        # variance_i = base_variance / normalized_weight_i
        observation_variances = covar_estimator.scalar ./ normalized_weights
        
        # Create diagonal covariance matrix with variable variances
        cov_matrix = Diagonal(observation_variances)
        
        # Log statistics
        min_var = minimum(observation_variances)
        max_var = maximum(observation_variances)
        mean_var = mean(observation_variances)
        @info "Quality-weighted covariance created" base_variance=covar_estimator.scalar mean_variance=round(mean_var, digits=6) min_variance=round(min_var, digits=6) max_variance=round(max_var, digits=6) variance_range_ratio=round(max_var/min_var, digits=2)
    else
        # Standard uniform covariance (original behavior)
        cov_matrix = covar_estimator.scalar * I(n_obs)
        @info "Uniform covariance created" variance=covar_estimator.scalar
    end
    
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
    
    # Quality weighting control
    # Set to false to ignore quality weights even if available (for testing/comparison)
    use_quality_weights = true
    
    # Covariance estimator (same as before)
    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = 0.01,  # (m³/m³)² - base observation error variance
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )
    
    @info "Generating SMAP observation vector from regridded data" input=regridded_smap_file output=obs_output_file use_quality_weights
    
    if !isfile(regridded_smap_file)
        error("Regridded SMAP file not found: $regridded_smap_file\nPlease run: julia regrid_smap_to_model_grid.jl")
    end
    
    # Create observation vector
    observation_vector, smap_data = make_smap_observation_vector_from_regridded(
        regridded_smap_file,
        covar_estimator;
        use_quality_weights
    )
    
    # Save observation vector with metadata
    isfile(obs_output_file) && @warn "Overwriting $obs_output_file"
    
    if smap_data.has_weights && use_quality_weights
        # Save with quality weights - they were applied to covariance
        @info "Saving observation vector with quality-weighted covariance"
        JLD2.jldsave(
            obs_output_file;
            observation_vector,
            valid_indices = smap_data.valid_indices,
            weights = smap_data.weights,  # Quality-based weights for each observation
            sm_std = smap_data.sm_std,    # Standard deviation for uncertainty
            n_obs = smap_data.n_obs,      # Number of observations per pixel
            grid_size = smap_data.grid_size,
            total_points = smap_data.total_points,
            sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges,
            short_names = ["sm_surface"],
            has_quality_weights = true,
            quality_weights_applied = true,  # Weights incorporated into covariance
            weighting_scheme = get(smap_data.metadata, "weighting_scheme", :unknown),
            quality_flag_threshold = get(smap_data.metadata, "quality_flag_threshold", 0),
            base_variance = covar_estimator.scalar,
            created = now(),
            description = "SMAP observations with quality-weighted covariance (404×202 lat-lon grid)"
        )
    elseif smap_data.has_weights && !use_quality_weights
        # Weights available but not applied (for comparison)
        @info "Saving observation vector with weights metadata (not applied to covariance)"
        JLD2.jldsave(
            obs_output_file;
            observation_vector,
            valid_indices = smap_data.valid_indices,
            weights = smap_data.weights,
            sm_std = smap_data.sm_std,
            n_obs = smap_data.n_obs,
            grid_size = smap_data.grid_size,
            total_points = smap_data.total_points,
            sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges,
            short_names = ["sm_surface"],
            has_quality_weights = true,
            quality_weights_applied = false,  # Weights NOT incorporated into covariance
            weighting_scheme = get(smap_data.metadata, "weighting_scheme", :unknown),
            quality_flag_threshold = get(smap_data.metadata, "quality_flag_threshold", 0),
            base_variance = covar_estimator.scalar,
            created = now(),
            description = "SMAP observations from quality-weighted data (weights available but not applied to covariance)"
        )
    else
        # Save without weights (backward compatible)
        @info "Saving observation vector with uniform covariance (no quality weights)"
        JLD2.jldsave(
            obs_output_file;
            observation_vector,
            valid_indices = smap_data.valid_indices,
            grid_size = smap_data.grid_size,
            total_points = smap_data.total_points,
            sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges,
            short_names = ["sm_surface"],
            has_quality_weights = false,
            quality_weights_applied = false,
            base_variance = covar_estimator.scalar,
            created = now(),
            description = "SMAP observations with uniform covariance (404×202 lat-lon grid)"
        )
    end
    
    @info "SMAP observation vector created" filepath=obs_output_file n_observations=length(smap_data.observation_values) has_weights=smap_data.has_weights weights_applied=(smap_data.has_weights && use_quality_weights)
    @info "Next steps:" 
    @info "  1. observation_map.jl will use valid_indices to extract matching simulation data"
    if smap_data.has_weights && use_quality_weights
        @info "  2. Quality weights INCORPORATED into observation covariance matrix"
        @info "     - Higher quality observations have lower variance (more influence on calibration)"
        @info "     - Check log above for variance range statistics"
    elseif smap_data.has_weights && !use_quality_weights
        @info "  2. Quality weights available but NOT applied (set use_quality_weights=true to enable)"
    end
    @info "  $(smap_data.has_weights ? "3" : "2"). Run calibration with: qsub PBS_calibration.pbs"
end
