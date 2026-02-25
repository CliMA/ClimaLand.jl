"""
Regrid SMAP L3 temporal mean soil moisture from EASE-Grid to ClimaAnalysis's lat-lon output grid.

This creates a static file with SMAP observations on the same 404×202 lat-lon grid
that ClimaAnalysis uses for diagnostic output, eliminating the need for nearest-neighbor
matching during calibration.
"""

using Dates
using JLD2
using Statistics

include(joinpath(@__DIR__, "smap_utils.jl"))
include(joinpath(@__DIR__, "observation_utils.jl"))
include(joinpath(@__DIR__, "api.jl"))
include(joinpath(@__DIR__, "run_uspac_calibration.jl"))

"""
    regrid_smap_to_latlon(start_date, stop_date; quality_flag_threshold, use_quality_weights, weighting_scheme, n_lon, n_lat)

Load SMAP data, calculate temporal mean, and regrid to regular lat-lon grid.
Matches the grid that ClimaAnalysis uses for diagnostic output (404×202).

# Arguments
- `start_date::DateTime`: Start date for SMAP data
- `stop_date::DateTime`: End date for SMAP data
- `quality_flag_threshold::Int`: Maximum quality flag to accept (0-3)
- `use_quality_weights::Bool`: If true, use quality-weighted averaging; if false, use binary threshold
- `weighting_scheme::Symbol`: Weighting scheme (:inverse, :exponential, :quadratic, :binary)
- `n_lon::Int`: Number of longitude points (default: 404)
- `n_lat::Int`: Number of latitude points (default: 202)

# Returns
- If `use_quality_weights=false`: (sm_latlon, lons_grid, lats_grid)
- If `use_quality_weights=true`: (sm_latlon, lons_grid, lats_grid, weights_latlon, sm_std_latlon, n_obs_latlon)
"""
function regrid_smap_to_latlon(
    start_date::DateTime,
    stop_date::DateTime;
    quality_flag_threshold::Int = 0,
    use_quality_weights::Bool = false,
    weighting_scheme::Symbol = :inverse,
    n_lon::Int = 404,
    n_lat::Int = 202
)
    # Create regular lat-lon grid matching ClimaAnalysis output
    lons_grid = range(-180.0, 180.0, length=n_lon)
    lats_grid = range(-90.0, 90.0, length=n_lat)
    
    @info "Created lat-lon grid" n_lon n_lat total_points=n_lon*n_lat
    
    # Load SMAP data and calculate temporal mean
    smap_files = find_smap_files(start_date, stop_date)
    
    # Choose between quality-weighted or threshold-only averaging
    if use_quality_weights
        @info "Using quality-weighted temporal averaging" threshold=quality_flag_threshold scheme=weighting_scheme
        result = calculate_smap_weighted_temporal_mean(
            smap_files;
            quality_flag_threshold,
            weighting_scheme
        )
        lons_smap = result.lons
        lats_smap = result.lats
        sm_smap = result.sm_mean
        weights_smap = result.total_weight
        sm_std_smap = result.sm_std
        n_obs_smap = result.n_obs
    else
        @info "Using threshold-only temporal averaging" threshold=quality_flag_threshold
        lons_smap, lats_smap, sm_smap = calculate_smap_temporal_mean(
            smap_files;
            quality_flag_threshold
        )
        weights_smap = nothing
        sm_std_smap = nothing
        n_obs_smap = nothing
    end
    
    # Flatten SMAP arrays (they may be multi-dimensional from EASE grid)
    lons_smap_vec = vec(lons_smap)
    lats_smap_vec = vec(lats_smap)
    sm_smap_vec = vec(sm_smap)
    weights_smap_vec = use_quality_weights ? vec(weights_smap) : nothing
    sm_std_smap_vec = use_quality_weights ? vec(sm_std_smap) : nothing
    n_obs_smap_vec = use_quality_weights ? vec(n_obs_smap) : nothing
    
    @info "SMAP temporal mean calculated" n_smap_pixels=length(sm_smap_vec) valid_pixels=count(!isnan, sm_smap_vec)
    
    # Create 2D arrays for target grid
    sm_latlon = fill(NaN, n_lon, n_lat)
    weights_latlon = use_quality_weights ? fill(NaN, n_lon, n_lat) : nothing
    sm_std_latlon = use_quality_weights ? fill(NaN, n_lon, n_lat) : nothing
    n_obs_latlon = use_quality_weights ? zeros(Int, n_lon, n_lat) : nothing
    
    # Regrid SMAP to lat-lon grid
    for j in 1:n_lat
        for i in 1:n_lon
            lon_target = lons_grid[i]
            lat_target = lats_grid[j]
            
            # Find nearest SMAP pixel within threshold
            dist_sq = (lons_smap_vec .- lon_target).^2 .+ (lats_smap_vec .- lat_target).^2
            min_dist_sq = minimum(dist_sq)
            
            # Distance threshold: 0.15° (~16.7 km)
            if min_dist_sq <= 0.15^2
                nearest_idx = argmin(dist_sq)
                sm_val = sm_smap_vec[nearest_idx]
                if !isnan(sm_val)
                    sm_latlon[i, j] = sm_val
                    
                    # Also regrid quality metrics if using weighted approach
                    if use_quality_weights
                        weights_latlon[i, j] = weights_smap_vec[nearest_idx]
                        sm_std_latlon[i, j] = sm_std_smap_vec[nearest_idx]
                        n_obs_latlon[i, j] = n_obs_smap_vec[nearest_idx]
                    end
                end
            end
        end
        
        # Progress update
        if j % 50 == 0
            @info "Progress: $j/$n_lat latitudes processed"
        end
    end
    
    n_valid = count(!isnan, sm_latlon)
    @info "SMAP regridded to lat-lon" n_points=length(sm_latlon) n_valid_points=n_valid coverage_pct=round(100*n_valid/length(sm_latlon), digits=2)
    
    if use_quality_weights
        mean_weight = isempty(weights_latlon[.!isnan.(weights_latlon)]) ? NaN : mean(weights_latlon[.!isnan.(weights_latlon)])
        mean_nobs = isempty(n_obs_latlon[n_obs_latlon .> 0]) ? 0 : mean(n_obs_latlon[n_obs_latlon .> 0])
        @info "Quality metrics" mean_total_weight=round(mean_weight, digits=2) mean_n_obs=round(mean_nobs, digits=1)
        return sm_latlon, lons_grid, lats_grid, weights_latlon, sm_std_latlon, n_obs_latlon
    else
        return sm_latlon, lons_grid, lats_grid
    end
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    output_file = joinpath(@__DIR__, "smap_regridded_latlon.jld2")
    
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    
    # Assuming single period calibration
    start_date, stop_date = sample_date_ranges[1]
    
    # Use ClimaAnalysis's standard lat-lon grid dimensions
    n_lon = 404
    n_lat = 202
    
    # Quality control settings
    # Option 1: Simple threshold (backward compatible)
    # quality_flag_threshold = 0
    # use_quality_weights = false
    
    # Option 2: Quality weighting (NEW - recommended for more data coverage)
    quality_flag_threshold = 2  # Include all except poor quality (flag=3)
    use_quality_weights = true
    weighting_scheme = :inverse  # Options: :inverse, :exponential, :quadratic, :binary
    
    @info "Regridding SMAP to lat-lon grid" n_lon n_lat start_date stop_date output_file use_weights=use_quality_weights
    
    result_tuple = regrid_smap_to_latlon(
        start_date,
        stop_date;
        quality_flag_threshold,
        use_quality_weights,
        weighting_scheme,
        n_lon,
        n_lat
    )
    
    # Unpack results (handle both weighted and non-weighted cases)
    if use_quality_weights
        sm_latlon, lons_grid, lats_grid, weights_latlon, sm_std_latlon, n_obs_latlon = result_tuple
    else
        sm_latlon, lons_grid, lats_grid = result_tuple
    end
    
    # Save regridded SMAP field
    @info "Saving regridded SMAP data to $output_file"
    
    if use_quality_weights
        # Save with quality metrics
        JLD2.jldsave(
            output_file;
            sm_latlon,
            weights_latlon,
            sm_std_latlon,
            n_obs_latlon,
            lons = collect(lons_grid),
            lats = collect(lats_grid),
            start_date,
            stop_date,
            quality_flag_threshold,
            use_quality_weights,
            weighting_scheme,
            created = now(),
            description = "SMAP L3 soil moisture regridded to 404×202 lat-lon grid with quality weights (threshold=$quality_flag_threshold, scheme=$weighting_scheme)"
        )
        @info "Regridding complete!" output_file n_valid=count(!isnan, sm_latlon) method="quality-weighted"
    else
        # Save simple version (backward compatible)
        JLD2.jldsave(
            output_file;
            sm_latlon,
            lons = collect(lons_grid),
            lats = collect(lats_grid),
            start_date,
            stop_date,
            quality_flag_threshold,
            use_quality_weights,
            created = now(),
            description = "SMAP L3 soil moisture regridded to 404×202 lat-lon grid (threshold=$quality_flag_threshold, no quality weighting)"
        )
        @info "Regridding complete!" output_file n_valid=count(!isnan, sm_latlon) method="threshold-only"
    end
end
