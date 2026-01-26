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
    regrid_smap_to_latlon(start_date, stop_date; quality_flag_threshold, n_lon, n_lat)

Load SMAP data, calculate temporal mean, and regrid to regular lat-lon grid.
Matches the grid that ClimaAnalysis uses for diagnostic output (404×202).
"""
function regrid_smap_to_latlon(
    start_date::DateTime,
    stop_date::DateTime;
    quality_flag_threshold::Int = 0,
    n_lon::Int = 404,
    n_lat::Int = 202
)
    # Create regular lat-lon grid matching ClimaAnalysis output
    lons_grid = range(-180.0, 180.0, length=n_lon)
    lats_grid = range(-90.0, 90.0, length=n_lat)
    
    @info "Created lat-lon grid" n_lon n_lat total_points=n_lon*n_lat
    
    # Load SMAP data and calculate temporal mean
    smap_files = find_smap_files(start_date, stop_date)
    lons_smap, lats_smap, sm_smap = calculate_smap_temporal_mean(
        smap_files;
        quality_flag_threshold
    )
    
    # Flatten SMAP arrays (they may be multi-dimensional from EASE grid)
    lons_smap_vec = vec(lons_smap)
    lats_smap_vec = vec(lats_smap)
    sm_smap_vec = vec(sm_smap)
    
    @info "SMAP temporal mean calculated" n_smap_pixels=length(sm_smap_vec) valid_pixels=count(!isnan, sm_smap_vec)
    
    # Create 2D arrays for target grid
    sm_latlon = fill(NaN, n_lon, n_lat)
    
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
    
    return sm_latlon, lons_grid, lats_grid
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
    
    @info "Regridding SMAP to lat-lon grid" n_lon n_lat start_date stop_date output_file
    
    sm_latlon, lons_grid, lats_grid = regrid_smap_to_latlon(
        start_date,
        stop_date;
        quality_flag_threshold = 0,
        n_lon = n_lon,
        n_lat = n_lat
    )
    
    # Save regridded SMAP field
    @info "Saving regridded SMAP data to $output_file"
    
    JLD2.jldsave(
        output_file;
        sm_latlon,
        lons = collect(lons_grid),
        lats = collect(lats_grid),
        start_date,
        stop_date,
        created = now(),
        description = "SMAP L3 soil moisture regridded to 404×202 lat-lon grid (ClimaAnalysis output grid) using nearest-neighbor with 0.15° threshold"
    )
    
    @info "Regridding complete!" output_file n_valid=count(!isnan, sm_latlon)
end
