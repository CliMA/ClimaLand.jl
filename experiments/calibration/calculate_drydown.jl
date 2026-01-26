# Calculate spatially-explicit soil moisture variability from SMAP L3_SM_P_E data

# The SMAP-derived coefficient of variation (CV) serves as a proxy for landscape-scale soil moisture dynamics. 
# While SMAP data are coarser (9 km) than ERA5 forcing (0.25°), the temporal CV at each SMAP pixel captures the characteristic 
# variability of soil moisture within that landscape patch. High CV indicates highly dynamic soil moisture 
# (e.g., sandy soils, rain-fed agriculture), while low CV indicates buffered conditions (e.g., deep soils, irrigated areas). 
# This temporal variability metric complements the spatial CWD gradient, providing information about soil moisture dynamics 
# that is not captured by climate forcing alone. The CV is particularly valuable for calibrating plant hydraulic responses, as it represents 
# the range of soil water conditions plants experience at a given location.

using NCDatasets
using HDF5
using Dates
using Statistics
using ClimaCore
import ClimaComms
using Interpolations
using JLD2
import ClimaLand

include("smap_utils.jl")

"""
    calculate_sm_variability_from_smap(smap_files::Vector{String}, start_date::DateTime, stop_date::DateTime)

Calculate soil moisture variability metrics from SMAP timeseries.

Computes both standard deviation and coefficient of variation (CV = σ/μ) of soil moisture
over the time period, which captures both absolute and relative variability.

Returns:
- lons: longitude grid (degrees East, EASE-Grid 2.0)
- lats: latitude grid (degrees North, EASE-Grid 2.0)
- sm_std_spatial: standard deviation of soil moisture (cm³/cm³) at each pixel
- sm_cv_spatial: coefficient of variation (dimensionless) at each pixel
- sm_mean_spatial: mean soil moisture (cm³/cm³) for context
"""
function calculate_sm_variability_from_smap(
    smap_files::Vector{String},
    start_date::DateTime,
    stop_date::DateTime;
    chunk_size::Int = 500  # Process 500 lat lines at a time
)
    @info "Processing $(length(smap_files)) SMAP files"
    
    # Read first file to get dimensions and variable names
    # SMAP L3_SM_P_E has hierarchical structure with groups
    ds = NCDataset(smap_files[1])
    
    # SMAP L3_SM_P_E structure: 
    # Group: "Soil_Moisture_Retrieval_Data_AM" or "Soil_Moisture_Retrieval_Data_PM"
    # Variables in group: longitude, latitude, soil_moisture
    
    # Try to find the correct group and variables
    group_name = if haskey(ds.group, "Soil_Moisture_Retrieval_Data_AM")
        "Soil_Moisture_Retrieval_Data_AM"
    elseif haskey(ds.group, "Soil_Moisture_Retrieval_Data_PM")
        "Soil_Moisture_Retrieval_Data_PM"
    else
        @info "Available groups:" keys(ds.group)
        error("Cannot find SMAP soil moisture group")
    end
    
    @info "Using SMAP group: $group_name"
    
    # Access the group
    sm_group = ds.group[group_name]
    
    # Get coordinates from the group
    lons = Array(sm_group["longitude"])
    lats = Array(sm_group["latitude"])
    nlat, nlon = size(lons)
    
    close(ds)
    
    @info "Grid dimensions: $nlat x $nlon"
    
    # Initialize output arrays
    sm_std_spatial = fill(NaN, nlat, nlon)
    sm_cv_spatial = fill(NaN, nlat, nlon)
    sm_mean_spatial = fill(NaN, nlat, nlon)
    
    # Process in latitude chunks to reduce memory
    n_chunks = ceil(Int, nlat / chunk_size)
    
    for chunk_idx in 1:n_chunks
        lat_start = (chunk_idx - 1) * chunk_size + 1
        lat_end = min(chunk_idx * chunk_size, nlat)
        
        @info "Processing chunk $chunk_idx/$n_chunks (lats $lat_start:$lat_end)"
        
        # Accumulate time series for this chunk only
        # Initialize as Float64 array (will use NaN for missing values)
        sm_timeseries = fill(NaN, lat_end - lat_start + 1, nlon, length(smap_files))
        
        for (t_idx, smap_file) in enumerate(smap_files)
            ds = NCDataset(smap_file)
            sm_group = ds.group[group_name]
            
            # Read and convert Missing to NaN in one step
            sm_data = Array(sm_group["soil_moisture"][lat_start:lat_end, :])
            sm_timeseries[:, :, t_idx] .= coalesce.(sm_data, NaN)
            
            close(ds)
            
            # Free memory every 100 files
            if t_idx % 100 == 0
                @info "  Processed $t_idx/$(length(smap_files)) files"
                GC.gc()
            end
        end
        
        # Calculate statistics for this chunk (no need to convert, already Float64)
        @info "  Calculating statistics for chunk $chunk_idx..."
        Threads.@threads for i in 1:(lat_end - lat_start + 1)
            for j in 1:nlon
                ts = sm_timeseries[i, j, :]
                valid = filter(!isnan, ts)
                
                if length(valid) >= 30  # Require at least 30 valid observations
                    sm_mean_spatial[lat_start + i - 1, j] = mean(valid)
                    sm_std_spatial[lat_start + i - 1, j] = std(valid)
                    
                    if sm_mean_spatial[lat_start + i - 1, j] > 0
                        sm_cv_spatial[lat_start + i - 1, j] = 
                            sm_std_spatial[lat_start + i - 1, j] / sm_mean_spatial[lat_start + i - 1, j]
                    end
                end
            end
        end
        
        # Force garbage collection between chunks
        sm_timeseries = nothing
        GC.gc()
    end
    
    # Print statistics
    valid_cv = filter(!isnan, sm_cv_spatial)
    @info "Calculated CV statistics:" n_valid=length(valid_cv) mean_cv=mean(valid_cv) median_cv=median(valid_cv)
    
    return lons, lats, sm_std_spatial, sm_cv_spatial, sm_mean_spatial
end

"""
    regrid_sm_variability_to_model_space(lons, lats, sm_var_spatial, surface_space, FT)

Regrid the SM variability field from lat/lon grid to model space.
"""
function regrid_sm_variability_to_model_space(lons, lats, sm_var_spatial, surface_space, FT)
    # This function should interpolate from the SMAP grid to the model grid
    # Similar to regrid_cwd_to_model_space but for SM variability
    
    # For now, a simple nearest-neighbor approach:
    # Get model coordinates
    coords = ClimaCore.Fields.coordinate_field(surface_space)
    
    # Create a field to store regridded values
    sm_var_field = ClimaCore.Fields.Field(FT(0), surface_space)
    
    # For each model grid point, find nearest SMAP pixel
    for i in eachindex(coords)
        model_lat = coords.lat[i]
        model_lon = coords.lon[i]
        
        # Find nearest grid point in SMAP data
        lat_idx = argmin(abs.(lats .- model_lat))
        lon_idx = argmin(abs.(lons .- model_lon))
        
        sm_var_field[i] = FT(sm_var_spatial[lon_idx, lat_idx])
    end
    
    return sm_var_field
end

# Main execution
function main(start_year::Int = 2015, stop_year::Int = 2015; metric::Symbol = :cv)
    context = ClimaComms.context()
    
    start_date = DateTime(start_year, 4, 1)  # SMAP starts April 2015
    stop_date = DateTime(stop_year, 12, 31)
    
    # Get data directory
    data_dir = get(ENV, "SMAP_DATA_PATH", "/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E")
    
    # Find SMAP files
    smap_files = find_smap_files(start_date, stop_date; data_dir=data_dir)
    
    # If no files found, attempt download
    if isempty(smap_files)
        @warn "No SMAP files found. Attempting download..."
        
        try
            # Download for the date range
            smap_files = download_smap_daterange(
                start_date, 
                stop_date, 
                data_dir
            )
        catch e
            @error "Failed to download SMAP data" exception=(e, catch_backtrace())
            error("""
            No SMAP files found and download failed.
            
            Please either:
            1. Download manually from https://nsidc.org/data/spl3smp_e/versions/6
            2. Check your Earthdata credentials
            3. Verify SMAP_DATA_PATH is set correctly: $data_dir
            """)
        end
    end
    
    if isempty(smap_files)
        error("No SMAP files available after download attempt")
    end
    
    @info "Found $(length(smap_files)) SMAP files for $start_year-$stop_year"
    
    # Calculate SM variability
    @info "Processing $(length(smap_files)) SMAP files for soil moisture variability"
    lons, lats, sm_var_spatial = calculate_sm_variability_from_smap(
        smap_files, 
        start_date, 
        stop_date;
        metric = metric
    )
    
    # Select which metric to use as primary output
    sm_var_spatial = metric == :cv ? sm_cv_spatial : sm_std_spatial
    metric_name = metric == :cv ? "CV" : "Std Dev"
    
    println("=" ^ 60)
    println("Soil moisture variability ($metric_name) calculated on SMAP EASE-Grid for $start_year-$stop_year")
    println("Shape: ", size(sm_var_spatial))
    println("=" ^ 60)
    
    # Save all metrics
    if start_year == stop_year
        output_file = joinpath(@__DIR__, "sm_variability_smap_$(start_year).jld2")
    else
        output_file = joinpath(@__DIR__, "sm_variability_smap_$(start_year)_$(stop_year).jld2")
    end
    
    jldsave(
        output_file; 
        lons, 
        lats, 
        sm_std_spatial, 
        sm_cv_spatial, 
        sm_mean_spatial,
        sm_var_spatial,  # Primary metric
        metric,
        start_year, 
        stop_year
    )
    @info "Saved soil moisture variability to" output_file
    
    return lons, lats, sm_var_spatial
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Default: use coefficient of variation for normalized variability
    lons, lats, sm_var_spatial = main(metric = :cv)
end