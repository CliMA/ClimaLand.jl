# SMAP L3_SM_P_E soil moisture data utilities
# Data: SPL3SMP_E v6 - Enhanced L3 Radiometer Global Daily 9 km EASE-Grid Soil Moisture

using Dates
using NCDatasets
using Statistics
using Glob
using Interpolations  # Add to your using statements at the top

"""
    find_smap_files(start_date::DateTime, stop_date::DateTime; data_dir::String)

Find SMAP L3_SM_P_E files for given date range.
Looks in data_dir (default from environment variable SMAP_DATA_PATH).

# Arguments
- `start_date::DateTime`: Start date for data search
- `stop_date::DateTime`: End date for data search
- `data_dir::String`: Directory containing SMAP HDF5 files (default from ENV or "/net/sampo/data1/smap/L3_SM_P_E")

# Returns
- `Vector{String}`: Sorted list of file paths
"""
function find_smap_files(
    start_date::DateTime, 
    stop_date::DateTime;
    data_dir::String = get(ENV, "SMAP_DATA_PATH", "/net/sampo/data1/smap/L3_SM_P_E")
)
    if !isdir(data_dir)
        @warn "SMAP data directory not found: $data_dir"
        return String[]
    end
    
    files = String[]
    for date in start_date:Day(1):stop_date
        date_str = Dates.format(date, "yyyymmdd")
        # Pattern matches: SMAP_L3_SM_P_E_YYYYMMDD_Rxxxxx_001.h5
        pattern = "SMAP_L3_SM_P_E_$(date_str)_R*.h5"
        matching = glob(pattern, data_dir)
        append!(files, matching)
    end
    
    @info "Found $(length(files)) SMAP files for $(Dates.format(start_date, "yyyy-mm-dd")) to $(Dates.format(stop_date, "yyyy-mm-dd"))"
    
    return sort(files)
end

"""
    load_smap_sm(filepath::String, lat::Float64, lon::Float64; variable::String="soil_moisture")

Load soil moisture data from a SMAP HDF5 file for a specific location.

# Arguments
- `filepath::String`: Path to SMAP HDF5 file
- `lat::Float64`: Latitude in degrees (-90 to 90)
- `lon::Float64`: Longitude in degrees (-180 to 180)
- `variable::String`: Variable name to read (default: "soil_moisture")

# Returns
- `Union{Float64, Nothing}`: Soil moisture value [m³/m³] or nothing if unavailable
"""
function load_smap_sm(filepath::String, lat::Float64, lon::Float64; variable::String="soil_moisture")
    try
        NCDataset(filepath, "r") do ds
            # SMAP L3_SM_P_E uses EASE-Grid 2.0 projection (9 km resolution)
            # The data is stored in groups: /Soil_Moisture_Retrieval_Data_AM or /Soil_Moisture_Retrieval_Data_PM
            
            # Try AM retrieval first, then PM
            for period in ["AM", "PM"]
                group_name = "Soil_Moisture_Retrieval_Data_$(period)"
                
                if haskey(ds.group, group_name)
                    group = ds.group[group_name]
                    
                    if haskey(group, variable)
                        # Get the soil moisture data
                        sm_data = group[variable][:]
                        
                        # Get latitude and longitude grids
                        if haskey(group, "latitude") && haskey(group, "longitude")
                            lats = group["latitude"][:]
                            lons = group["longitude"][:]
                            
                            # Find nearest grid cell
                            # Calculate distances
                            distances = @. sqrt((lats - lat)^2 + (lons - lon)^2)
                            min_idx = argmin(distances)
                            
                            # Get soil moisture value
                            sm_value = sm_data[min_idx]
                            
                            # Check for fill value / missing data
                            # SMAP uses -9999.0 as fill value
                            if sm_value < -1000 || isnan(sm_value)
                                continue  # Try next period
                            end
                            
                            # Check quality flag if available
                            if haskey(group, "retrieval_qual_flag")
                                qual_flag = group["retrieval_qual_flag"][min_idx]
                                # Quality flag: 0 = good, use only good retrievals
                                if qual_flag != 0
                                    continue
                                end
                            end
                            
                            return Float64(sm_value)
                        end
                    end
                end
            end
        end
        
        return nothing
        
    catch e
        @warn "Error reading SMAP file: $filepath" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    calculate_sm_statistics(files::Vector{String}, lat::Float64, lon::Float64)

Calculate soil moisture statistics from a time series of SMAP files.

# Arguments
- `files::Vector{String}`: List of SMAP file paths
- `lat::Float64`: Latitude in degrees
- `lon::Float64`: Longitude in degrees

# Returns
- `Dict{String, Float64}`: Dictionary with keys:
  - "mean": Mean soil moisture [m³/m³]
  - "std": Standard deviation [m³/m³]
  - "min": Minimum value [m³/m³]
  - "max": Maximum value [m³/m³]
  - "median": Median value [m³/m³]
  - "p25": 25th percentile [m³/m³]
  - "p75": 75th percentile [m³/m³]
  - "count": Number of valid observations
"""
function calculate_sm_statistics(files::Vector{String}, lat::Float64, lon::Float64)
    sm_values = Float64[]
    
    for file in files
        sm = load_smap_sm(file, lat, lon)
        if !isnothing(sm)
            push!(sm_values, sm)
        end
    end
    
    if isempty(sm_values)
        @warn "No valid soil moisture data found for location (lat=$lat, lon=$lon)"
        return Dict{String, Float64}(
            "mean" => NaN,
            "std" => NaN,
            "min" => NaN,
            "max" => NaN,
            "median" => NaN,
            "p25" => NaN,
            "p75" => NaN,
            "count" => 0.0
        )
    end
    
    return Dict{String, Float64}(
        "mean" => mean(sm_values),
        "std" => std(sm_values),
        "min" => minimum(sm_values),
        "max" => maximum(sm_values),
        "median" => median(sm_values),
        "p25" => quantile(sm_values, 0.25),
        "p75" => quantile(sm_values, 0.75),
        "count" => Float64(length(sm_values))
    )
end

"""
    load_smap_timeseries(files::Vector{String}, lat::Float64, lon::Float64)

Load a time series of soil moisture values from SMAP files.

# Arguments
- `files::Vector{String}`: List of SMAP file paths (should be sorted by date)
- `lat::Float64`: Latitude in degrees
- `lon::Float64`: Longitude in degrees

# Returns
- `Tuple{Vector{DateTime}, Vector{Float64}}`: Tuple of (dates, soil_moisture_values)
"""
function load_smap_timeseries(files::Vector{String}, lat::Float64, lon::Float64)
    dates = DateTime[]
    sm_values = Float64[]
    
    for file in files
        # Extract date from filename: SMAP_L3_SM_P_E_YYYYMMDD_Rxxxxx_001.h5
        m = match(r"SMAP_L3_SM_P_E_(\d{8})_", basename(file))
        if isnothing(m)
            @warn "Could not parse date from filename: $(basename(file))"
            continue
        end
        
        date = DateTime(m.captures[1], "yyyymmdd")
        sm = load_smap_sm(file, lat, lon)
        
        if !isnothing(sm)
            push!(dates, date)
            push!(sm_values, sm)
        end
    end
    
    return dates, sm_values
end

"""
    load_smap_grid(filepath::String; quality_flag_threshold=0, period_preference=["AM", "PM"])

Load full SMAP soil moisture grid for a given date.
Returns (lons, lats, soil_moisture) on EASE-Grid 2.0 (9 km).

# Arguments
- `filepath::String`: Path to SMAP HDF5 file
- `quality_flag_threshold::Int`: Maximum quality flag to accept (default: 0 = good only)
- `period_preference::Vector{String}`: Order of preference for AM/PM retrievals

# Returns
- `Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}`: (lons, lats, soil_moisture)
  Each is a 2D array matching the EASE-Grid 2.0 projection (typically 1624 x 3856 for 9km)
  Missing/bad data is filled with NaN.
"""
function load_smap_grid(
    filepath::String; 
    quality_flag_threshold::Int = 0,
    period_preference::Vector{String} = ["AM", "PM"]
)
    try
        NCDataset(filepath, "r") do ds
            for period in period_preference
                group_name = "Soil_Moisture_Retrieval_Data_$(period)"
                
                if haskey(ds.group, group_name)
                    group = ds.group[group_name]
                    
                    # Load coordinate grids
                    if haskey(group, "latitude") && haskey(group, "longitude")
                        lats_raw = group["latitude"][:]
                        lons_raw = group["longitude"][:]
                        # Replace Missing with NaN for coordinates
                        lats = replace(lats_raw, missing => NaN) |> x -> convert(Vector{Float64}, x)
                        lons = replace(lons_raw, missing => NaN) |> x -> convert(Vector{Float64}, x)
                        
                        # Load soil moisture
                        if haskey(group, "soil_moisture")
                            sm_raw = group["soil_moisture"][:]
                            # Replace Missing with NaN before converting to Float64
                            sm = replace(sm_raw, missing => NaN) |> x -> convert(Array{Float64}, x)
                            
                            # Apply quality filtering
                            if haskey(group, "retrieval_qual_flag")
                                qual = group["retrieval_qual_flag"][:]
                                # Set poor quality retrievals to NaN
                                sm[qual .> quality_flag_threshold] .= NaN
                            end
                            
                            # Filter fill values (SMAP uses -9999.0)
                            sm[sm .< -1000] .= NaN
                            
                            @info "Loaded SMAP grid" period size(sm) valid_pixels=count(.!isnan.(sm))
                            
                            return (lons, lats, sm)
                        end
                    end
                end
            end
            
            error("No valid SMAP data found in file: $filepath")
        end
        
    catch e
        @error "Error reading SMAP grid from file: $filepath" exception=(e, catch_backtrace())
        rethrow(e)
    end
end

"""
    load_smap_grid_timeseries(files::Vector{String}; quality_flag_threshold=0)

Load a time series of SMAP soil moisture grids from multiple files.
Returns spatial grids and dates.

# Arguments
- `files::Vector{String}`: Sorted list of SMAP file paths
- `quality_flag_threshold::Int`: Maximum quality flag to accept

# Returns
- `Tuple`: (lons, lats, sm_timeseries, dates)
  - `lons::Array{Float64,2}`: Longitude grid
  - `lats::Array{Float64,2}`: Latitude grid  
  - `sm_timeseries::Array{Float64,3}`: Soil moisture [lon, lat, time]
  - `dates::Vector{DateTime}`: Date for each time slice
"""
function load_smap_grid_timeseries(
    files::Vector{String}; 
    quality_flag_threshold::Int = 0
)
    if isempty(files)
        error("No SMAP files provided")
    end
    
    # Load first file to get grid dimensions
    lons_ref, lats_ref, sm_first = load_smap_grid(files[1]; quality_flag_threshold)
    
    # SMAP data is 1D flattened array - Pre-allocate 2D array [n_pixels, time]
    n_pixels = length(sm_first)
    nt = length(files)
    sm_timeseries = Array{Float64,2}(undef, n_pixels, nt)
    dates = Vector{DateTime}(undef, nt)
    
    # Load first file data
    sm_timeseries[:, 1] = sm_first
    
    # Extract date from first filename
    m = match(r"SMAP_L3_SM_P_E_(\d{8})_", basename(files[1]))
    dates[1] = DateTime(m.captures[1], "yyyymmdd")
    
    # Load remaining files
    for (i, file) in enumerate(files[2:end])
        idx = i + 1
        
        # Extract date
        m = match(r"SMAP_L3_SM_P_E_(\d{8})_", basename(file))
        if isnothing(m)
            @warn "Could not parse date from filename: $(basename(file)), skipping"
            sm_timeseries[:, idx] .= NaN
            dates[idx] = DateTime(0)
            continue
        end
        dates[idx] = DateTime(m.captures[1], "yyyymmdd")
        
        # Load soil moisture
        try
            _, _, sm = load_smap_grid(file; quality_flag_threshold)
            
            # Verify grid consistency (1D flattened array)
            if length(sm) != n_pixels
                @warn "Grid size mismatch in file: $(basename(file)), expected $n_pixels pixels, got $(length(sm))"
                sm_timeseries[:, idx] .= NaN
            else
                sm_timeseries[:, idx] = sm
            end
        catch e
            @warn "Failed to load SMAP grid from: $(basename(file))" exception=(e, catch_backtrace())
            sm_timeseries[:, idx] .= NaN
        end
    end
    
    @info "Loaded SMAP timeseries" n_files=nt n_pixels date_range=(minimum(dates), maximum(dates))
    
    return (lons_ref, lats_ref, sm_timeseries, dates)
end

"""
    calculate_smap_temporal_mean(files::Vector{String}; quality_flag_threshold=0)

Calculate temporal mean soil moisture from multiple SMAP files using streaming approach.
This avoids loading the entire timeseries into memory at once.

# Arguments
- `files::Vector{String}`: List of SMAP file paths
- `quality_flag_threshold::Int`: Maximum quality flag to accept

# Returns
- `Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}`: (lons, lats, sm_mean)
  - `sm_mean`: Temporal mean (1D flattened), with NaN where data was missing for all times
"""
function calculate_smap_temporal_mean(
    files::Vector{String}; 
    quality_flag_threshold::Int = 0
)
    if isempty(files)
        error("No SMAP files provided")
    end
    
    # Load first file to get grid dimensions
    lons_ref, lats_ref, sm_first = load_smap_grid(files[1]; quality_flag_threshold)
    n_pixels = length(sm_first)
    
    # Initialize accumulators for online mean calculation (Welford's algorithm)
    sum_sm = zeros(Float64, n_pixels)
    count_valid = zeros(Int, n_pixels)
    
    @info "Streaming SMAP temporal mean calculation" n_files=length(files) n_pixels
    
    # Process each file incrementally
    for (idx, file) in enumerate(files)
        try
            lons, lats, sm = load_smap_grid(file; quality_flag_threshold)
            
            # Verify grid consistency
            if length(sm) != n_pixels
                @warn "Grid size mismatch, skipping file" file expected=n_pixels got=length(sm)
                continue
            end
            
            # Accumulate sum and count for non-NaN values
            for i in 1:n_pixels
                if !isnan(sm[i])
                    sum_sm[i] += sm[i]
                    count_valid[i] += 1
                end
            end
            
            # Progress update every 100 files
            if idx % 100 == 0
                avg_coverage = mean(count_valid) / idx
                @info "Progress: $idx/$(length(files)) files" avg_valid_pixels_per_file=round(Int, avg_coverage * n_pixels)
            end
            
        catch e
            @warn "Error loading SMAP file, skipping" file exception=e
            continue
        end
    end
    
    # Calculate mean: sum / count, NaN where count == 0
    sm_mean = [count_valid[i] > 0 ? sum_sm[i] / count_valid[i] : NaN for i in 1:n_pixels]
    
    n_valid_pixels = count(count_valid .> 0)
    mean_obs_per_pixel = mean(count_valid[count_valid .> 0])
    
    @info "Calculated SMAP temporal mean" n_files=length(files) n_valid_pixels mean_obs_per_pixel=round(mean_obs_per_pixel, digits=1) coverage_pct=round(100*n_valid_pixels/n_pixels, digits=2)
    
    return (lons_ref, lats_ref, sm_mean)
end

"""
    regrid_smap_to_coords(
        lons_smap::Array{Float64,2}, 
        lats_smap::Array{Float64,2}, 
        sm_smap::Array{Float64,2},
        target_lons::Vector{Float64},
        target_lats::Vector{Float64}
    )

Regrid SMAP data to target coordinates using nearest-neighbor or linear interpolation.

# Arguments
- `lons_smap::Array{Float64,2}`: SMAP longitude grid (2D)
- `lats_smap::Array{Float64,2}`: SMAP latitude grid (2D)
- `sm_smap::Array{Float64,2}`: SMAP soil moisture grid (2D)
- `target_lons::Vector{Float64}`: Target longitudes (1D vector)
- `target_lats::Vector{Float64}`: Target latitudes (1D vector)

# Returns
- `Vector{Float64}`: Regridded soil moisture values at target coordinates
"""
function regrid_smap_to_coords(
    lons_smap::Array{Float64,2}, 
    lats_smap::Array{Float64,2}, 
    sm_smap::Array{Float64,2},
    target_lons::Vector{Float64},
    target_lats::Vector{Float64}
)
    n_targets = length(target_lons)
    @assert length(target_lats) == n_targets "Longitude and latitude vectors must have same length"
    
    sm_regrid = Vector{Float64}(undef, n_targets)
    
    # Flatten SMAP grids for nearest-neighbor search
    lons_flat = vec(lons_smap)
    lats_flat = vec(lats_smap)
    sm_flat = vec(sm_smap)
    
    # Build KDTree for efficient nearest-neighbor search (if you have NearestNeighbors.jl)
    # For now, use simple distance search
    for i in 1:n_targets
        target_lon = target_lons[i]
        target_lat = target_lats[i]
        
        # Handle longitude wrapping
        target_lon_wrapped = mod(target_lon + 180, 360) - 180
        
        # Find nearest SMAP pixel
        distances = @. sqrt((lons_flat - target_lon_wrapped)^2 + (lats_flat - target_lat)^2)
        min_idx = argmin(distances)
        
        sm_regrid[i] = sm_flat[min_idx]
    end
    
    @info "Regridded SMAP to target coordinates" n_targets valid_targets=count(.!isnan.(sm_regrid))
    
    return sm_regrid
end

"""
    get_smap_spatial_subset(
        lons_smap::Array{Float64,2}, 
        lats_smap::Array{Float64,2}, 
        sm_smap::Array{Float64,2};
        lon_range::Tuple{Float64,Float64} = (-180.0, 180.0),
        lat_range::Tuple{Float64,Float64} = (-90.0, 90.0)
    )

Extract a spatial subset of SMAP data within specified lon/lat bounds.

# Arguments
- `lons_smap::Array{Float64,2}`: SMAP longitude grid
- `lats_smap::Array{Float64,2}`: SMAP latitude grid
- `sm_smap::Array{Float64,2}`: SMAP soil moisture grid
- `lon_range::Tuple{Float64,Float64}`: (min_lon, max_lon)
- `lat_range::Tuple{Float64,Float64}`: (min_lat, max_lat)

# Returns
- `Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}`: Subset (lons, lats, sm)
"""
function get_smap_spatial_subset(
    lons_smap::Array{Float64,2}, 
    lats_smap::Array{Float64,2}, 
    sm_smap::Array{Float64,2};
    lon_range::Tuple{Float64,Float64} = (-180.0, 180.0),
    lat_range::Tuple{Float64,Float64} = (-90.0, 90.0)
)
    # Find indices within bounds
    mask = (lons_smap .>= lon_range[1]) .& (lons_smap .<= lon_range[2]) .&
           (lats_smap .>= lat_range[1]) .& (lats_smap .<= lat_range[2])
    
    # Get bounding box indices
    rows, cols = findall(mask) |> x -> (getindex.(x, 1), getindex.(x, 2))
    row_min, row_max = extrema(rows)
    col_min, col_max = extrema(cols)
    
    # Extract subset
    lons_subset = lons_smap[row_min:row_max, col_min:col_max]
    lats_subset = lats_smap[row_min:row_max, col_min:col_max]
    sm_subset = sm_smap[row_min:row_max, col_min:col_max]
    
    @info "Extracted SMAP spatial subset" original_size=size(sm_smap) subset_size=size(sm_subset) lon_range lat_range
    
    return (lons_subset, lats_subset, sm_subset)
end

"""
    create_valid_observation_mask(lons_smap, lats_smap, sm_smap, land_mask)

Create a combined mask of:
1. Land pixels (from model's land mask)
2. Valid SMAP observations (not NaN)

Returns indices and coordinates of valid observations only.
"""
function create_valid_observation_mask(
    lons_smap::Array{Float64,2},
    lats_smap::Array{Float64,2},
    sm_smap::Array{Float64,2},
    land_mask::Vector{Bool}  # From your model domain
)
    # Flatten SMAP grids
    lons_flat = vec(lons_smap)
    lats_flat = vec(lats_smap)
    sm_flat = vec(sm_smap)
    
    # Create combined mask: land AND valid SMAP data
    valid_smap = .!isnan.(sm_flat)
    valid_mask = land_mask .& valid_smap
    
    # Get indices and values of valid observations
    valid_indices = findall(valid_mask)
    
    obs_vector = sm_flat[valid_mask]
    obs_lons = lons_flat[valid_mask]
    obs_lats = lats_flat[valid_mask]
    
    @info "Valid observations" total_land=sum(land_mask) valid_smap=sum(valid_smap) valid_both=sum(valid_mask) coverage_pct=round(100*sum(valid_mask)/sum(land_mask), digits=2)
    
    return (
        observation_vector = obs_vector,
        valid_indices = valid_indices,
        lons = obs_lons,
        lats = obs_lats,
        mask = valid_mask
    )
end

"""
    create_valid_observation_mask(
        lons_smap, lats_smap, sm_smap, model_coords
    )

Map model grid points to SMAP pixels and create observation vector.

Returns SMAP pixel assignments so model output can be spatially averaged 
to match SMAP's 9km resolution during forward model evaluation.

# Returns
- observation_vector: SMAP observations (one per unique SMAP pixel with valid data)
- smap_pixel_groups: Vector of vectors - each inner vector contains model grid indices 
                     that fall within one SMAP pixel
- smap_pixel_values: SMAP observation value for each pixel group
- coverage_stats: Coverage statistics
"""
# Overload for 1D vectors (from streaming temporal mean)
function create_valid_observation_mask(
    lons_smap::Vector{Float64},
    lats_smap::Vector{Float64},
    sm_smap::Vector{Float64},
    model_coords::NamedTuple
)
    n_model = length(model_coords.lons)
    
    @info "Assigning model points to SMAP pixels (1D)..." n_model n_smap=length(sm_smap)
    
    # For each model point, find nearest SMAP pixel
    smap_assignments = fill(0, n_model)  # Which SMAP pixel does each model point belong to?
    
    # Assign each model point to nearest SMAP pixel (with distance threshold)
    # SMAP resolution is ~9km, so use conservative threshold of ~0.15° (~16.7 km at equator)
    # This ensures we only match model points that actually overlap with SMAP pixels
    max_distance_deg = 0.15  # degrees
    max_dist_sq = max_distance_deg^2
    
    for i in 1:n_model
        # Skip ocean points
        !model_coords.land_mask[i] && continue
        
        lon_model = model_coords.lons[i]
        lat_model = model_coords.lats[i]
        
        # Find nearest SMAP pixel (Euclidean distance in lat-lon space)
        dist_sq = (lons_smap .- lon_model).^2 .+ (lats_smap .- lat_model).^2
        min_dist_sq = minimum(dist_sq)
        
        # Only assign if within threshold distance
        if min_dist_sq <= max_dist_sq
            nearest_idx = argmin(dist_sq)
            smap_assignments[i] = nearest_idx
        end
        # Otherwise leave as 0 (no assignment)
    end
    
    # Group model indices by SMAP pixel
    smap_pixel_groups = Dict{Int, Vector{Int}}()
    
    for i in 1:n_model
        smap_idx = smap_assignments[i]
        smap_idx == 0 && continue  # Skip ocean
        
        if !haskey(smap_pixel_groups, smap_idx)
            smap_pixel_groups[smap_idx] = Int[]
        end
        push!(smap_pixel_groups[smap_idx], i)
    end
    
    @info "Model points grouped into SMAP pixels" n_unique_smap_pixels=length(smap_pixel_groups)
    
    # Filter to only SMAP pixels with valid data
    valid_smap_pixels = Int[]
    valid_model_groups = Vector{Vector{Int}}()
    valid_pixel_coords = Vector{NamedTuple{(:lons, :lats), Tuple{Vector{Float64}, Vector{Float64}}}}()
    valid_smap_values = Float64[]
    
    for (smap_idx, model_indices) in smap_pixel_groups
        sm_val = sm_smap[smap_idx]
        
        # Only include if SMAP has valid data (not NaN)
        if !isnan(sm_val)
            push!(valid_smap_pixels, smap_idx)
            push!(valid_model_groups, model_indices)
            push!(valid_smap_values, sm_val)
            
            # Store lat/lon coordinates for this SMAP pixel's model points
            pixel_lons = model_coords.lons[model_indices]
            pixel_lats = model_coords.lats[model_indices]
            push!(valid_pixel_coords, (lons = pixel_lons, lats = pixel_lats))
        end
    end
    
    @info "Valid SMAP pixels identified" n_valid=length(valid_smap_pixels)
    
    # Calculate coverage statistics
    n_land = sum(model_coords.land_mask)
    n_model_in_valid_smap = sum(length(group) for group in valid_model_groups)
    
    coverage_stats = Dict(
        "total_model_points" => n_model,
        "land_model_points" => n_land,
        "unique_smap_pixels" => length(smap_pixel_groups),
        "valid_smap_pixels" => length(valid_smap_pixels),
        "model_points_in_valid_smap" => n_model_in_valid_smap,
        "land_coverage_pct" => round(100 * n_model_in_valid_smap / max(n_land, 1), digits=2),
        "avg_model_points_per_smap_pixel" => round(n_model_in_valid_smap / max(length(valid_smap_pixels), 1), digits=2)
    )
    
    @info "Coverage statistics" total_model_points=n_model land_model_points=n_land unique_smap_pixels=length(smap_pixel_groups) valid_smap_pixels=length(valid_smap_pixels) model_points_in_valid_smap=n_model_in_valid_smap land_coverage_pct=coverage_stats["land_coverage_pct"] avg_model_points_per_smap_pixel=coverage_stats["avg_model_points_per_smap_pixel"]
    
    return (
        observation_vector = valid_smap_values,  # SMAP observations (one per pixel)
        smap_pixel_groups = valid_model_groups,  # Model indices to average per SMAP pixel (legacy, for debugging)
        smap_pixel_coords = valid_pixel_coords,  # Model lat/lon coords for each SMAP pixel
        smap_pixel_indices = valid_smap_pixels,  # SMAP linear indices
        coverage_stats = coverage_stats
    )
end

# Original 2D matrix version
function create_valid_observation_mask(
    lons_smap::Array{Float64,2},
    lats_smap::Array{Float64,2},
    sm_smap::Array{Float64,2},
    model_coords::NamedTuple
)
    n_model = length(model_coords.lons)
    
    # For each model point, find nearest SMAP pixel
    smap_assignments = fill(0, n_model)  # Which SMAP pixel does each model point belong to?
    
    # Flatten SMAP grid for easier indexing
    smap_flat_idx = LinearIndices(sm_smap)
    lons_smap_flat = vec(lons_smap)
    lats_smap_flat = vec(lats_smap)
    sm_smap_flat = vec(sm_smap)
    
    @info "Assigning model points to SMAP pixels..." n_model
    
    # Assign each model point to nearest SMAP pixel (with distance threshold)
    # SMAP resolution is ~9km, so use conservative threshold of ~0.15° (~16.7 km at equator)
    # This ensures we only match model points that actually overlap with SMAP pixels
    max_distance_deg = 0.15  # degrees
    max_dist_sq = max_distance_deg^2
    
    for i in 1:n_model
        # Skip ocean points
        !model_coords.land_mask[i] && continue
        
        lon_model = model_coords.lons[i]
        lat_model = model_coords.lats[i]
        
        # Find nearest SMAP pixel (Euclidean distance in lat-lon space)
        # Note: For better accuracy, could use great circle distance
        dist_sq = (lons_smap_flat .- lon_model).^2 .+ (lats_smap_flat .- lat_model).^2
        min_dist_sq = minimum(dist_sq)
        
        # Only assign if within threshold distance
        if min_dist_sq <= max_dist_sq
            nearest_idx = argmin(dist_sq)
            smap_assignments[i] = nearest_idx
        end
        # Otherwise leave as 0 (no assignment)
    end
    
    # Group model indices by SMAP pixel
    smap_pixel_groups = Dict{Int, Vector{Int}}()
    
    for i in 1:n_model
        smap_idx = smap_assignments[i]
        smap_idx == 0 && continue  # Skip ocean
        
        if !haskey(smap_pixel_groups, smap_idx)
            smap_pixel_groups[smap_idx] = Int[]
        end
        push!(smap_pixel_groups[smap_idx], i)
    end
    
    @info "Model points grouped into SMAP pixels" n_unique_smap_pixels=length(smap_pixel_groups)
    
    # Filter to only SMAP pixels with valid data
    valid_smap_pixels = Int[]
    valid_model_groups = Vector{Vector{Int}}()
    valid_pixel_coords = Vector{NamedTuple{(:lons, :lats), Tuple{Vector{Float64}, Vector{Float64}}}}()
    valid_smap_values = Float64[]
    
    for (smap_idx, model_indices) in smap_pixel_groups
        sm_val = sm_smap_flat[smap_idx]
        
        # Only include if SMAP has valid data (not NaN)
        if !isnan(sm_val)
            push!(valid_smap_pixels, smap_idx)
            push!(valid_model_groups, model_indices)
            push!(valid_smap_values, sm_val)
            
            # Store lat/lon coordinates for this SMAP pixel's model points
            pixel_lons = model_coords.lons[model_indices]
            pixel_lats = model_coords.lats[model_indices]
            push!(valid_pixel_coords, (lons = pixel_lons, lats = pixel_lats))
        end
    end
    
    @info "Valid SMAP pixels identified" n_valid=length(valid_smap_pixels)
    
    # Calculate coverage statistics
    n_land = sum(model_coords.land_mask)
    n_model_in_valid_smap = sum(length(group) for group in valid_model_groups)
    
    coverage_stats = Dict(
        "total_model_points" => n_model,
        "land_model_points" => n_land,
        "unique_smap_pixels" => length(smap_pixel_groups),
        "valid_smap_pixels" => length(valid_smap_pixels),
        "model_points_in_valid_smap" => n_model_in_valid_smap,
        "land_coverage_pct" => round(100 * n_model_in_valid_smap / max(n_land, 1), digits=2),
        "avg_model_points_per_smap_pixel" => round(n_model_in_valid_smap / max(length(valid_smap_pixels), 1), digits=2)
    )
    
    @info "Coverage statistics" total_model_points=n_model land_model_points=n_land unique_smap_pixels=length(smap_pixel_groups) valid_smap_pixels=length(valid_smap_pixels) model_points_in_valid_smap=n_model_in_valid_smap land_coverage_pct=coverage_stats["land_coverage_pct"] avg_model_points_per_smap_pixel=coverage_stats["avg_model_points_per_smap_pixel"]
    
    return (
        observation_vector = valid_smap_values,  # SMAP observations (one per pixel)
        smap_pixel_groups = valid_model_groups,  # Model indices to average per SMAP pixel (legacy, for debugging)
        smap_pixel_coords = valid_pixel_coords,  # Model lat/lon coords for each SMAP pixel
        smap_pixel_indices = valid_smap_pixels,  # SMAP linear indices
        coverage_stats = coverage_stats
    )
end

"""
    extract_model_coords_from_space(surface_space, land_mask_field)

Extract coordinate vectors and land mask from ClimaCore surface space.
"""
function extract_model_coords_from_space(surface_space, land_mask_field)
    coords = ClimaCore.Fields.coordinate_field(surface_space)
    
    # Extract and flatten coordinates
    lons_field = coords.long
    lats_field = coords.lat
    
    lons_vec = vec(parent(lons_field))
    lats_vec = vec(parent(lats_field))
    land_mask_vec = vec(parent(land_mask_field)) .> 0.5  # Convert to boolean
    
    @info "Extracted model coordinates" n_points=length(lons_vec) n_land=sum(land_mask_vec)
    
    return (
        lons = lons_vec,
        lats = lats_vec,
        land_mask = land_mask_vec
    )
end