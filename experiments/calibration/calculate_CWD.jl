# Calculate spatially-explicit CWD from ERA5 land forcing artifact
#
# Climatic Water Deficit (CWD) Calculation
# =========================================
# 
# CWD represents the difference between potential evapotranspiration (PET) and 
# actual evapotranspiration (AET), accumulated over a time period. It indicates
# the amount of water that could be evapotranspired if water were not limiting.
# CWD is a key metric for characterizing water stress in ecosystems.
#
# Formula:
# -------
# CWD = PET - AET
#
# Where:
# - PET = Potential evapotranspiration (water demand)
# - AET = Actual evapotranspiration (limited by water availability)
# - AET is approximated as min(PET, Precipitation) in water-limited conditions
#
# Methodology:
# -----------
# 1. For each location and timestep:
#    - Calculate PET from temperature
#    - Estimate AET as the minimum of PET and available water (precipitation)
#    - Accumulate the deficit: CWD += (PET - AET)
#
# 2. Sum CWD over the desired time period (annual or multi-year average)
#
# PET Estimation:
# --------------
# We use a simplified temperature-based formula:
#   PET (mm/day) = 0.05 * T_C * (1 + 0.01 * T_C)  for T_C > 0°C
#   PET = 0 for T_C ≤ 0°C
# 
# This is scaled to hourly timesteps for ERA5 hourly data.
#
# Note: This is a simplified approach suitable for global applications.
# More accurate PET methods (Penman-Monteith) would require additional 
# ERA5 variables (radiation, humidity, wind speed).
#
# References:
# ----------
# - Stephenson, N. L. (1998). "Actual evapotranspiration and deficit: biologically 
#   meaningful correlates of vegetation distribution across spatial scales."
#   Journal of Biogeography, 25(5), 855-870.
#   https://doi.org/10.1046/j.1365-2699.1998.00233.x
#   (Defines CWD = PET - AET as a measure of unmet atmospheric water demand 
#   and demonstrates its utility for predicting vegetation distribution)

using NCDatasets
using Dates
using Statistics
using ClimaCore
using ClimaCore: Fields
import ClimaComms
using Interpolations
using JLD2
import ClimaLand
import CUDA

function calculate_cwd_from_era5(
    era5_paths::Vector{String},
    start_date::DateTime,
    stop_date::DateTime,
)
    @info "Reading ERA5 data from $(length(era5_paths)) file(s)"
    
    # Initialize accumulators
    lons = nothing
    lats = nothing
    CWD_sum = nothing  # Sum of annual CWD across all years
    year_count = 0
    
    # Process each file (year)
    for era5_path in era5_paths
        @info "Processing file: $era5_path"
        ds = NCDataset(era5_path)
        
        # Get coordinates (same for all files)
        if lons === nothing
            lons = Array(ds["lon"])  # degrees East
            lats = Array(ds["lat"])  # degrees North
            nlon, nlat = length(lons), length(lats)
            CWD_sum = zeros(Float64, nlon, nlat)
        end
        
        precip_var = ds["mtpr"]  # m/s - (lon, lat, time)
        temp_var = ds["t2m"]     # K - (lon, lat, time)
        
        nlon, nlat, ntime = size(precip_var)
        
        @info "Processing spatial CWD" nlon nlat ntime
        
        # Timestep in hours (ERA5 is hourly)
        dt_hours = 1.0
        dt_seconds = dt_hours * 3600.0
        
        # Temporary storage for this year's CWD
        CWD_this_year = zeros(Float64, nlon, nlat)
        
        # Process each spatial location
        Threads.@threads for j in 1:nlat
            for i in 1:nlon
                # Extract time series at this location
                precip_ts = precip_var[i, j, :]  # m/s
                temp_ts = temp_var[i, j, :]      # K
                
                # Skip if all missing
                if all(ismissing, precip_ts) || all(ismissing, temp_ts)
                    continue
                end
                
                # Calculate CWD as cumulative PET - AET
                cwd_cumulative = 0.0  # mm
                
                for t in 1:ntime
                    if ismissing(precip_ts[t]) || ismissing(temp_ts[t])
                        continue
                    end
                    
                    # Convert precipitation from m/s to mm per timestep
                    P_mm = precip_ts[t] * dt_seconds * 1000.0  # m/s → mm
                    
                    # Estimate PET using simplified temperature-based formula
                    temp_C = temp_ts[t] - 273.15
                    if temp_C > 0
                        PET_mm_day = 0.05 * temp_C * (1.0 + 0.01 * temp_C)
                        PET_mm = PET_mm_day / 24.0  # Convert to hourly
                    else
                        PET_mm = 0.0
                    end
                    
                    # Estimate AET as minimum of PET and available water
                    # (simplified assumption: all precipitation is immediately available)
                    AET_mm = min(PET_mm, P_mm)
                    
                    # Accumulate deficit: CWD = PET - AET
                    cwd_cumulative += (PET_mm - AET_mm)
                end
                
                # Store this year's total CWD
                CWD_this_year[i, j] = cwd_cumulative
            end
        end
        
        # Accumulate across years
        CWD_sum .+= CWD_this_year
        year_count += 1
        
        close(ds)
    end
    
    # Average across all years to get climatological mean annual CWD
    CWD_spatial = CWD_sum ./ year_count
    
    # Convert to Float32 for storage efficiency
    CWD_spatial = Float32.(CWD_spatial)
    
    # Print statistics
    valid_cwd = filter(!isnan, CWD_spatial)
    @info "CWD statistics (mm/year) averaged over $year_count year(s):" min=minimum(valid_cwd) max=maximum(valid_cwd) mean=mean(valid_cwd) median=median(valid_cwd)
    
    return lons, lats, CWD_spatial
end

function regrid_cwd_to_model_space(
    lons_cwd::Vector{Float32},
    lats_cwd::Vector{Float32},
    CWD_spatial::Matrix{Float32},
    surface_space,
    FT,
)
    # Get coordinates from surface space
    coords = Fields.coordinate_field(surface_space)
    
    # Extract coordinate data to CPU using parent() and Array()
    # parent() gets the underlying data layout, then Array() moves to CPU
    lon_values = Array(parent(coords.long))
    lat_values = Array(parent(coords.lat))
    
    # Create interpolation object on CPU
    itp = Interpolations.interpolate(
        (lons_cwd, lats_cwd),
        CWD_spatial,
        Interpolations.Gridded(Interpolations.Linear())
    )
    
    # Add extrapolation with Flat boundary condition
    itp_extrap = Interpolations.extrapolate(itp, Interpolations.Flat())
    
    # Interpolate on CPU - element-wise over CPU arrays
    CWD_values_cpu = itp_extrap.(lon_values, lat_values)
    
    # Create a new Field from the CPU data
    # similar() creates a field with the same space and element type
    CWD_field = similar(coords.long, FT)
    
    # Transfer to GPU if needed, then fill
    # Use ClimaComms to detect device and transfer appropriately
    if parent(CWD_field) isa CUDA.CuArray
        # GPU backend: transfer CPU array to GPU
        CWD_values_gpu = CUDA.CuArray(FT.(CWD_values_cpu))
        parent(CWD_field) .= CWD_values_gpu
    else
        # CPU backend: direct assignment
        parent(CWD_field) .= FT.(CWD_values_cpu)
    end
    
    return CWD_field
end

# Main execution - run this once to generate CWD field
function main(start_year::Int = 2008, stop_year::Int = 2008)
    # Get the ERA5 artifact path for the specified year range
    # This will use the path from Overrides.toml on the HPC
    context = ClimaComms.context()
    
    start_date = DateTime(start_year, 1, 1)
    stop_date = DateTime(stop_year, 12, 31)
    
    era5_paths = ClimaLand.Artifacts.find_era5_year_paths(
        start_date,
        stop_date;
        context = context
    )
    
    @info "Found $(length(era5_paths)) ERA5 file(s) for $start_year-$stop_year"
    
    # Calculate CWD from ERA5
    lons, lats, CWD_spatial = calculate_cwd_from_era5(
        era5_paths,
        start_date,
        stop_date,
    )
    
    println("=" ^ 60)
    println("CWD calculated on ERA5 grid for $start_year-$stop_year")
    println("Shape: ", size(CWD_spatial))
    println("=" ^ 60)
    
    # Save the raw CWD data for later use
    # Use year range in filename
    if start_year == stop_year
        output_file = joinpath(@__DIR__, "cwd_era5_$(start_year).jld2")
    else
        output_file = joinpath(@__DIR__, "cwd_era5_$(start_year)_$(stop_year).jld2")
    end
    jldsave(output_file; lons, lats, CWD_spatial, start_year, stop_year)
    @info "Saved CWD to" output_file
    
    return lons, lats, CWD_spatial
end

# Run it
if abspath(PROGRAM_FILE) == @__FILE__
    lons, lats, CWD_spatial = main()
end