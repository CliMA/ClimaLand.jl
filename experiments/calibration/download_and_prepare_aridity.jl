#!/usr/bin/env julia
"""
Download and prepare Global Aridity Index for use in ClimaLand calibration pipeline

This script:
1. Downloads the reference Global Aridity Index from Trabucco & Zomer (2019)
2. Converts from GeoTIFF to a format we can work with
3. Regrids to the ClimaLand spectral element grid
4. Saves as a ClimaCore Field compatible with the calibration pipeline
5. Replaces the ERA5-derived CWD with reference aridity

Reference:
Trabucco, A., and Zomer, R.J. (2019) Global Aridity Index and Potential 
Evapotranspiration (ET0) Climate Database v2. figshare. Dataset.
https://doi.org/10.6084/m9.figshare.7504448.v3
"""

using Downloads
using ArchGDAL
using JLD2
using ClimaCore
using ClimaCore: Fields
import ClimaComms
import ClimaLand
using Statistics
using Interpolations
using Dates

# Try to load CUDA, but don't fail if it's not available
try
    import CUDA
catch
    @warn "CUDA not available, using CPU only"
end

# Configuration
FT = Float32
nelements = (101, 15)

function download_aridity_index()
    """Download the Global Aridity Index GeoTIFF if not already present"""
    
    aridity_file = joinpath(@__DIR__, "Global-AI_ET0_v3_annual.tif")
    
    if isfile(aridity_file)
        @info "✓ Aridity Index file already exists" aridity_file
        return aridity_file
    end
    
    @error "="^70
    @error "Global Aridity Index file not found!"
    @error "="^70
    @error ""
    @error "Please download manually:"
    @error "1. Visit: https://figshare.com/articles/dataset/7504448"
    @error "2. Click on 'Global-AI_ET0_v3_annual.tif'"  
    @error "3. Click 'Download' button"
    @error "4. Save to: $(dirname(aridity_file))"
    @error ""
    @error "OR use browser/curl:"
    @error "  Open: https://figshare.com/ndownloader/files/34377245"
    @error "  Save the downloaded file as 'Global-AI_ET0_v3_annual.tif'"
    @error ""
    @error "="^70
    
    error("Aridity Index file not found: $aridity_file")
end

function read_aridity_geotiff(filepath::String)
    """Read aridity index from GeoTIFF and return lons, lats, aridity"""
    
    @info "Reading GeoTIFF with ArchGDAL..."
    
    ArchGDAL.read(filepath) do dataset
        # Read the raster band
        band = ArchGDAL.getband(dataset, 1)
        aridity_data = ArchGDAL.read(band)
        
        # Get geotransform for coordinates
        transform = ArchGDAL.getgeotransform(dataset)
        
        # Extract grid information
        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)
        
        @info "Raster dimensions" width height
        @info "Geotransform" transform
        
        # Calculate lat/lon arrays
        # Geotransform: [x_origin, pixel_width, 0, y_origin, 0, pixel_height]
        x_origin = transform[1]
        pixel_width = transform[2]
        y_origin = transform[4]
        pixel_height = transform[6]  # Usually negative
        
        lons = Float32[x_origin + (i - 0.5) * pixel_width for i in 1:width]
        lats = Float32[y_origin + (j - 0.5) * pixel_height for j in 1:height]
        
        # Handle nodata values
        nodata = ArchGDAL.getnodatavalue(band)
        @info "NoData value" nodata
        
        # Convert to Float32 and handle nodata
        aridity = Float32.(aridity_data)
        if nodata !== nothing
            aridity[aridity .== Float32(nodata)] .= NaN32
        end
        
        # Data is scaled: divide by 10000 to get actual aridity index values
        # (This GeoTIFF uses integer encoding: value = aridity * 10000)
        aridity = aridity ./ 10000.0f0
        
        # Aridity values should be 0-10+ (unitless, P/ET0 ratio)
        # Negative values or very large values are invalid
        aridity[aridity .< 0] .= NaN32
        aridity[aridity .> 100] .= NaN32  # Unrealistic values
        
        @info "Aridity range" min=minimum(filter(!isnan, aridity)) max=maximum(filter(!isnan, aridity))
        @info "NaN count" sum(isnan.(aridity))
        
        return lons, lats, aridity
    end
end

function regrid_aridity_to_model_space(
    lons_aridity::Vector{Float32},
    lats_aridity::Vector{Float32},
    aridity_data::Matrix{Float32},
    surface_space,
    FT,
)
    """Regrid aridity from lat/lon grid to ClimaLand spectral element grid"""
    
    @info "Regridding aridity to model spectral element grid..."
    
    # Get coordinates from surface space
    coords = Fields.coordinate_field(surface_space)
    
    # Extract coordinate data to CPU
    lon_values = Array(parent(coords.long))
    lat_values = Array(parent(coords.lat))
    
    # Convert longitudes if needed (GeoTIFF might be 0-360, we need -180 to 180)
    if minimum(lons_aridity) >= 0 && maximum(lons_aridity) > 180
        @info "Converting longitudes from 0-360 to -180-180"
        lons_aridity = ifelse.(lons_aridity .> 180, lons_aridity .- 360, lons_aridity)
        # Need to reorder data to match
        sorted_idx = sortperm(lons_aridity)
        lons_aridity = lons_aridity[sorted_idx]
        aridity_data = aridity_data[sorted_idx, :]
    end
    
    # Ensure latitudes are sorted (should go from -90 to +90 or +90 to -90)
    if lats_aridity[1] > lats_aridity[end]
        @info "Reversing latitude order (top-to-bottom → bottom-to-top)"
        reverse!(lats_aridity)
        aridity_data = reverse(aridity_data, dims=2)
    end
    
    # IMPORTANT: Keep NaN for ocean points - they have meaning!
    # Aridity = 0 means hyperarid land (Sahara), NOT ocean
    # Ocean points are NaN in the GeoTIFF and should stay NaN
    # 
    # CRITICAL: Linear interpolation near NaN values produces NaN at nearby points!
    # This causes coastal land points to incorrectly become NaN.
    # Solution: Use nearest-neighbor interpolation to avoid NaN contamination at boundaries
    
    @info "Creating interpolation object (using nearest-neighbor to avoid NaN contamination)..."
    
    # Use nearest-neighbor instead of linear to prevent NaN from spreading
    # to coastal land points during interpolation
    itp = Interpolations.interpolate(
        (lons_aridity, lats_aridity),
        aridity_data,
        Interpolations.Gridded(Interpolations.Constant())  # Nearest-neighbor
    )
    
    # Add extrapolation - NaN will be used for points outside the data bounds
    itp_extrap = Interpolations.extrapolate(itp, NaN32)
    
    @info "Interpolating to model coordinates..."
    
    # Interpolate on CPU - element-wise over CPU arrays
    aridity_values_cpu = itp_extrap.(lon_values, lat_values)
    
    # Create a new Field from the CPU data
    aridity_field = similar(coords.long, FT)
    
    # Transfer to GPU if needed, then fill
    if isdefined(Main, :CUDA) && parent(aridity_field) isa CUDA.CuArray
        # GPU backend: transfer CPU array to GPU
        aridity_values_gpu = CUDA.CuArray(FT.(aridity_values_cpu))
        parent(aridity_field) .= aridity_values_gpu
    else
        # CPU backend: direct assignment
        parent(aridity_field) .= FT.(aridity_values_cpu)
    end
    
    @info "Regridding complete"
    
    return aridity_field
end

function main()
    @info "="^70
    @info "Global Aridity Index Download and Preparation"
    @info "="^70
    
    # Step 1: Download if needed
    aridity_file = download_aridity_index()
    
    # Step 2: Read the GeoTIFF
    lons_aridity, lats_aridity, aridity_data = read_aridity_geotiff(aridity_file)
    
    # Step 3: Create model grid
    @info "Creating ClimaLand spectral element grid..."
    context = ClimaComms.context()
    domain = ClimaLand.Domains.global_domain(FT; nelements)
    surface_space = domain.space.surface
    
    # Step 4: Regrid to model space
    aridity_field = regrid_aridity_to_model_space(
        lons_aridity, lats_aridity, aridity_data, surface_space, FT
    )
    
    # Step 5: Print statistics
    aridity_values = Array(parent(aridity_field))
    n_total = length(aridity_values)
    n_nan = count(isnan, aridity_values)
    n_valid = n_total - n_nan
    valid_values = filter(!isnan, aridity_values)
    
    @info "="^70
    @info "Aridity Field Statistics (on model grid)"
    @info "  Total points: $(n_total)"
    @info "  Ocean points (NaN): $(n_nan) ($(round(100*n_nan/n_total, digits=1))%)"
    @info "  Land points (valid): $(n_valid) ($(round(100*n_valid/n_total, digits=1))%)"
    if n_valid > 0
        @info "  Land aridity range: [$(minimum(valid_values)), $(maximum(valid_values))]"
        @info "  Land aridity mean: $(mean(valid_values))"
        @info "  Land aridity median: $(median(valid_values))"
    end
    @info "="^70
    
    # Step 6: Save in format compatible with calibration pipeline
    # Save as aridity.jld2 (Global Aridity Index where higher = wetter)
    output_file = joinpath(@__DIR__, "aridity.jld2")
    backup_file = joinpath(@__DIR__, "aridity_backup.jld2")
    
    # Backup old file if it exists
    if isfile(output_file)
        @info "Backing up old aridity file to aridity_backup.jld2"
        cp(output_file, backup_file, force=true)
    end
    
    # Save with same structure as before for compatibility
    JLD2.jldsave(
        output_file;
        lons = lons_aridity,
        lats = lats_aridity,
        CWD_field = aridity_field,  # Actually aridity, but keeping name for compatibility
        start_year = 1970,  # Reference period for Global AI
        stop_year = 2000,
        created = now(),
        description = "Global Aridity Index (Trabucco & Zomer 2019) regridded to model grid. NOTE: This is aridity (P/ET0). Higher values = more humid."
    )
    
    @info "✓ Saved aridity field to" output_file
    @info "  File size: $(filesize(output_file) / 1e6) MB"
    
    @info "="^70
    @info "IMPORTANT NOTES:"
    @info "="^70
    @info "1. This file contains ARIDITY INDEX (P/ET0), not CWD"
    @info "   - Higher aridity = MORE water available (wetter)"
    @info "   - CWD was opposite: higher CWD = LESS water (drier)"
    @info ""
    @info "2. Ocean masking:"
    @info "   - Ocean points are NaN (NOT zero!)"
    @info "   - Aridity = 0 means hyperarid land (Sahara, Atacama)"
    @info "   - Zero is a valid data value, NaN indicates no land"
    @info ""
    @info "3. Classification (for LAND points only):"
    @info "   - Aridity < 0.05: Hyper Arid"
    @info "   - Aridity 0.05-0.20: Arid"
    @info "   - Aridity 0.20-0.50: Semi-Arid"
    @info "   - Aridity 0.50-0.65: Dry Sub-humid"
    @info "   - Aridity > 0.65: Humid"
    @info "="^70
    
    return aridity_field
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    aridity_field = main()
end
