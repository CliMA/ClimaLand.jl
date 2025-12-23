"""
Compute Growing Season Length (GSL) from LAI and temperature, plus extract annual potential GPP.

GSL computation follows Zhou et al. (2025) definition with hybrid approach:
- For regions with clear LAI seasonality: GSL from LAI dynamics
- For non-seasonal regions (tropics, deserts): GSL from temperature (days with T > 0°C)

This avoids GSL = 0 which would break the optimal LAI model equations.

Output: NetCDF file with GSL (days) and A0_annual (mol CO2 m^-2 yr^-1) on (lon, lat) grid.
"""

using NCDatasets
using Statistics

# Configuration
LAI_FILE = joinpath(@__DIR__, "lai_1M_average.nc")
A0_FILE = joinpath(@__DIR__, "a0_annual_1M_average.nc")
TAIR_FILE = joinpath(@__DIR__, "tair_1M_average.nc")
OUTPUT_FILE = joinpath(@__DIR__, "gsl_a0_annual.nc")

# GSL parameters
LAI_THRESHOLD_FRACTION = 0.2  # LAI must be above 20% of annual max to count as "growing"
SEASONALITY_CV_THRESHOLD = 0.15  # Coefficient of variation below this = non-seasonal
FREEZING_TEMP = 273.15  # K (0°C) - Zhou et al. threshold for growing season
MIN_GSL = 30.0  # Minimum GSL in days (for numerical stability in ice-covered regions)

function compute_gsl_from_temperature(tair_timeseries::Vector{T}) where T
    """
    Compute GSL from monthly temperature following Zhou et al. (2025):
    GSL = number of days with T > 0°C.

    For monthly data, we count months above freezing and convert to days.
    """
    # Handle NaN
    valid = .!isnan.(tair_timeseries)
    if !any(valid)
        return T(NaN)
    end

    # If we have 24 months (2 years), compute mean annual cycle
    if length(tair_timeseries) == 24
        tair_monthly = [mean([tair_timeseries[i], tair_timeseries[i+12]]) for i in 1:12]
    elseif length(tair_timeseries) == 12
        tair_monthly = tair_timeseries
    else
        tair_monthly = tair_timeseries[1:min(12, length(tair_timeseries))]
    end

    # Replace NaN with very cold temp (won't count as growing)
    tair_monthly = [isnan(x) ? T(200) : x for x in tair_monthly]

    # Count months above freezing
    months_above_freezing = sum(tair_monthly .> FREEZING_TEMP)

    # Convert to days
    gsl = months_above_freezing * T(365.25 / 12)

    return gsl
end

function compute_gsl_from_lai(lai_timeseries::Vector{T};
                               threshold_frac=LAI_THRESHOLD_FRACTION,
                               cv_threshold=SEASONALITY_CV_THRESHOLD) where T
    """
    Compute GSL from monthly LAI time series.
    Returns (gsl, is_seasonal) tuple.

    Args:
        lai_timeseries: Vector of monthly LAI values
        threshold_frac: Fraction of max LAI above which counts as growing season
        cv_threshold: Below this coefficient of variation, region is non-seasonal

    Returns:
        (gsl_days, is_seasonal): GSL in days and whether LAI shows seasonality
    """
    # Handle NaN - if all NaN, return NaN
    valid = .!isnan.(lai_timeseries)
    if !any(valid)
        return T(NaN), false
    end

    # If we have 24 months (2 years), compute mean annual cycle
    if length(lai_timeseries) == 24
        lai_monthly = [mean([lai_timeseries[i], lai_timeseries[i+12]]) for i in 1:12]
    elseif length(lai_timeseries) == 12
        lai_monthly = lai_timeseries
    else
        lai_monthly = lai_timeseries[1:min(12, length(lai_timeseries))]
    end

    # Replace NaN with 0 for calculations
    lai_monthly = [isnan(x) ? T(0) : x for x in lai_monthly]

    # Check for seasonality using coefficient of variation
    lai_mean = mean(lai_monthly)
    lai_std = std(lai_monthly)

    # If mean LAI is essentially zero, not seasonal (desert/ice)
    if lai_mean < T(0.01)
        return T(0), false  # Will use temperature-based GSL
    end

    cv = lai_std / lai_mean

    # Non-seasonal regions
    if cv < cv_threshold
        return T(365), false  # Will use temperature-based GSL
    end

    # Seasonal: compute threshold-based GSL
    lai_max = maximum(lai_monthly)
    threshold = threshold_frac * lai_max
    months_above = sum(lai_monthly .> threshold)
    gsl = months_above * T(365.25 / 12)

    return gsl, true
end

function compute_hybrid_gsl(lai_ts::Vector{T}, tair_ts::Vector{T}) where T
    """
    Hybrid GSL computation:
    - If LAI shows clear seasonality → use LAI-based GSL
    - Otherwise → use temperature-based GSL (Zhou et al. definition)

    This ensures:
    - Tropics get GSL ≈ 365 (warm year-round)
    - Deserts get GSL ≈ 365 (warm year-round, just water-limited)
    - Tundra gets short GSL based on temperature
    - Temperate gets LAI-based GSL reflecting actual phenology
    """
    lai_gsl, is_seasonal = compute_gsl_from_lai(lai_ts)

    if is_seasonal && !isnan(lai_gsl) && lai_gsl > 0
        # Use LAI-based GSL for regions with clear seasonality
        return lai_gsl
    else
        # Use temperature-based GSL for non-seasonal or no-vegetation regions
        temp_gsl = compute_gsl_from_temperature(tair_ts)
        return temp_gsl
    end
end

function main()
    println("Loading LAI data from: $LAI_FILE")
    ds_lai = NCDataset(LAI_FILE)
    lai = ds_lai["lai"][:, :, :]  # (time, lon, lat)
    lon = ds_lai["lon"][:]
    lat = ds_lai["lat"][:]
    close(ds_lai)

    println("Loading A0_annual data from: $A0_FILE")
    ds_a0 = NCDataset(A0_FILE)
    a0_annual = ds_a0["a0_annual"][:, :, :]  # (time, lon, lat)
    close(ds_a0)

    println("Loading temperature data from: $TAIR_FILE")
    ds_tair = NCDataset(TAIR_FILE)
    tair = ds_tair["tair"][:, :, :]  # (time, lon, lat)
    close(ds_tair)

    nlon, nlat = length(lon), length(lat)
    ntime = size(lai, 1)
    println("Grid size: $nlon x $nlat, $ntime time steps")

    # Allocate output arrays
    gsl = fill(NaN, nlon, nlat)
    a0_out = fill(NaN, nlon, nlat)

    # Track which method was used
    n_lai_based = 0
    n_temp_based = 0

    # Compute GSL and A0 for each grid cell
    println("Computing GSL (hybrid) and A0_annual...")
    for i in 1:nlon
        for j in 1:nlat
            # Extract time series for this pixel
            lai_ts = Float64.(lai[:, i, j])
            a0_ts = Float64.(a0_annual[:, i, j])
            tair_ts = Float64.(tair[:, i, j])

            # Compute hybrid GSL
            lai_gsl, is_seasonal = compute_gsl_from_lai(lai_ts)
            if is_seasonal && !isnan(lai_gsl) && lai_gsl > 0
                gsl[i, j] = max(lai_gsl, MIN_GSL)
                n_lai_based += 1
            else
                gsl[i, j] = max(compute_gsl_from_temperature(tair_ts), MIN_GSL)
                n_temp_based += 1
            end

            # Use last time point for A0_annual (model spins up from uniform initial value)
            a0_last = a0_ts[end]
            if !isnan(a0_last)
                a0_out[i, j] = a0_last
            end
        end
    end

    # Statistics
    valid_gsl = filter(!isnan, vec(gsl))
    valid_a0 = filter(!isnan, vec(a0_out))

    println("\nGSL computation method:")
    println("  LAI-based (seasonal): $n_lai_based pixels")
    println("  Temperature-based (non-seasonal): $n_temp_based pixels")

    println("\nGSL statistics:")
    println("  Valid points: $(length(valid_gsl))")
    println("  Range: $(minimum(valid_gsl)) - $(maximum(valid_gsl)) days")
    println("  Mean: $(round(mean(valid_gsl), digits=1)) days")

    # Check for zeros
    n_zeros = sum(valid_gsl .== 0)
    println("  Points with GSL=0: $n_zeros")

    println("\nA0_annual statistics:")
    println("  Valid points: $(length(valid_a0))")
    println("  Range: $(round(minimum(valid_a0), digits=1)) - $(round(maximum(valid_a0), digits=1)) mol CO2 m^-2 yr^-1")
    println("  Mean: $(round(mean(valid_a0), digits=1)) mol CO2 m^-2 yr^-1")

    # Write output NetCDF
    println("\nWriting output to: $OUTPUT_FILE")
    ds_out = NCDataset(OUTPUT_FILE, "c")

    # Define dimensions
    defDim(ds_out, "lon", nlon)
    defDim(ds_out, "lat", nlat)

    # Define coordinate variables
    lon_var = defVar(ds_out, "lon", Float64, ("lon",))
    lon_var.attrib["units"] = "degrees_east"
    lon_var.attrib["long_name"] = "longitude"
    lon_var[:] = lon

    lat_var = defVar(ds_out, "lat", Float64, ("lat",))
    lat_var.attrib["units"] = "degrees_north"
    lat_var.attrib["long_name"] = "latitude"
    lat_var[:] = lat

    # Define data variables
    gsl_var = defVar(ds_out, "gsl", Float64, ("lon", "lat"))
    gsl_var.attrib["units"] = "days"
    gsl_var.attrib["long_name"] = "Growing Season Length"
    gsl_var.attrib["description"] = "Hybrid GSL following Zhou et al. (2025): LAI-based for seasonal regions (CV > $(SEASONALITY_CV_THRESHOLD)), temperature-based (days with T > 0C) for non-seasonal regions. Minimum GSL = $(MIN_GSL) days for numerical stability."
    gsl_var[:, :] = gsl

    a0_var = defVar(ds_out, "a0_annual", Float64, ("lon", "lat"))
    a0_var.attrib["units"] = "mol CO2 m^-2 yr^-1"
    a0_var.attrib["long_name"] = "Annual Potential GPP"
    a0_var.attrib["description"] = "Annual potential GPP computed with fAPAR=1 and beta=1 (no moisture stress). Last time point from simulation (after spin-up)."
    a0_var[:, :] = a0_out

    # Global attributes
    ds_out.attrib["title"] = "Growing Season Length and Annual Potential GPP"
    ds_out.attrib["source"] = "Computed from ClimaLand.jl optimal LAI simulation"
    ds_out.attrib["history"] = "Created by compute_gsl_a0.jl with hybrid GSL method"
    ds_out.attrib["references"] = "Zhou et al. (2025) Global Change Biology - GSL defined as days with T > 0C"

    close(ds_out)
    println("Done!")
end

main()
