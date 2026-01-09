"""
Compute Growing Season Length (GSL) from LAI and temperature, plus extract annual potential GPP,
mean annual precipitation, and growing season VPD.

GSL computation follows Zhou et al. (2025) definition with hybrid approach:
- For regions with clear LAI seasonality: GSL from LAI dynamics
- For non-seasonal regions (tropics, deserts): GSL from temperature (days with T > 0°C)

This avoids GSL = 0 which would break the optimal LAI model equations.

Additional variables for water limitation in the optimal LAI model:
- precip_annual: Mean annual precipitation (mm yr⁻¹)
- vpd_gs: Average VPD during growing season (Pa)

These are needed for the f₀×P/A₀ water limitation term in LAI_max (Zhou et al. 2025, Eq. 11).

Output: NetCDF file with GSL (days), A0_annual (mol CO2 m^-2 yr^-1), precip_annual (mm yr^-1),
        and vpd_gs (Pa) on (lon, lat) grid.
"""

using NCDatasets
using Statistics

# Configuration
LAI_FILE = joinpath(@__DIR__, "lai_1M_average.nc")
A0_FILE = joinpath(@__DIR__, "a0_annual_1M_average.nc")
TAIR_FILE = joinpath(@__DIR__, "tair_1M_average.nc")
PRECIP_FILE = joinpath(@__DIR__, "precip_1M_average.nc")
VPD_FILE = joinpath(@__DIR__, "vpd_1M_average.nc")
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

function compute_annual_precipitation(precip_timeseries::Vector{T}) where T
    """
    Compute mean annual precipitation from monthly precipitation rates.

    Args:
        precip_timeseries: Vector of monthly precipitation rates (kg m⁻² s⁻¹)
            Note: ERA5 often uses negative values for precipitation (downward flux convention).

    Returns:
        precip_annual: Mean annual precipitation (mm yr⁻¹), always positive

    Note: 1 kg m⁻² = 1 mm water depth
    """
    # Handle NaN
    valid = .!isnan.(precip_timeseries)
    if !any(valid)
        return T(NaN)
    end

    # If we have 24 months (2 years), compute mean annual cycle
    if length(precip_timeseries) == 24
        precip_monthly = [mean([precip_timeseries[i], precip_timeseries[i+12]]) for i in 1:12]
    elseif length(precip_timeseries) == 12
        precip_monthly = precip_timeseries
    else
        precip_monthly = precip_timeseries[1:min(12, length(precip_timeseries))]
    end

    # Replace NaN with 0 for calculations and take absolute value
    # (ERA5 uses negative values for precipitation by convention)
    precip_monthly = [isnan(x) ? T(0) : abs(x) for x in precip_monthly]

    # Convert from kg m⁻² s⁻¹ to mm yr⁻¹
    # Sum monthly values, each multiplied by seconds in that month
    # For simplicity, use average seconds per month = 365.25 * 24 * 3600 / 12
    seconds_per_month = T(365.25 * 24 * 3600 / 12)

    # Annual precipitation = sum of (monthly_rate × seconds_per_month)
    # Since kg m⁻² = mm, this gives mm yr⁻¹
    precip_annual = sum(precip_monthly) * seconds_per_month

    return precip_annual
end

function get_growing_season_mask(lai_ts::Vector{T}, tair_ts::Vector{T};
                                  threshold_frac=LAI_THRESHOLD_FRACTION,
                                  cv_threshold=SEASONALITY_CV_THRESHOLD) where T
    """
    Determine which months are in the growing season.

    Returns a 12-element boolean vector indicating growing season months.
    Uses LAI-based criterion for seasonal regions, temperature-based for non-seasonal.
    """
    # Get mean annual cycle
    if length(lai_ts) == 24
        lai_monthly = [mean([lai_ts[i], lai_ts[i+12]]) for i in 1:12]
        tair_monthly = [mean([tair_ts[i], tair_ts[i+12]]) for i in 1:12]
    elseif length(lai_ts) == 12
        lai_monthly = lai_ts
        tair_monthly = tair_ts
    else
        lai_monthly = lai_ts[1:min(12, length(lai_ts))]
        tair_monthly = tair_ts[1:min(12, length(tair_ts))]
    end

    # Replace NaN with 0/very cold for calculations
    lai_monthly = [isnan(x) ? T(0) : x for x in lai_monthly]
    tair_monthly = [isnan(x) ? T(200) : x for x in tair_monthly]

    # Check for seasonality using coefficient of variation
    lai_mean = mean(lai_monthly)
    lai_std = std(lai_monthly)

    # If mean LAI is essentially zero (desert/ice), use temperature
    if lai_mean < T(0.01)
        return tair_monthly .> FREEZING_TEMP
    end

    cv = lai_std / lai_mean

    # Non-seasonal regions: use temperature criterion
    if cv < cv_threshold
        return tair_monthly .> FREEZING_TEMP
    end

    # Seasonal: use LAI threshold
    lai_max = maximum(lai_monthly)
    threshold = threshold_frac * lai_max
    return lai_monthly .> threshold
end

function compute_vpd_growing_season(vpd_timeseries::Vector{T}, lai_ts::Vector{T}, tair_ts::Vector{T}) where T
    """
    Compute average VPD during growing season.

    Args:
        vpd_timeseries: Vector of monthly VPD values (Pa)
        lai_ts: Vector of monthly LAI values (for determining growing season)
        tair_ts: Vector of monthly temperature values (K)

    Returns:
        vpd_gs: Mean VPD during growing season (Pa)
    """
    # Handle NaN
    valid = .!isnan.(vpd_timeseries)
    if !any(valid)
        return T(NaN)
    end

    # Get mean annual cycle for VPD
    if length(vpd_timeseries) == 24
        vpd_monthly = [mean([vpd_timeseries[i], vpd_timeseries[i+12]]) for i in 1:12]
    elseif length(vpd_timeseries) == 12
        vpd_monthly = vpd_timeseries
    else
        vpd_monthly = vpd_timeseries[1:min(12, length(vpd_timeseries))]
    end

    # Get growing season mask
    gs_mask = get_growing_season_mask(lai_ts, tair_ts)

    # Replace NaN with 0 for calculations
    vpd_monthly = [isnan(x) ? T(0) : x for x in vpd_monthly]

    # Compute mean VPD during growing season
    gs_vpd = vpd_monthly[gs_mask]
    if isempty(gs_vpd)
        # No growing season months (very cold region), return annual mean
        return mean(vpd_monthly)
    end

    return mean(gs_vpd)
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

    println("Loading precipitation data from: $PRECIP_FILE")
    ds_precip = NCDataset(PRECIP_FILE)
    precip = ds_precip["precip"][:, :, :]  # (time, lon, lat) in kg m⁻² s⁻¹
    close(ds_precip)

    println("Loading VPD data from: $VPD_FILE")
    ds_vpd = NCDataset(VPD_FILE)
    vpd = ds_vpd["vpd"][:, :, :]  # (time, lon, lat) in Pa
    close(ds_vpd)

    nlon, nlat = length(lon), length(lat)
    ntime = size(lai, 1)
    println("Grid size: $nlon x $nlat, $ntime time steps")

    # Allocate output arrays
    gsl = fill(NaN, nlon, nlat)
    a0_out = fill(NaN, nlon, nlat)
    precip_annual = fill(NaN, nlon, nlat)
    vpd_gs = fill(NaN, nlon, nlat)

    # Track which method was used
    n_lai_based = 0
    n_temp_based = 0

    # Compute GSL, A0, precip_annual, and vpd_gs for each grid cell
    println("Computing GSL (hybrid), A0_annual, precip_annual, and vpd_gs...")
    for i in 1:nlon
        for j in 1:nlat
            # Extract time series for this pixel
            lai_ts = Float64.(lai[:, i, j])
            a0_ts = Float64.(a0_annual[:, i, j])
            tair_ts = Float64.(tair[:, i, j])
            precip_ts = Float64.(precip[:, i, j])
            vpd_ts = Float64.(vpd[:, i, j])

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

            # Compute mean annual precipitation
            precip_annual[i, j] = compute_annual_precipitation(precip_ts)

            # Compute average VPD during growing season
            vpd_gs[i, j] = compute_vpd_growing_season(vpd_ts, lai_ts, tair_ts)
        end
    end

    # Statistics
    valid_gsl = filter(!isnan, vec(gsl))
    valid_a0 = filter(!isnan, vec(a0_out))
    valid_precip = filter(!isnan, vec(precip_annual))
    valid_vpd = filter(!isnan, vec(vpd_gs))

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

    println("\nPrecip_annual statistics:")
    println("  Valid points: $(length(valid_precip))")
    println("  Range: $(round(minimum(valid_precip), digits=1)) - $(round(maximum(valid_precip), digits=1)) mm yr^-1")
    println("  Mean: $(round(mean(valid_precip), digits=1)) mm yr^-1")

    println("\nVPD_gs (growing season) statistics:")
    println("  Valid points: $(length(valid_vpd))")
    println("  Range: $(round(minimum(valid_vpd), digits=1)) - $(round(maximum(valid_vpd), digits=1)) Pa")
    println("  Mean: $(round(mean(valid_vpd), digits=1)) Pa")

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

    precip_var = defVar(ds_out, "precip_annual", Float64, ("lon", "lat"))
    precip_var.attrib["units"] = "mm yr^-1"
    precip_var.attrib["long_name"] = "Mean Annual Precipitation"
    precip_var.attrib["description"] = "Mean annual precipitation computed from monthly averages. Used for water limitation term (f0*P/A0) in LAI_max calculation following Zhou et al. (2025)."
    precip_var[:, :] = precip_annual

    vpd_var = defVar(ds_out, "vpd_gs", Float64, ("lon", "lat"))
    vpd_var.attrib["units"] = "Pa"
    vpd_var.attrib["long_name"] = "Growing Season Vapor Pressure Deficit"
    vpd_var.attrib["description"] = "Average VPD during growing season months. Growing season is determined using the same hybrid LAI/temperature criterion as GSL. Used for water limitation term in LAI_max calculation following Zhou et al. (2025)."
    vpd_var[:, :] = vpd_gs

    # Global attributes
    ds_out.attrib["title"] = "Growing Season Length, Annual Potential GPP, Precipitation, and VPD for Optimal LAI Model"
    ds_out.attrib["source"] = "Computed from ClimaLand.jl optimal LAI simulation and ERA5 climate data"
    ds_out.attrib["history"] = "Created by compute_gsl_a0.jl with hybrid GSL method"
    ds_out.attrib["references"] = "Zhou et al. (2025) Global Change Biology - GSL defined as days with T > 0C, water limitation via f0*P/A0 term"

    close(ds_out)
    println("Done!")
end

main()
