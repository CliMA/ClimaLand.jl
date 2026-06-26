# Site-level "drivers" for the per-site calibration CSV.
#
# Reads the fluxnet forcing CSV directly (not via ClimaLand's TVI pipeline) so
# the aggregation matches what one would compute by hand from the half-hourly
# record. Spatial fields (porosity, MODIS max LAI, IGBP) are looked up at the
# site's lat/lon from the same artifacts the forward model uses.
#
# Returns a NamedTuple consumed by `append_site_csv` in calibrate_site.jl.

import NCDatasets: NCDataset
import Dates: DateTime, year
import Statistics: mean
import DelimitedFiles: readdlm

const _FLUXNET_MISSING = -9999.0

"""
    _read_fluxnet_column(site_ID, col_name) -> (utc_datetimes, values_vec)

Read one column from the half-hourly fluxnet CSV (the same file ClimaLand's
forcing pipeline reads), returning UTC timestamps and a Float64 vector with
fluxnet missing values (-9999) replaced by `NaN`. Returns `nothing` if the
column is absent so callers can fall back gracefully.
"""
function _read_fluxnet_column(site_ID, time_offset, col_name)
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)
    col_idx = findfirst(vec(columns) .== col_name)
    col_idx === nothing && return nothing
    ts_idx = findfirst(vec(columns) .== "TIMESTAMP_START")
    te_idx = findfirst(vec(columns) .== "TIMESTAMP_END")

    raw = Float64.(data[:, col_idx])
    raw[raw .== _FLUXNET_MISSING] .= NaN

    function _to_utc(ts)
        local_dt = DateTime(string(Int64(ts)), "yyyymmddHHMM")
        return local_dt + Dates.Hour(time_offset)
    end
    starts = _to_utc.(data[:, ts_idx])
    ends = _to_utc.(data[:, te_idx])
    midpoints = starts .+ (ends .- starts) ./ 2
    return midpoints, raw
end

"""
    _aggregate_in_windows(times, vals, windows; reducer)

Apply `reducer` (e.g. `mean`, `sum`) to `vals` over the union of `windows`
(each `(start, stop)`), skipping NaNs. Returns the reducer's result, or NaN
if no valid points fall inside any window.
"""
function _aggregate_in_windows(times, vals, windows; reducer)
    mask = falses(length(times))
    for (ws, we) in windows
        @. mask = mask | ((times >= ws) & (times < we))
    end
    sel = vals[mask]
    sel = sel[.!isnan.(sel)]
    isempty(sel) && return NaN
    return reducer(sel)
end

"""
    _porosity_at_site(lat, long) -> Float64

Look up surface (top layer) porosity at the site lat/lon from the Gupta et al.
(2020) NetCDF used by ClimaLand's spatially-varying soil parameters. Returns
the global default 0.47 if anything goes wrong (artifact missing, point at
ocean mask, etc.) so the CSV always has a number.
"""
function _porosity_at_site(lat, long)
    # Use the highres (0.1°) Gupta map — it's the artifact the simulation
    # itself downloads, so we never have to fetch a second one just for the
    # CSV row. Variable name in the NetCDF is the literal "ν".
    try
        soil_params_dir =
            ClimaLand.Artifacts.soil_params_artifact_folder_path(; lowres = false)
        nc_path = joinpath(
            soil_params_dir,
            "porosity_map_gupta_etal2020_0.1x0.1x4.nc",
        )
        return _nc_value_at_site(nc_path, "ν", lat, long; default = 0.47)
    catch err
        @warn "Could not read porosity at site, using default" exception = err
        return 0.47
    end
end

"""
    _nc_value_at_site(nc_path, varname, lat, long; default)

Open `nc_path`, find the closest grid cell to (`lat`, `long`), and return the
value of `varname` at that cell (top layer if 3D). Returns `default` if the
variable is missing or the value is NaN/missing.
"""
function _nc_value_at_site(nc_path, varname, lat, long; default = NaN)
    NCDataset(nc_path) do ds
        haskey(ds, varname) || return default
        lats = ds["lat"][:]
        lons = ds["lon"][:]
        i = argmin(abs.(lons .- long))
        j = argmin(abs.(lats .- lat))
        v = ds[varname]
        val = if ndims(v) == 2
            v[i, j]
        elseif ndims(v) == 3
            v[i, j, 1]   # top layer
        else
            v[i, j, 1, 1]
        end
        (val === missing || (val isa Number && isnan(val))) && return default
        return Float64(val)
    end
end

"""
    compute_site_drivers(site_ID, coords, samples; soc_default = 5.0)

Aggregate fluxnet forcing variables, MODIS max LAI, soil porosity, IGBP class
and a placeholder SOC value over the union of the calibration `samples`
(year-long windows). All means / sums are computed by skipping NaN
(fluxnet -9999) values. Returns a NamedTuple suitable for the per-site CSV.

`soc_default` is the constant SOC initial condition the model uses
(`Y.soilco2.SOC = 5.0 kg C m⁻³` in `src/simulations/initial_conditions.jl`).
Reported here so the CSV has a column even though the value is currently
constant across sites — once a spatial SOC initial-condition map lands in
ClimaLand, swap this for a real lookup.
"""
function compute_site_drivers(site_ID, coords, samples; soc_default = 5.0)
    time_offset = coords.time_offset
    lat = Float64(coords.lat)
    long = Float64(coords.long)
    n_years = max(length(samples), 1)

    # Fluxnet aggregations
    function _mean_col(col, transform = identity)
        ret = _read_fluxnet_column(site_ID, time_offset, col)
        ret === nothing && return NaN
        (times, vals) = ret
        vals = transform.(vals)
        return _aggregate_in_windows(times, vals, samples; reducer = mean)
    end

    mean_TA = _mean_col("TA_F", x -> x + 273.15)         # °C → K
    mean_VPD = _mean_col("VPD_F", x -> x * 100.0)         # hPa → Pa
    mean_SW = _mean_col("SW_IN_F")                        # W m⁻²
    mean_TS_top = _mean_col("TS_F_MDS_1", x -> x + 273.15)

    # Annual cumulative precip: sum(P_F mm/half-hour) divided by year count
    annual_cum_P = let ret = _read_fluxnet_column(site_ID, time_offset, "P_F")
        if ret === nothing
            NaN
        else
            (times, vals) = ret
            total = _aggregate_in_windows(times, vals, samples; reducer = sum)
            total / n_years
        end
    end

    # Max LAI from MODIS at site (use first sample's mid-year as the lookup
    # year — MODIS climatology is roughly year-agnostic via the site grid cell)
    max_LAI = try
        ref_date =
            isempty(samples) ? DateTime(2010, 7, 1) :
            samples[1][1] + (samples[1][2] - samples[1][1]) ÷ 2
        FluxnetSimulations.get_maxLAI_at_site(ref_date, lat, long)
    catch err
        @warn "Could not read max LAI for $site_ID" exception = err
        NaN
    end

    porosity = _porosity_at_site(lat, long)

    igbp = try
        FluxnetSimulations.get_site_igbp(site_ID)
    catch err
        @warn "Could not read IGBP for $site_ID" exception = err
        "NA"
    end

    return (;
        mean_TA,
        mean_VPD,
        annual_cum_P,
        mean_SW,
        mean_TS_top,
        max_LAI,
        porosity,
        SOC_top = soc_default,
        IGBP = igbp,
    )
end
