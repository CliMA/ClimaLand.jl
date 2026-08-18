# Multi-site forcing for CalLMIP Phase 1b on a shared time axis.
#
# The 21 calibration sites' met records cover different periods (intersection
# is only 2009-2010), but the multi-point `TimeVaryingInput` needs one time
# axis for all columns. Sites are therefore PADDED onto a common 30-min UTC
# axis by year-recycling their own record (the protocol §8 spin-up/transient
# convention: cycle the site met), and each padded site carries a `native`
# mask so downstream code can restrict observations to real data.
#
# Padded segments recycle CO2air along with everything else; this only feeds
# state carry-over in the calibration window. The per-site DELIVERABLE runs
# use native windows only and never see padded forcing.

include(joinpath(@__DIR__, "..", "integrated", "fluxnet", "callmip_forcing.jl"))
include(joinpath(@__DIR__, "sites.jl"))

"""
    modis_lai_var(met_nc_path)

The variable name (`"LAI"` or `"LAI_alternative"`) holding the MODIS product in
a PLUMBER2 met file, identified by the `source` attribute; which product sits
under which name varies by site (13/21 calibration sites keep MODIS under
`"LAI"`, 8 under `"LAI_alternative"`).
"""
function modis_lai_var(met_nc_path)
    NCDataset(met_nc_path) do ds
        for v in ("LAI", "LAI_alternative")
            haskey(ds, v) &&
                get(ds[v].attrib, "source", "") == "MODIS" &&
                return v
        end
        error("no MODIS LAI variable in $met_nc_path")
    end
end

"""
    load_site(id; lai_var = :modis)

Read one calibration site's met record using the verified UTC offset from
`CALIBRATION_SITES` and cross-check the file's coordinates, tower height, and
period against the table. LAI defaults to the MODIS product at every site
(decision D7 as revised by Renato 2026-08-18, matching Phase 1a DK-Sor),
resolved per file via [`modis_lai_var`](@ref); pass an explicit variable name
to override.
"""
function load_site(id; lai_var = :modis)
    info = site_info(id)
    path = ClimaLand.Artifacts.callmip_phase1_forcing_path(id)
    lv = lai_var === :modis ? modis_lai_var(path) : lai_var
    site = read_callmip_met(path, info.utc_offset; lai_var = lv)
    isapprox(site.lat, info.lat; atol = 1e-3) &&
    isapprox(site.long, info.long; atol = 1e-3) ||
        error("$id: file coords ($(site.lat), $(site.long)) != site table")
    # reference_height is Float32 in the files; compare with tolerance
    isapprox(site.atmos_h, info.zref; atol = 1e-3) ||
        error("$id: tower height mismatch")
    yr = Dates.year(site.utc_dates[1] + Dates.Hour(info.utc_offset))
    yr == info.met_years[1] || error("$id: met start year $yr != table")
    return resample_to_30min(site)
end

"""
    resample_to_30min(site)

US-MMS reports hourly met (the other 20 calibration sites are half-hourly);
insert linear midpoints so every site shares the 30-min grid required by the
shared-axis matrix `TimeVaryingInput`. Inserted values equal what linear
interpolation of the hourly series would give at those times, so the forcing a
simulation sees is unchanged.
"""
function resample_to_30min(site)
    dt = site.utc_dates[2] - site.utc_dates[1]
    dt == Dates.Minute(30) && return site
    dt == Dates.Minute(60) || error("unsupported met sampling interval $dt")
    n = length(site.utc_dates)
    function mid(v::Vector{Float64})
        w = Vector{Float64}(undef, 2n - 1)
        w[1:2:end] = v
        w[2:2:end] = (v[1:(end - 1)] .+ v[2:end]) ./ 2
        return w
    end
    mid(x) = x # scalars, dates, and `nothing` pass through
    utc = collect(first(site.utc_dates):Dates.Minute(30):last(site.utc_dates))
    return merge(map(mid, site), (; utc_dates = utc))
end

"""
    union_axis(sites; window = nothing)

The 30-min UTC axis spanning all sites' records, optionally cropped to
`window = (start_datetime, stop_datetime)`.
"""
function union_axis(sites; window = nothing)
    t0 = minimum(s -> first(s.utc_dates), sites)
    t1 = maximum(s -> last(s.utc_dates), sites)
    if !isnothing(window)
        t0, t1 = max(t0, window[1]), min(t1, window[2])
    end
    return collect(t0:Dates.Minute(30):t1)
end

"""
    pad_to_axis(site, axis)

Extend `site` (from [`load_site`](@ref)) to the 30-min UTC `axis`. Times inside
the native record keep their values; times outside reuse the same
month/day/time from an INTERIOR native year (first/last years are excluded so
UTC offsets can't make the source ragged), cycling through those years.
Feb 29 maps to Feb 28 when the source year is not a leap year. The result
carries `native::BitVector` marking real (non-recycled) samples.
"""
function pad_to_axis(site, axis)
    idx = Dict(t => i for (i, t) in enumerate(site.utc_dates))
    yr_lo = Dates.year(first(site.utc_dates)) + 1
    yr_hi = Dates.year(last(site.utc_dates)) - 1
    yr_hi >= yr_lo || error("record too short for interior-year recycling")
    n_years = yr_hi - yr_lo + 1
    srcidx = Vector{Int}(undef, length(axis))
    native = falses(length(axis))
    for (k, t) in enumerate(axis)
        i = get(idx, t, 0)
        if i > 0
            srcidx[k] = i
            native[k] = true
        else
            ysrc = yr_lo + mod(Dates.year(t) - yr_lo, n_years)
            d = Dates.day(t)
            m = Dates.month(t)
            m == 2 && d == 29 && !Dates.isleapyear(ysrc) && (d = 28)
            srcidx[k] =
                idx[DateTime(ysrc, m, d, Dates.hour(t), Dates.minute(t))]
        end
    end
    remap(v::Vector) = v[srcidx]
    remap(x) = x # scalars and `nothing` pass through
    return merge(map(remap, site), (; utc_dates = collect(axis), native))
end

"""
    load_calibration_sites(site_ids = CALIBRATION_SITE_IDS;
                           window = nothing, lai_var = :modis)

Load and pad the requested sites onto their common (optionally cropped) axis,
ready for `prescribed_forcing_callmip` / `build_callmip_simulation`. Returns
the padded site vector; site order defines the ensemble column order.
"""
function load_calibration_sites(
    site_ids = CALIBRATION_SITE_IDS;
    window = nothing,
    lai_var = :modis,
)
    sites = [load_site(id; lai_var) for id in site_ids]
    axis = union_axis(sites; window)
    return [pad_to_axis(s, axis) for s in sites]
end
