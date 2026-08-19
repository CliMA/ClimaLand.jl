"""
Parametrized CalLMIP Phase 1b NetCDF writer (per site, Cal or Val).

Port of experiments/callmip_dksor/write_callmip_netcdf.jl with per-site
metadata: time axis "days since <y0-1>-12-31 00:00" (value 1.0 = Jan 1 of the
site's met start year), site lat/lon, Phase1b/Scen1 global attributes.
Variable mapping, units, and the hfg sign convention are identical to Phase 1a
(verified there against the CLASSIC.v2 template).

Include this file and call `write_callmip_nc_phase1b(...)`.
"""

using NCDatasets
using Dates
using JLD2
using Statistics

const FILL_VALUE_1B = Float64(1.0e38)

function column_integral_1b(
    profile::Matrix{Float64},
    z_centers::Vector{Float64};
    z_max::Float64 = -Inf,
)
    idx = sortperm(z_centers; rev = true)
    z_s = z_centers[idx]
    dz = abs.(diff([-0.0; z_s]))
    result = zeros(Float64, size(profile, 2))
    for (i, (zi, dzi)) in enumerate(zip(z_s, dz))
        zi >= z_max || continue
        result .+= profile[idx[i], :] .* dzi
    end
    return result
end

"""
    write_callmip_nc_phase1b(diag_path, nc_path;
                             site_id, lat, long, met_years,
                             stage, cal_val, unc_path = nothing)

Read a per-site diagnostics JLD2 (dates, surface_data, column_data, z_soil)
and write one CalLMIP Phase 1b NetCDF file. `stage` is "Prior"/"Posterior",
`cal_val` is "Calibration"/"Validation".
"""
function write_callmip_nc_phase1b(
    diag_path::String,
    nc_path::String;
    site_id::String,
    lat::Float64,
    long::Float64,
    met_years::Tuple{Int, Int},
    stage::String,
    cal_val::String,
    unc_path = nothing,
)
    isfile(diag_path) || error("Diagnostics not found: $diag_path")
    d = JLD2.load(diag_path)
    dates = d["dates"]
    surface_data = d["surface_data"]
    column_data = d["column_data"]
    z_soil = d["z_soil"]

    y0, y1 = met_years
    out_start, out_stop = Date(y0, 1, 1), Date(y1, 12, 31)
    n_days = (out_stop - out_start).value + 1
    time_days = collect(1.0:Float64(n_days))
    out_dates = out_start:Day(1):out_stop
    date_to_idx = Dict(dt => i for (i, dt) in enumerate(dates))

    function get_surface(var; scale = 1.0, sign = 1.0)
        out = fill(FILL_VALUE_1B, n_days)
        src = get(surface_data, var, Float64[])
        isempty(src) && return out
        for (t, dt) in enumerate(out_dates)
            i = get(date_to_idx, dt, nothing)
            isnothing(i) && continue
            v = src[i]
            isfinite(v) && (out[t] = sign * v * scale)
        end
        return out
    end
    function get_column_integral(var; scale = 1.0, z_max = -Inf)
        out = fill(FILL_VALUE_1B, n_days)
        mat = get(column_data, var, Matrix{Float64}(undef, 0, 0))
        (isempty(mat) || isempty(z_soil)) && return out
        col = column_integral_1b(mat, z_soil; z_max) .* scale
        for (t, dt) in enumerate(out_dates)
            i = get(date_to_idx, dt, nothing)
            isnothing(i) && continue
            isfinite(col[i]) && (out[t] = col[i])
        end
        return out
    end

    nep = get_surface("nee"; scale = 12e-3, sign = -1.0)
    hfls = get_surface("lhf")
    hfss = get_surface("shf")
    gpp_nc = get_surface("gpp"; scale = 12e-3)
    reco = get_surface("er"; scale = 12e-3)
    tran = get_surface("trans")
    ts = get_surface("ct")
    lai_nc = get_surface("lai")
    cLiveBioAbove = get_surface("cveg")

    # hfg = (−soilrn) − soillhf − soilshf; evspsblsoi = soillhf/Lv
    # (sign convention verified in Phase 1a: annual-mean hfg ≈ 0)
    Lv = 2.5e6
    soilrn_ = get_surface("soilrn")
    soillhf_ = get_surface("soillhf")
    soilshf_ = get_surface("soilshf")
    hfg = fill(FILL_VALUE_1B, n_days)
    evspsblsoi = fill(FILL_VALUE_1B, n_days)
    for t in 1:n_days
        soillhf_[t] != FILL_VALUE_1B && (evspsblsoi[t] = soillhf_[t] / Lv)
        (
            soilrn_[t] != FILL_VALUE_1B &&
            soillhf_[t] != FILL_VALUE_1B &&
            soilshf_[t] != FILL_VALUE_1B
        ) && (hfg[t] = -soilrn_[t] - soillhf_[t] - soilshf_[t])
    end

    mrso = get_column_integral("swc"; scale = 1000.0)
    mrsos = get_column_integral("swc"; scale = 1000.0, z_max = -0.1)
    cSoil = get_column_integral("soc")

    # posterior uncertainties from an EKI-member diagnostics file (optional)
    nep_unc = fill(FILL_VALUE_1B, n_days)
    hfls_unc = fill(FILL_VALUE_1B, n_days)
    hfss_unc = fill(FILL_VALUE_1B, n_days)
    if !isnothing(unc_path) && isfile(unc_path)
        unc = JLD2.load(unc_path)
        ens_idx = Dict(dt => i for (i, dt) in enumerate(unc["dates"]))
        function ens_std(key, scale)
            out = fill(FILL_VALUE_1B, n_days)
            arr = get(unc, key, nothing)
            isnothing(arr) && return out
            for (t, dt) in enumerate(out_dates)
                i = get(ens_idx, dt, nothing)
                isnothing(i) && continue
                fc = filter(isfinite, arr[:, i])
                isempty(fc) || (out[t] = std(fc) * abs(scale))
            end
            return out
        end
        nep_unc = ens_std("member_nee", 12e-3)
        hfls_unc = ens_std("member_lhf", 1.0)
        hfss_unc = ens_std("member_shf", 1.0)
    elseif stage == "Posterior"
        @warn "$site_id: no ensemble diagnostics — _post_unc variables all fill"
    end

    NCDataset(nc_path, "c") do ds
        ds.attrib["Model"] = "ClimaLand.v1"
        ds.attrib["CalLMIP_Phase"] = "Phase1b"
        ds.attrib["Calibration_Scenario"] = "Scen1"
        ds.attrib["Calibration_stage"] = stage
        ds.attrib["Cal_Val"] = cal_val
        ds.attrib["Conventions"] = "COARDS"
        ds.attrib["site"] = site_id
        ds.attrib["created"] = string(today())

        defDim(ds, "lon", 1)
        defDim(ds, "lat", 1)
        defDim(ds, "time", n_days)

        v_lat = defVar(ds, "lat", Float64, ())
        v_lat[] = lat
        v_lon = defVar(ds, "lon", Float64, ())
        v_lon[] = long
        v_latitude = defVar(ds, "latitude", Float64, ("lat",))
        v_latitude[:] = [lat]
        v_latitude.attrib["standard_name"] = "Latitude"
        v_latitude.attrib["long_name"] = "latitude"
        v_latitude.attrib["units"] = "degrees_north"
        v_longitude = defVar(ds, "longitude", Float64, ("lon",))
        v_longitude[:] = [long]
        v_longitude.attrib["standard_name"] = "Longitude"
        v_longitude.attrib["long_name"] = "longitude"
        v_longitude.attrib["units"] = "degrees_east"

        v_time = defVar(ds, "time", Float64, ("time",))
        v_time[:] = time_days
        v_time.attrib["long_name"] = "time"
        v_time.attrib["units"] = "days since $(y0 - 1)-12-31 00:00"
        v_time.attrib["calendar"] = "standard"

        function defdata(name, data, long_name, units)
            v = defVar(
                ds,
                name,
                Float64,
                ("lon", "lat", "time");
                fillvalue = FILL_VALUE_1B,
            )
            v.attrib["long_name"] = long_name
            v.attrib["units"] = units
            v.attrib["coordinates"] = "latitude longitude"
            v[1, 1, :] = data
        end

        defdata("nep", nep, "Net ecosystem production", "kgC m-2 s-1")
        defdata("hfls", hfls, "Surface upward latent heat flux", "W m-2")
        defdata("hfss", hfss, "Surface upward sensible heat flux", "W m-2")
        defdata("gpp", gpp_nc, "Gross primary production", "kgC m-2 s-1")
        defdata("reco", reco, "Ecosystem respiration", "kgC m-2 s-1")
        defdata("tran", tran, "Transpiration", "kg m-2 s-1")
        defdata("evspsblsoi", evspsblsoi, "Bare soil evaporation", "kg m-2 s-1")
        defdata("hfg", hfg, "Downward heat flux at soil surface", "W m-2")
        defdata("ts", ts, "Surface temperature", "K")
        defdata("mrso", mrso, "Total soil moisture content", "kg m-2")
        defdata("mrsos", mrsos, "Soil moisture in top 10cm", "kg m-2")
        defdata("lai", lai_nc, "Leaf area index", "m2 m-2")
        defdata("cSoil", cSoil, "Carbon in soil", "kg m-2")
        defdata(
            "cLiveBioAbove",
            cLiveBioAbove,
            "Carbon in above-ground live biomass",
            "kg m-2",
        )
        if stage == "Posterior"
            defdata(
                "nep_post_unc",
                nep_unc,
                "Posterior uncertainty of nep",
                "kgC m-2 s-1",
            )
            defdata(
                "hfls_post_unc",
                hfls_unc,
                "Posterior uncertainty of hfls",
                "W m-2",
            )
            defdata(
                "hfss_post_unc",
                hfss_unc,
                "Posterior uncertainty of hfss",
                "W m-2",
            )
        end
    end
    return nc_path
end
