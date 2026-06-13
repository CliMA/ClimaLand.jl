"""
Write CalLMIP Phase 1b NetCDF output files for DK-Sor.

Reads callmip_diagnostics.jld2 from prior and posterior simulation runs
and writes 2 NetCDF files matching the official CalLMIP template exactly:

  ClimaLand.v1_Phase1b_Scen1_DK-Sor_Cal_Prior.nc
  ClimaLand.v1_Phase1b_Scen1_DK-Sor_Cal_Posterior.nc

Template spec (verified against CLASSIC.v2_Phase1b_Scen1_DK-Sor_Cal_Prior.nc):
  - Time axis: "days since 1997-01-01 00:00:00", standard calendar
    Day 1 = 1997-01-01, day 6574 = 2014-12-31
  - Dims: time, lat(=1), lon(=1)
  - _FillValue = 1e38 for all data variables
  - Conventions = "COARDS"

Variable mapping (model → CalLMIP):
  nee     → nep        : -(nee × 12e-3)           [kgC m⁻² s⁻¹]  sign flip: nep = -NEE
  lhf     → hfls       : lhf                        [W m⁻²]
  shf     → hfss       : shf                        [W m⁻²]
  gpp     → gpp        : gpp × 12e-3               [kgC m⁻² s⁻¹]
  er      → reco       : er × 12e-3                [kgC m⁻² s⁻¹]
  trans   → tran       : trans                      [kg m⁻² s⁻¹]
  *       → evspsblsoi : fill value (soillhf not in LandModel diagnostics)
  *       → hfg        : fill value (soilrn/soillhf/soilshf not in LandModel diagnostics)
  ct      → ts         : ct                         [K]
  swc     → mrso       : ∫swc·dz·1000 (full col)  [kg m⁻²]
  swc     → mrsos      : ∫swc·dz·1000 (z>-0.1m)   [kg m⁻²]
  lai     → lai        : lai                        [m² m⁻²]
  soc     → cSoil      : ∫soc·dz                   [kg m⁻²]
  cveg    → cLiveBioAbove : cveg                    [kg m⁻²]

Posterior file also includes uncertainty variables (CalLMIP Protocol Section 6.2):
  nep_post_unc, hfls_post_unc, hfss_post_unc  [same units as parent variable]

Usage:
  julia --project=.buildkite \
        experiments/callmip_dksor/write_callmip_netcdf.jl
"""

using NCDatasets
using Dates
using JLD2
using Statistics

import ClimaLand

const CLIMALAND_DIR = pkgdir(ClimaLand)
const SIM_OUTDIR    = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_callmip_sims")
const NC_OUTDIR     = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_netcdf")

# CalLMIP time axis: "days since 1997-01-01 00:00:00" (first year of FLUXNET data)
# Day 1 = 1997-01-01, day 6574 = 2014-12-31
const TIME_REF  = Date(1997, 1, 1)
const OUT_START = Date(1997, 1, 1)
const OUT_STOP  = Date(2014, 12, 31)
const N_DAYS    = (OUT_STOP - OUT_START).value + 1   # 6574

# Site lat/lon
const SITE_LAT = 55.49
const SITE_LON = 11.64

# Fill value matching CalLMIP template (double precision)
const FILL_VALUE = Float64(1.0e38)

# ── Column integration helper ─────────────────────────────────────────────────

"""
    column_integral(profile, z_centers, z_max) -> Vector{Float64}

Integrate a (n_z × n_days) soil profile over depth using layer thicknesses
derived from z_centers (m, negative downward). Optionally restrict to
layers above z_max (e.g., z_max = -0.1 for top 10 cm).
"""
function column_integral(profile::Matrix{Float64}, z_centers::Vector{Float64};
                          z_max::Float64 = -Inf)
    n_z, n_days = size(profile)
    idx  = sortperm(z_centers; rev = true)   # shallowest first
    z_s  = z_centers[idx]
    dz   = diff([-0.0; z_s])
    dz   = abs.(dz)

    result = zeros(Float64, n_days)
    for (i, (zi, dzi)) in enumerate(zip(z_s, dz))
        zi >= z_max || continue
        result .+= profile[idx[i], :] .* dzi
    end
    return result
end

# ── Write one NetCDF file ─────────────────────────────────────────────────────

"""
    write_callmip_nc(diag_path, unc_path, nc_path, stage)

Read JLD2 diagnostics and write a CalLMIP-compliant NetCDF file.
`stage` is "Prior" or "Posterior".
`unc_path` is either the path to a posterior ensemble diagnostics JLD2
(for posterior uncertainty variables) or `nothing` (for prior).
"""
function write_callmip_nc(diag_path::String, unc_path, nc_path::String, stage::String)
    isfile(diag_path) || error("Diagnostics not found: $diag_path")

    d = JLD2.load(diag_path)
    dates        = d["dates"]           # Vector{Date}
    surface_data = d["surface_data"]    # Dict{String, Vector{Float64}}
    column_data  = d["column_data"]     # Dict{String, Matrix{Float64}}
    z_soil       = d["z_soil"]          # Vector{Float64}, meters (negative down)

    # Time axis: days since 1997-01-01 (day 1 = 1997-01-01)
    time_days = Float64[1.0 + (dt - TIME_REF).value for dt in OUT_START:Day(1):OUT_STOP]
    @assert length(time_days) == N_DAYS
    @assert time_days[1] == 1.0       # 1997-01-01 = day 1
    @assert time_days[end] == N_DAYS  # 2014-12-31 = day 6574

    # Map date → index in model output
    date_to_idx = Dict(dt => i for (i, dt) in enumerate(dates))

    function get_surface(var; scale = 1.0, offset = 0.0, sign = 1.0)
        out = fill(FILL_VALUE, N_DAYS)
        src = get(surface_data, var, Float64[])
        isempty(src) && return out
        for (t, dt) in enumerate(OUT_START:Day(1):OUT_STOP)
            i = get(date_to_idx, dt, nothing)
            isnothing(i) && continue
            v = src[i]
            isfinite(v) || continue
            out[t] = sign * v * scale + offset
        end
        return out
    end

    function get_column_integral(var; scale = 1.0, z_max = -Inf)
        out = fill(FILL_VALUE, N_DAYS)
        mat = get(column_data, var, Matrix{Float64}(undef, 0, 0))
        (isempty(mat) || isempty(z_soil)) && return out
        col = column_integral(mat, z_soil; z_max) .* scale
        for (t, dt) in enumerate(OUT_START:Day(1):OUT_STOP)
            i = get(date_to_idx, dt, nothing)
            isnothing(i) && continue
            v = col[i]
            isfinite(v) && (out[t] = v)
        end
        return out
    end

    # ── Variable computations ─────────────────────────────────────────────────
    nep          = get_surface("nee"; scale = 12e-3, sign = -1.0)
    hfls         = get_surface("lhf")
    hfss         = get_surface("shf")
    gpp_nc       = get_surface("gpp"; scale = 12e-3)
    reco         = get_surface("er";  scale = 12e-3)
    tran         = get_surface("trans")
    ts           = get_surface("ct")
    lai_nc       = get_surface("lai")
    cLiveBioAbove = get_surface("cveg")

    # hfg and evspsblsoi require soil flux diagnostics (soilrn, soillhf, soilshf)
    # which are not available in get_possible_diagnostics(LandModel) — set to fill.
    hfg        = fill(FILL_VALUE, N_DAYS)
    evspsblsoi = fill(FILL_VALUE, N_DAYS)

    mrso  = get_column_integral("swc"; scale = 1000.0)
    mrsos = get_column_integral("swc"; scale = 1000.0, z_max = -0.1)
    cSoil = get_column_integral("soc")

    # ── Posterior uncertainties (CalLMIP Section 6.2): std across ensemble ────
    # Loaded from ensemble diagnostics saved by emulate_sample.jl.
    # Each ensemble member produces its own diagnostics; we compute std over members.
    nep_unc  = fill(FILL_VALUE, N_DAYS)
    hfls_unc = fill(FILL_VALUE, N_DAYS)
    hfss_unc = fill(FILL_VALUE, N_DAYS)
    if !isnothing(unc_path) && isfile(unc_path)
        unc = JLD2.load(unc_path)   # keys: "member_nee", "member_lhf", "member_shf"
        # member arrays: (N_members × N_days_model) after alignment
        function ens_std(key, scale, sign)
            out = fill(FILL_VALUE, N_DAYS)
            arr = get(unc, key, nothing)
            isnothing(arr) && return out
            for (t, dt) in enumerate(OUT_START:Day(1):OUT_STOP)
                i = get(date_to_idx, dt, nothing)
                isnothing(i) && continue
                col = arr[:, i]
                finite_col = filter(isfinite, col)
                isempty(finite_col) && continue
                out[t] = std(finite_col) * abs(scale)
            end
            return out
        end
        nep_unc  = ens_std("member_nee", 12e-3, -1.0)
        hfls_unc = ens_std("member_lhf", 1.0,   1.0)
        hfss_unc = ens_std("member_shf", 1.0,   1.0)
        @info "Posterior uncertainties computed from $(size(unc["member_nee"], 1)) ensemble members"
    elseif stage == "Posterior"
        @warn "No ensemble diagnostics found at $unc_path — _post_unc variables will be all fill"
    end

    # ── Write NetCDF ──────────────────────────────────────────────────────────
    NCDataset(nc_path, "c") do ds
        # Global attributes
        ds.attrib["Model"]                = "ClimaLand"
        ds.attrib["Model_version"]        = "v1"
        ds.attrib["CalLMIP_Phase"]        = "Phase1b"
        ds.attrib["Calibration_Scenario"] = "Scen1"
        ds.attrib["Calibration_stage"]    = stage
        ds.attrib["Cal_Val"]              = "Calibration"
        ds.attrib["Conventions"]          = "COARDS"
        ds.attrib["site"]                 = "DK-Sor"
        ds.attrib["created"]              = string(today())

        # Dimensions
        defDim(ds, "lon",  1)
        defDim(ds, "lat",  1)
        defDim(ds, "time", Inf)   # unlimited

        # Scalar lat/lon (matching template structure)
        v_lat_scalar = defVar(ds, "lat", Float64, ())
        v_lat_scalar[] = SITE_LAT

        v_lon_scalar = defVar(ds, "lon", Float64, ())
        v_lon_scalar[] = SITE_LON

        v_latitude = defVar(ds, "latitude", Float64, ("lat",))
        v_latitude[:] = [SITE_LAT]
        v_latitude.attrib["standard_name"] = "Latitude"
        v_latitude.attrib["long_name"]     = "latitude"
        v_latitude.attrib["units"]         = "degrees_north"

        v_longitude = defVar(ds, "longitude", Float64, ("lon",))
        v_longitude[:] = [SITE_LON]
        v_longitude.attrib["standard_name"] = "Longitude"
        v_longitude.attrib["long_name"]     = "longitude"
        v_longitude.attrib["units"]         = "degrees_east"

        v_time = defVar(ds, "time", Float64, ("time",))
        v_time[:] = time_days
        v_time.attrib["long_name"] = "time"
        v_time.attrib["units"]     = "days since 1997-01-01 00:00:00"
        v_time.attrib["calendar"]  = "standard"

        function defdata(name, data, long_name, units)
            v = defVar(ds, name, Float64, ("time", "lat", "lon");
                       fillvalue = FILL_VALUE)
            v[:, 1, 1]              = data
            v.attrib["long_name"]   = long_name
            v.attrib["units"]       = units
            v.attrib["_FillValue"]  = FILL_VALUE
            v.attrib["coordinates"] = "latitude longitude"
        end

        defdata("nep",           nep,           "Net ecosystem production",            "kg m-2 s-1")
        defdata("hfls",          hfls,          "Surface upward latent heat flux",     "W m-2")
        defdata("hfss",          hfss,          "Surface upward sensible heat flux",   "W m-2")
        defdata("gpp",           gpp_nc,        "Gross primary production",            "kg m-2 s-1")
        defdata("reco",          reco,          "Ecosystem respiration",               "kg m-2 s-1")
        defdata("tran",          tran,          "Transpiration",                       "kg m-2 s-1")
        defdata("evspsblsoi",    evspsblsoi,    "Bare soil evaporation",               "kg m-2 s-1")
        defdata("hfg",           hfg,           "Downward heat flux at soil surface",  "W m-2")
        defdata("ts",            ts,            "Surface temperature",                 "K")
        defdata("mrso",          mrso,          "Total soil moisture content",         "kg m-2")
        defdata("mrsos",         mrsos,         "Soil moisture in top 10cm",           "kg m-2")
        defdata("lai",           lai_nc,        "Leaf area index",                     "1")
        defdata("cSoil",         cSoil,         "Carbon in soil",                      "kg m-2")
        defdata("cLiveBioAbove", cLiveBioAbove, "Carbon in above-ground live biomass", "kg m-2")

        # Posterior uncertainty variables (CalLMIP Protocol Section 6.2)
        if stage == "Posterior"
            defdata("nep_post_unc",  nep_unc,  "Posterior uncertainty of nep",  "kg m-2 s-1")
            defdata("hfls_post_unc", hfls_unc, "Posterior uncertainty of hfls", "W m-2")
            defdata("hfss_post_unc", hfss_unc, "Posterior uncertainty of hfss", "W m-2")
        end
    end

    @info "Written: $nc_path"

    # Quick validation
    NCDataset(nc_path, "r") do ds
        t = ds["time"][:]
        @assert length(t) == N_DAYS "Expected $N_DAYS time steps, got $(length(t))"
        @assert t[1] == 1.0         "Expected time[1] = 1.0, got $(t[1])"
        n_vars = length(keys(ds)) - 5   # subtract time, lat, lon, latitude, longitude
        @info "  $n_vars data variables | $(N_DAYS) time steps | stage=$stage"
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────
mkpath(NC_OUTDIR)

# Prior: no ensemble uncertainty
prior_diag = joinpath(SIM_OUTDIR, "prior", "callmip_diagnostics.jld2")
prior_nc   = joinpath(NC_OUTDIR, "ClimaLand.v1_Phase1b_Scen1_DK-Sor_Cal_Prior.nc")
write_callmip_nc(prior_diag, nothing, prior_nc, "Prior")

# Posterior: ensemble uncertainty from emulate_sample.jl output
post_diag  = joinpath(SIM_OUTDIR, "posterior", "callmip_diagnostics.jld2")
post_unc   = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_ces/posterior_ensemble_diags.jld2")
post_nc    = joinpath(NC_OUTDIR, "ClimaLand.v1_Phase1b_Scen1_DK-Sor_Cal_Posterior.nc")
write_callmip_nc(post_diag, post_unc, post_nc, "Posterior")

@info "CalLMIP NetCDF output complete."
@info "Files in: $NC_OUTDIR"
@info "Validate with: ncdump -h <file>.nc"
