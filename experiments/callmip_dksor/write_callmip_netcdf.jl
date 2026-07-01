"""
Write CalLMIP Phase 1a NetCDF output files for DK-Sor.

Phase 1a (Protocol §3.2) is the single-site (DK-Sor) test calibration; the full
21-site set is Phase 1b. Per §7 Note 1, the "Phase1b" filename token is reserved
for Phase 1b — Phase 1a uses the "Phase1a" token and no scenario number (cf. the
registered ME-org outputs CLASSIC.v2_Phase1a_Cal_{Prior,Posterior}).

Reads callmip_diagnostics.jld2 from prior and posterior simulation runs
and writes 2 NetCDF files (ME-org upload Name = ClimaLand.v1_Phase1a_Cal_<stage>):

  ClimaLand.v1_Phase1a_DK-Sor_Cal_Prior.nc
  ClimaLand.v1_Phase1a_DK-Sor_Cal_Posterior.nc

Template spec (variable names/units verified against CLASSIC.v2_Phase1b template):
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

    # hfg and evspsblsoi from the canonical soil-surface flux diagnostics (soilrn,
    # soillhf, soilshf — now exposed for LandModel in src/diagnostics):
    #   hfg        = −soilrn − soillhf − soilshf           [W m⁻²]
    #   evspsblsoi = soillhf / Lv                          [kg m⁻² s⁻¹]
    # NB sign: the `soilrn` diagnostic returns p.soil.R_n, which ClimaLand stores as
    # the *negative* of the net radiation into the soil (src/integrated/land.jl). The
    # soil-surface energy balance is Rnet = LHF + SHF + G, with LHF/SHF positive upward
    # and G (ground heat flux, = hfg) positive downward into the soil, so
    #   hfg = Rnet − LHF − SHF = (−soilrn) − soillhf − soilshf.
    # Verified: this gives an annual-mean hfg ≈ +1.6 W m⁻² (≈ 0, as expected); the
    # un-negated form gave ≈ −64 W m⁻² and did not close the surface budget.
    Lv         = 2.5e6                                     # latent heat of vaporization (J/kg)
    soilrn_    = get_surface("soilrn")
    soillhf_   = get_surface("soillhf")
    soilshf_   = get_surface("soilshf")
    hfg        = fill(FILL_VALUE, N_DAYS)
    evspsblsoi = fill(FILL_VALUE, N_DAYS)
    for t in 1:N_DAYS
        soillhf_[t] != FILL_VALUE && (evspsblsoi[t] = soillhf_[t] / Lv)
        (soilrn_[t] != FILL_VALUE && soillhf_[t] != FILL_VALUE && soilshf_[t] != FILL_VALUE) &&
            (hfg[t] = -soilrn_[t] - soillhf_[t] - soilshf_[t])
    end

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
        unc = JLD2.load(unc_path)   # keys: "member_nee", "member_lhf", "member_shf", "dates"
        # member arrays are (N_members × N_days) indexed by the ENSEMBLE file's own
        # `dates` axis (run_posterior_ensemble.jl) — map via that, not the
        # posterior-mean diagnostic dates, so alignment is independent of the mean run.
        ens_dates = unc["dates"]
        ens_idx   = Dict(d => i for (i, d) in enumerate(ens_dates))
        function ens_std(key, scale, sign)
            out = fill(FILL_VALUE, N_DAYS)
            arr = get(unc, key, nothing)
            isnothing(arr) && return out
            for (t, dt) in enumerate(OUT_START:Day(1):OUT_STOP)
                i = get(ens_idx, dt, nothing)
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
        # Phase 1a = single-site (DK-Sor) test calibration (Protocol §3.2). Per §7
        # Note 1 the "Phase1b" token is ONLY for the full multi-site Phase 1b.
        ds.attrib["Model"]                = "ClimaLand.v1"
        ds.attrib["CalLMIP_Phase"]        = "Phase1a"
        ds.attrib["Calibration_Scenario"] = "Scen1"   # NEE/LHF/SHF (Scenario 1 variables)
        ds.attrib["Calibration_stage"]    = stage
        ds.attrib["Cal_Val"]              = "Calibration"
        ds.attrib["Conventions"]          = "COARDS"
        ds.attrib["site"]                 = "DK-Sor"
        ds.attrib["created"]              = string(today())

        # Dimensions
        defDim(ds, "lon",  1)
        defDim(ds, "lat",  1)
        # Fixed length N_DAYS (=6574). An unlimited dim (Inf) left `v[:,1,1]=data`
        # seeing extent 0 → DimensionMismatch (0,1,1) vs 6574; all data arrays and
        # time_days are exactly N_DAYS, so a fixed dimension is correct.
        defDim(ds, "time", N_DAYS)

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
        # Reference 1996-12-31 so that numeric value 1.0 = 1997-01-01 (matches the
        # CalLMIP template); time_days = 1.0 + (dt - 1997-01-01) gives 1.0…6574.0.
        v_time.attrib["units"]     = "days since 1996-12-31 00:00"  # exact template string
        v_time.attrib["calendar"]  = "standard"

        function defdata(name, data, long_name, units)
            # NCDatasets is column-major; declaring (lon,lat,time) here yields the
            # on-disk/row-major order (time,lat,lon) the CalLMIP template uses.
            v = defVar(ds, name, Float64, ("lon", "lat", "time");
                       fillvalue = FILL_VALUE)
            # Set attributes BEFORE writing data; the fillvalue kwarg above already
            # creates _FillValue, so re-setting it after data exists errors (-122).
            v.attrib["long_name"]   = long_name
            v.attrib["units"]       = units
            v.attrib["coordinates"] = "latitude longitude"
            v[1, 1, :]              = data
        end

        # Units match the CalLMIP template exactly (carbon fluxes are "kgC m-2 s-1").
        defdata("nep",           nep,           "Net ecosystem production",            "kgC m-2 s-1")
        defdata("hfls",          hfls,          "Surface upward latent heat flux",     "W m-2")
        defdata("hfss",          hfss,          "Surface upward sensible heat flux",   "W m-2")
        defdata("gpp",           gpp_nc,        "Gross primary production",            "kgC m-2 s-1")
        defdata("reco",          reco,          "Ecosystem respiration",               "kgC m-2 s-1")
        defdata("tran",          tran,          "Transpiration",                       "kg m-2 s-1")
        defdata("evspsblsoi",    evspsblsoi,    "Bare soil evaporation",               "kg m-2 s-1")
        defdata("hfg",           hfg,           "Downward heat flux at soil surface",  "W m-2")
        defdata("ts",            ts,            "Surface temperature",                 "K")
        defdata("mrso",          mrso,          "Total soil moisture content",         "kg m-2")
        defdata("mrsos",         mrsos,         "Soil moisture in top 10cm",           "kg m-2")
        defdata("lai",           lai_nc,        "Leaf area index",                     "m2 m-2")
        defdata("cSoil",         cSoil,         "Carbon in soil",                      "kg m-2")
        defdata("cLiveBioAbove", cLiveBioAbove, "Carbon in above-ground live biomass", "kg m-2")

        # Posterior uncertainty variables (CalLMIP Protocol Section 6.2)
        if stage == "Posterior"
            defdata("nep_post_unc",  nep_unc,  "Posterior uncertainty of nep",  "kgC m-2 s-1")
            defdata("hfls_post_unc", hfls_unc, "Posterior uncertainty of hfls", "W m-2")
            defdata("hfss_post_unc", hfss_unc, "Posterior uncertainty of hfss", "W m-2")
        end
    end

    @info "Written: $nc_path"

    # Quick validation
    NCDataset(nc_path, "r") do ds
        t = ds["time"].var[:]   # raw numeric (avoid NCDatasets CF DateTime decoding)
        @assert length(t) == N_DAYS "Expected $N_DAYS time steps, got $(length(t))"
        @assert t[1] == 1.0         "Expected time[1] = 1.0, got $(t[1])"
        n_vars = length(keys(ds)) - 5   # subtract time, lat, lon, latitude, longitude
        @info "  $n_vars data variables | $(N_DAYS) time steps | stage=$stage"
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────
# Guard so this file can be `include`d to reuse write_callmip_nc / column_integral
# (e.g. single-year pipeline test) without running the full prior+posterior write.
if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(NC_OUTDIR)

    # Prior: no ensemble uncertainty
    prior_diag = joinpath(SIM_OUTDIR, "prior", "callmip_diagnostics.jld2")
    prior_nc   = joinpath(NC_OUTDIR, "ClimaLand.v1_Phase1a_DK-Sor_Cal_Prior.nc")
    write_callmip_nc(prior_diag, nothing, prior_nc, "Prior")

    # Posterior: ensemble uncertainty propagated from the CES posterior by
    # run_posterior_ensemble.jl (member_nee/lhf/shf + dates).
    post_diag  = joinpath(SIM_OUTDIR, "posterior", "callmip_diagnostics.jld2")
    post_unc   = joinpath(SIM_OUTDIR, "ensemble_diagnostics.jld2")
    post_nc    = joinpath(NC_OUTDIR, "ClimaLand.v1_Phase1a_DK-Sor_Cal_Posterior.nc")
    write_callmip_nc(post_diag, post_unc, post_nc, "Posterior")

    @info "CalLMIP NetCDF output complete."
    @info "Files in: $NC_OUTDIR"
end
@info "Validate with: ncdump -h <file>.nc"
