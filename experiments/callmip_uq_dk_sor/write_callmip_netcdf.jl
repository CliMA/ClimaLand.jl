"""
Write CalLMIP Phase 1 compliant NetCDF output for the DK-Sor site.

Reads:
  • output_callmip_sims/iteration_000/member_001/callmip_diagnostics.jld2  (prior)
  • output_callmip_sims/iteration_000/member_002/callmip_diagnostics.jld2  (posterior)
  • output_posterior_analysis/posterior_uncertainty_bands.jld2              (optional,
      for NEE / Qle / Qh uncertainty bands in posterior Cal/Val files)

Writes 4 NetCDF files:
  ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Prior.nc
  ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Posterior.nc
  ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Prior.nc
  ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Posterior.nc

Variable list (ALMA names, daily frequency):
  NEE, Qle, Qh, GPP, Reco, TVeg, ESoil, Qg,
  AvgSurfT, SoilMoist, LAI, TotAbovBioMass, TotSoilCarb
  + NEE_p05 / NEE_p50 / NEE_p95 … (posterior files only, where available)

Units conventions
-----------------
  Carbon fluxes : mol CO₂ m⁻² s⁻¹  →  kg C m⁻² s⁻¹  (× 12e-3)
  ESoil         : W m⁻²             →  kg m⁻² s⁻¹    (÷ 2.5e6)
  SoilMoist     : m³ m⁻³ × dz × ρ  →  kg m⁻²         (ρ = 1000 kg m⁻³)
  TotSoilCarb   : kg C m⁻³ × dz    →  kg C m⁻²
  All others pass through without conversion.

Temporal split
--------------
  Cal : first CAL_END_YEAR years of simulation  (calibration period)
  Val : remaining years                          (temporal validation)
  Modify CAL_END_YEAR below if a different split is needed.

Usage
-----
    julia --project=experiments/callmip_uq_dk_sor \\
          experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl
"""

using Dates
import JLD2
import NCDatasets

# ── Configuration ──────────────────────────────────────────────────────────────
const MODEL_NAME  = "ClimaLand"
const MODEL_VER   = "CalLMIP1.0"
const SITE_ID     = "DK-Sor"
const EXPT_NO     = "1"
const CAL_END_YEAR = 2012          # last year of calibration period (inclusive)
                                    # → validation = years > CAL_END_YEAR

import ClimaLand
const climaland_dir         = abspath(joinpath(@__DIR__, "..", ".."))
const cal_dir               = joinpath(climaland_dir, "experiments", "calibrate_dk_sor")   # calibration artifacts
const exp_dir               = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")   # UQ/CalLMIP outputs
const callmip_sim_dir       = joinpath(exp_dir, "output_callmip_sims")
const posterior_analysis_dir = joinpath(exp_dir, "output_posterior_analysis")
const nc_output_dir         = joinpath(exp_dir, "callmip_output")
isdir(nc_output_dir) || mkpath(nc_output_dir)

# ── Load diagnostics ──────────────────────────────────────────────────────────
function load_diagnostics(member::Int)
    path = joinpath(callmip_sim_dir,
                    "iteration_000",
                    "member_$(lpad(member, 3, '0'))",
                    "callmip_diagnostics.jld2")
    isfile(path) || error("Diagnostics file not found: $path")
    JLD2.load(path)
end

@info "Loading prior diagnostics (member 001)…"
prior_d = load_diagnostics(1)
@info "Loading posterior diagnostics (member 002)…"
post_d  = load_diagnostics(2)

# Dates vector (daily, from 2004-01-01 after 2003 spinup)
dates = prior_d["dates"]           # Vector{Date}
n_days = length(dates)

@info "Model output period: $(first(dates)) – $(last(dates))  ($n_days days)"

# ── Temporal split masks ──────────────────────────────────────────────────────
cal_mask = year.(dates) .<= CAL_END_YEAR
val_mask = .!cal_mask
n_cal    = sum(cal_mask)
n_val    = sum(val_mask)
@info "Cal period: $(first(dates[cal_mask])) – $(last(dates[cal_mask]))  ($n_cal days)"
@info "Val period: $(first(dates[val_mask])) – $(last(dates[val_mask]))  ($n_val days)"

# ── Load posterior uncertainty bands (optional) ───────────────────────────────
has_uncertainty = false
unc_d = nothing
try
    unc_path = joinpath(posterior_analysis_dir, "posterior_uncertainty_bands.jld2")
    isfile(unc_path) || error("not found")
    unc_d = JLD2.load(unc_path)
    has_uncertainty = true
    @info "Loaded posterior uncertainty bands."
catch e
    @warn "Posterior uncertainty bands not available ($e). " *
          "Run analyze_posterior_ensemble.jl to add them to the posterior files."
end

# ── Unit conversion helpers ───────────────────────────────────────────────────
const mol_CO2_to_kgC = 12e-3             # mol CO₂ m⁻² s⁻¹ → kg C m⁻² s⁻¹
const Lv_water       = 2.5e6            # J kg⁻¹ (latent heat of vaporisation)
const rho_water      = 1000.0           # kg m⁻³

"""Integrate a column diagnostic over soil depth → kg m⁻² (or kg C m⁻²).

`col_data` must be a Matrix{Float64} with rows = z levels (bottom to surface)
and columns = time steps.  `z_soil` is the z-coordinate vector (m, negative).
`scale` is applied after integration (e.g. rho_water = 1000 for SoilMoist).
"""
function column_integral(col_data::Matrix{Float64}, z_soil::Vector{Float64};
                          scale::Float64 = 1.0)
    n_z, n_t = size(col_data)
    n_z == length(z_soil) || error("col_data rows $(n_z) ≠ z_soil length $(length(z_soil))")
    # Layer thicknesses via midpoint rule (cell centre spacing)
    dz = zeros(n_z)
    if n_z == 1
        dz[1] = abs(z_soil[1])
    else
        for i in 1:n_z
            lo = i == 1   ? 2*z_soil[1] - z_soil[2] : (z_soil[i-1] + z_soil[i]) / 2
            hi = i == n_z ? 0.0                      : (z_soil[i]   + z_soil[i+1]) / 2
            dz[i] = hi - lo
        end
    end
    result = zeros(n_t)
    for t in 1:n_t
        for i in 1:n_z
            result[t] += col_data[i, t] * dz[i] * scale
        end
    end
    return result
end

# ── Extract and convert all variables ────────────────────────────────────────
"""
Return a Dict with all ALMA-named fields (Vector{Float64}, length = n_days).
Also includes `_p05`, `_p50`, `_p95` variants for obs-optimized variables if
uncertainty bands were loaded (posterior only).
"""
function build_variable_dict(d::Dict, add_uncertainty::Bool = false)
    out = Dict{String, Vector{Float64}}()

    # Helper to safely get a surface diagnostic or fill with NaN
    function get_surf(key)
        haskey(d, key) ? Float64.(d[key]) : fill(NaN, n_days)
    end

    z_soil = haskey(d, "z_soil") ? Float64.(d["z_soil"]) : Float64[]

    # ── Carbon fluxes (mol CO₂ m⁻² s⁻¹ → kg C m⁻² s⁻¹) ──────────────────
    out["NEE"]  = get_surf("nee")  .* mol_CO2_to_kgC
    out["GPP"]  = get_surf("gpp")  .* mol_CO2_to_kgC
    out["Reco"] = get_surf("er")   .* mol_CO2_to_kgC

    # ── Heat fluxes (W m⁻², pass-through) ─────────────────────────────────
    out["Qle"] = get_surf("lhf")
    out["Qh"]  = get_surf("shf")

    # ── Water fluxes ───────────────────────────────────────────────────────
    # TVeg: transpiration already in kg m⁻² s⁻¹ (computed as vapor_flux × 1000)
    out["TVeg"]  = get_surf("trans")
    # ESoil: soil latent heat → kg m⁻² s⁻¹
    out["ESoil"] = get_surf("soillhf") ./ Lv_water

    # ── Ground heat flux (Qg = Rn_soil − LHF_soil − SHF_soil) ─────────────
    soilrn  = get_surf("soilrn")
    soillhf = get_surf("soillhf")
    soilshf = get_surf("soilshf")
    out["Qg"] = soilrn .- soillhf .- soilshf

    # ── Surface temperature (K) ────────────────────────────────────────────
    # Primary: canopy temperature (ct); fallback: top soil layer
    ct    = get_surf("ct")
    tsoil = if haskey(d, "tsoil") && d["tsoil"] isa Matrix
        # Top soil layer = last row (z closest to 0)
        Float64.(d["tsoil"][end, :])
    else
        fill(NaN, n_days)
    end
    avgt = copy(ct)
    for i in 1:n_days
        if isnan(avgt[i]) || avgt[i] < 150.0   # sanity: below 150 K is unphysical
            avgt[i] = tsoil[i]
        end
    end
    out["AvgSurfT"] = avgt

    # ── Total column soil moisture (kg m⁻²) ───────────────────────────────
    if haskey(d, "swc") && d["swc"] isa Matrix && !isempty(z_soil)
        out["SoilMoist"] = column_integral(
            Float64.(d["swc"]), z_soil; scale = rho_water)
    else
        out["SoilMoist"] = fill(NaN, n_days)
        @warn "swc column diagnostic not available; SoilMoist set to NaN."
    end

    # ── LAI (m² m⁻²) ───────────────────────────────────────────────────────
    out["LAI"] = get_surf("lai")

    # ── Total above-ground biomass carbon (kg C m⁻²) ──────────────────────
    out["TotAbovBioMass"] = get_surf("cveg")

    # ── Total soil organic carbon (kg C m⁻²) ──────────────────────────────
    if haskey(d, "soc") && d["soc"] isa Matrix && !isempty(z_soil)
        out["TotSoilCarb"] = column_integral(Float64.(d["soc"]), z_soil)
    else
        out["TotSoilCarb"] = fill(NaN, n_days)
        @warn "soc column diagnostic not available; TotSoilCarb set to NaN."
    end

    # ── Posterior uncertainty percentiles (NEE, Qle, Qh) ──────────────────
    if add_uncertainty && has_uncertainty
        for var in ("NEE", "Qle", "Qh")
            for pct in ("p05", "p50", "p95")
                key_out = "$(var)_$(pct)"
                key_in  = lowercase(var) * "_$(pct)"  # e.g. "nee_p05"
                if haskey(unc_d, key_in)
                    raw = Float64.(unc_d[key_in])
                    # Apply same conversion as point estimate
                    if var == "NEE"
                        out[key_out] = raw .* mol_CO2_to_kgC
                    else
                        out[key_out] = raw
                    end
                end
            end
        end
    end

    return out
end

@info "Building variable dictionaries…"
prior_vars = build_variable_dict(prior_d,  false)
post_vars  = build_variable_dict(post_d,   true)

# ── ALMA metadata ─────────────────────────────────────────────────────────────
const ALMA_META = Dict{String, NamedTuple}(
    "NEE"          => (long_name = "Net Ecosystem Exchange",
                       units     = "kg C m-2 s-1",
                       standard_name = "surface_net_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_all_land_processes"),
    "Qle"          => (long_name = "Latent Heat Flux",
                       units     = "W m-2",
                       standard_name = "surface_upward_latent_heat_flux"),
    "Qh"           => (long_name = "Sensible Heat Flux",
                       units     = "W m-2",
                       standard_name = "surface_upward_sensible_heat_flux"),
    "GPP"          => (long_name = "Gross Primary Production",
                       units     = "kg C m-2 s-1",
                       standard_name = ""),
    "Reco"         => (long_name = "Ecosystem Respiration",
                       units     = "kg C m-2 s-1",
                       standard_name = ""),
    "TVeg"         => (long_name = "Vegetation Transpiration",
                       units     = "kg m-2 s-1",
                       standard_name = "transpiration_flux"),
    "ESoil"        => (long_name = "Bare Soil Evaporation",
                       units     = "kg m-2 s-1",
                       standard_name = ""),
    "Qg"           => (long_name = "Ground Heat Flux",
                       units     = "W m-2",
                       standard_name = "downward_heat_flux_in_soil"),
    "AvgSurfT"     => (long_name = "Average Surface Temperature",
                       units     = "K",
                       standard_name = "surface_temperature"),
    "SoilMoist"    => (long_name = "Total Column Soil Water",
                       units     = "kg m-2",
                       standard_name = "soil_moisture_content"),
    "LAI"          => (long_name = "Leaf Area Index",
                       units     = "m2 m-2",
                       standard_name = "leaf_area_index"),
    "TotAbovBioMass" => (long_name = "Total Above-Ground Biomass Carbon",
                         units     = "kg C m-2",
                         standard_name = ""),
    "TotSoilCarb"  => (long_name = "Total Soil Organic Carbon",
                       units     = "kg C m-2",
                       standard_name = ""),
)
const ALMA_VARS = [
    "NEE", "Qle", "Qh", "GPP", "Reco",
    "TVeg", "ESoil", "Qg", "AvgSurfT", "SoilMoist",
    "LAI", "TotAbovBioMass", "TotSoilCarb",
]
const UNCERTAINTY_VARS = ["NEE", "Qle", "Qh"]

# ── NetCDF writer ─────────────────────────────────────────────────────────────
"""Write a single CalLMIP-compliant NetCDF file.

`cal_or_val`  : "Cal" or "Val"
`prior_or_post` : "Prior" or "Posterior"
`vars`         : Dict of variable name → value vector
`out_dates`    : date vector for this period
"""
function write_callmip_nc(cal_or_val, prior_or_post, vars, out_dates)
    fname = "$(MODEL_NAME).$(MODEL_VER)_Expt$(EXPT_NO)_$(SITE_ID)_$(cal_or_val)_$(prior_or_post).nc"
    fpath = joinpath(nc_output_dir, fname)

    n = length(out_dates)
    # Time axis: days since 2003-01-01 (common epoch for this site)
    epoch      = Date(2003, 1, 1)
    time_vals  = Float64.(Dates.value.(out_dates .- epoch))
    time_units = "days since $(epoch)"

    NCDatasets.Dataset(fpath, "c") do ds
        # ── Global attributes ────────────────────────────────────────────
        ds.attrib["title"]        = "CalLMIP Phase 1 Output — $(SITE_ID)"
        ds.attrib["model"]        = MODEL_NAME
        ds.attrib["version"]      = MODEL_VER
        ds.attrib["site"]         = SITE_ID
        ds.attrib["experiment"]   = "CalLMIP Phase 1 Experiment $(EXPT_NO)"
        ds.attrib["period"]       = cal_or_val
        ds.attrib["parameters"]   = prior_or_post
        ds.attrib["Conventions"]  = "CF-1.8"
        ds.attrib["created_by"]   = "ClimaLand calibration pipeline"
        ds.attrib["created"]      = string(now())
        if prior_or_post == "Prior"
            ds.attrib["parameter_source"] = "Prior distribution mean"
        else
            ds.attrib["parameter_source"] = "EKI optimal (max-posterior point estimate)"
        end

        # ── Dimensions ──────────────────────────────────────────────────
        NCDatasets.defDim(ds, "time", n)

        # ── Time coordinate ──────────────────────────────────────────────
        v = NCDatasets.defVar(ds, "time", Float64, ("time",))
        v[:] = time_vals
        v.attrib["units"]     = time_units
        v.attrib["calendar"]  = "standard"
        v.attrib["long_name"] = "time"
        v.attrib["axis"]      = "T"

        # ── ALMA variables ───────────────────────────────────────────────
        for vname in ALMA_VARS
            haskey(vars, vname) || begin
                @warn "  Variable $vname not in data dict — skipping."
                continue
            end
            data = vars[vname]
            @assert length(data) >= length(out_dates) "$(vname): data shorter than dates"
            # Subset to the requested period (using external mask indices)
            # vars dict was already built for the full sim period;
            # out_dates tells us which rows to keep (handled below via idx)
            period_mask = cal_or_val == "Cal" ? cal_mask : val_mask
            period_data = data[period_mask]

            dv = NCDatasets.defVar(ds, vname, Float64, ("time",),
                                   fillvalue = -9999.0)
            dv[:] = period_data
            meta = get(ALMA_META, vname, nothing)
            if !isnothing(meta)
                dv.attrib["long_name"]    = meta.long_name
                dv.attrib["units"]        = meta.units
                meta.standard_name != "" &&
                    (dv.attrib["standard_name"] = meta.standard_name)
            end
            dv.attrib["cell_methods"] = "time: mean"
        end

        # ── Uncertainty percentiles (posterior files only) ───────────────
        if prior_or_post == "Posterior"
            period_mask = cal_or_val == "Cal" ? cal_mask : val_mask
            for vname in UNCERTAINTY_VARS
                for pct in ("p05", "p50", "p95")
                    key = "$(vname)_$(pct)"
                    haskey(vars, key) || continue
                    data = vars[key][period_mask]
                    dv   = NCDatasets.defVar(ds, key, Float64, ("time",),
                                             fillvalue = -9999.0)
                    dv[:] = data
                    meta  = get(ALMA_META, vname, nothing)
                    pct_label = Dict("p05" => "5th",
                                     "p50" => "50th",
                                     "p95" => "95th")[pct]
                    long = isnothing(meta) ? key :
                           meta.long_name * " ($(pct_label) percentile, posterior ensemble)"
                    units = isnothing(meta) ? "unknown" : meta.units
                    dv.attrib["long_name"]    = long
                    dv.attrib["units"]        = units
                    dv.attrib["cell_methods"] = "time: mean"
                end
            end
        end
    end  # NCDatasets.Dataset

    sz = filesize(fpath)
    @info "  Written $(fname)  ($(sz ÷ 1024) kB)"
    return fpath
end

# ── Write all 4 files ─────────────────────────────────────────────────────────
@info "Writing CalLMIP NetCDF files to $(nc_output_dir)…"

write_callmip_nc("Cal", "Prior",     prior_vars, dates[cal_mask])
write_callmip_nc("Cal", "Posterior", post_vars,  dates[cal_mask])
write_callmip_nc("Val", "Prior",     prior_vars, dates[val_mask])
write_callmip_nc("Val", "Posterior", post_vars,  dates[val_mask])

@info "Done! All CalLMIP NetCDF files written."
@info "Output directory: $(nc_output_dir)"
