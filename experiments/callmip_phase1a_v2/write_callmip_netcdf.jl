"""
Write CalLMIP Phase 1a compliant NetCDF output for DK-Sor.

Produces 2 files (prior + posterior), each covering the full 1997-01-01 to
2014-12-31 period (6574 days):

  callmip_output/ClimaLand._Phase1a_Scen1_DK-Sor_Cal_Prior.nc
  callmip_output/ClimaLand._Phase1a_Scen1_DK-Sor_Cal_Posterior.nc

Variable naming: CMIP conventions (not ALMA).
All variables have dimensions (time, lat, lon) with lat=1, lon=1.
Lat and lon are scalar coordinate variables with long_name + units attributes.

Sign convention: nep = GPP - ER = -NEE (positive = carbon sink).
The CalLMIP flux obs file has NEE positive = source; we flip sign for nep.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/write_callmip_netcdf.jl
"""

using Dates
import JLD2
import NCDatasets

# ── Configuration ──────────────────────────────────────────────────────────────
const SITE_ID    = "DK-Sor"
const SITE_LAT   = 55.4859      # degrees N
const SITE_LON   = 11.6446      # degrees E
const N_DAYS     = 6574         # 1997-01-01 – 2014-12-31 inclusive
const START_DATE = Date(1997, 1, 1)
const STOP_DATE  = Date(2014, 12, 31)

const climaland_dir     = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir           = @__DIR__
const callmip_sim_dir   = joinpath(exp_dir, "output_callmip_sims")
const nc_output_dir     = joinpath(exp_dir, "callmip_output")
isdir(nc_output_dir) || mkpath(nc_output_dir)

# ── Unit conversion constants ──────────────────────────────────────────────────
const mol_CO2_to_kgC = 12e-3     # mol CO2/m²/s → kgC/m²/s
const Lv_water       = 2.5e6     # J/kg (latent heat of vapourisation)
const rho_water      = 1000.0    # kg/m³

# ── CMIP variable metadata ─────────────────────────────────────────────────────
const CMIP_META = Dict{String, NamedTuple}(
    "nep"          => (long_name = "Net Ecosystem Production (GPP - ER)",
                       units     = "kgC m-2 s-1"),
    "hfls"         => (long_name = "Latent Heat Flux",
                       units     = "W m-2"),
    "hfss"         => (long_name = "Sensible Heat Flux",
                       units     = "W m-2"),
    "gpp"          => (long_name = "Gross Primary Production",
                       units     = "kgC m-2 s-1"),
    "reco"         => (long_name = "Total Ecosystem Respiration (Ra + Rh)",
                       units     = "kgC m-2 s-1"),
    "tran"         => (long_name = "Transpiration",
                       units     = "kg m-2 s-1"),
    "evspsblsoi"   => (long_name = "Bare Soil Evaporation",
                       units     = "kg m-2 s-1"),
    "hfg"          => (long_name = "Ground Heat Flux",
                       units     = "W m-2"),
    "ts"           => (long_name = "Land Surface Skin Temperature",
                       units     = "K"),
    "lai"          => (long_name = "Leaf Area Index",
                       units     = "m2 m-2"),
    "cLiveBioAbove" => (long_name = "Aboveground Live Biomass Carbon",
                        units     = "kgC m-2"),
    "cSoil"        => (long_name = "Soil Organic Carbon (total column)",
                       units     = "kgC m-2"),
    "mrso"         => (long_name = "Total Column Soil Moisture",
                       units     = "kg m-2"),
    "mrsos"        => (long_name = "Top 10 cm Soil Moisture",
                       units     = "kg m-2"),
    "nep_post_unc" => (long_name = "Posterior Uncertainty on nep",
                       units     = "kgC m-2 s-1"),
    "hfls_post_unc" => (long_name = "Posterior Uncertainty on hfls",
                        units     = "W m-2"),
    "hfss_post_unc" => (long_name = "Posterior Uncertainty on hfss",
                        units     = "W m-2"),
)

# Output order (determines variable order in file)
const CMIP_VARS = [
    "nep", "hfls", "hfss", "gpp", "reco",
    "tran", "evspsblsoi", "hfg", "ts", "lai",
    "cLiveBioAbove", "cSoil", "mrso", "mrsos",
]

# ── Load simulation diagnostics ───────────────────────────────────────────────
function load_diagnostics(member::Int)
    path = joinpath(callmip_sim_dir,
                    "iteration_000",
                    "member_$(lpad(member, 3, '0'))",
                    "callmip_diagnostics.jld2")
    isfile(path) || error("Diagnostics not found: $path\nRun run_callmip_simulations.jl first.")
    JLD2.load(path)
end

@info "Loading prior diagnostics (member 001)…"
prior_d = load_diagnostics(1)
@info "Loading posterior diagnostics (member 002)…"
post_d  = load_diagnostics(2)

# ── Expand to full 1997-2014 date vector ─────────────────────────────────────
all_dates = collect(START_DATE:Day(1):STOP_DATE)
@assert length(all_dates) == N_DAYS "Expected $N_DAYS days, got $(length(all_dates))"

# ── Column integration helper ─────────────────────────────────────────────────
function column_integral(col_data, z_soil; scale = 1.0, top_m = nothing)
    isempty(col_data) && return fill(NaN, N_DAYS)
    n_z, n_t = size(col_data)
    n_zl = length(z_soil)
    dz = zeros(n_zl)
    if n_zl == 1
        dz[1] = abs(z_soil[1])
    else
        for i in 1:n_zl
            lo = i == 1      ? 2*z_soil[1] - z_soil[2]    : (z_soil[i-1]+z_soil[i])/2
            hi = i == n_zl   ? 0.0                         : (z_soil[i]+z_soil[i+1])/2
            dz[i] = hi - lo
        end
    end
    # For top-only integration (mrsos = top 10 cm)
    top_mask = isnothing(top_m) ? trues(n_zl) : (-z_soil .<= top_m)

    result = zeros(n_t)
    if n_z == n_zl
        for t in 1:n_t, i in 1:n_zl
            top_mask[i] || continue
            result[t] += col_data[i, t] * dz[i] * scale
        end
    elseif n_z == 1
        total = top_m === nothing ? sum(dz) : sum(dz[top_mask])
        for t in 1:n_t
            result[t] = col_data[1, t] * total * scale
        end
    else
        @warn "column_integral: mismatch n_z=$n_z vs n_zl=$n_zl — returning NaN"
        return fill(NaN, n_t)
    end
    return result
end

# ── Build variable dictionary from diagnostics ────────────────────────────────
function build_variable_dict(d::Dict, post_unc::Bool = false)
    out = Dict{String, Union{Vector{Float64}, Nothing}}()
    dates_sim = d["dates"]           # simulation dates (1997–2014 after spinup trim)
    sd = get(d, "surface_data", Dict{String,Any}())
    cd = get(d, "column_data",  Dict{String,Any}())
    z_soil = haskey(d, "z_soil") ? Float64.(d["z_soil"]) : Float64[]

    function get_surf(key)
        v = get(sd, key, nothing)
        (v !== nothing && length(v) > 0) ? Float64.(v) : nothing
    end

    # Align to full 1997-2014 date vector (fill missing days with NaN)
    function align_to_full(vals, sim_dates)
        isnothing(vals) && return fill(NaN, N_DAYS)
        d2v = Dict(zip(sim_dates, vals))
        [get(d2v, dt, NaN) for dt in all_dates]
    end

    # ── Carbon fluxes (mol CO2/m²/s → kgC/m²/s) ───────────────────────────
    nee_raw = align_to_full(get_surf("nee"), dates_sim)
    er_raw  = align_to_full(get_surf("er"),  dates_sim)
    gpp_raw = align_to_full(get_surf("gpp"), dates_sim)

    # Blowup guard: mask |flux| > 1e-3 mol/m²/s
    nee_raw[abs.(nee_raw) .> 1e-3] .= NaN
    er_raw[abs.(er_raw)   .> 1e-3] .= NaN
    gpp_raw[abs.(gpp_raw) .> 1e-3] .= NaN

    # nep = -NEE (positive = sink), convert mol CO2/m²/s → kgC/m²/s
    out["nep"]  = -nee_raw .* mol_CO2_to_kgC
    out["gpp"]  = gpp_raw  .* mol_CO2_to_kgC
    out["reco"] = er_raw   .* mol_CO2_to_kgC

    # ── Heat fluxes ────────────────────────────────────────────────────────
    out["hfls"] = align_to_full(get_surf("lhf"),  dates_sim)
    out["hfss"] = align_to_full(get_surf("shf"),  dates_sim)

    # ── Transpiration (kg/m²/s) ────────────────────────────────────────────
    out["tran"] = align_to_full(get_surf("trans"), dates_sim)

    # ── Bare soil evaporation ──────────────────────────────────────────────
    # Not available in current ClimaLand.jl — write NaN
    out["evspsblsoi"] = nothing

    # ── Ground heat flux ───────────────────────────────────────────────────
    # Approximate from energy balance if soil fluxes available, else NaN
    soilrn  = get_surf("soilrn")
    soillhf = get_surf("soillhf")
    soilshf = get_surf("soilshf")
    if !isnothing(soilrn) && !isnothing(soillhf) && !isnothing(soilshf)
        qg_raw = soilrn .- soillhf .- soilshf
        out["hfg"] = align_to_full(qg_raw, dates_sim)
    else
        out["hfg"] = nothing
    end

    # ── Surface temperature ────────────────────────────────────────────────
    out["ts"] = align_to_full(get_surf("ct"), dates_sim)

    # ── LAI ───────────────────────────────────────────────────────────────
    out["lai"] = align_to_full(get_surf("lai"), dates_sim)

    # ── Aboveground biomass carbon (not available) ─────────────────────────
    out["cLiveBioAbove"] = nothing

    # ── Soil organic carbon (total column, kgC/m²) ─────────────────────────
    if haskey(cd, "soc") && !isempty(cd["soc"]) && !isempty(z_soil)
        soc_col = Float64.(cd["soc"])
        soc_sim = column_integral(soc_col, z_soil)
        out["cSoil"] = align_to_full(soc_sim, dates_sim)
    else
        out["cSoil"] = nothing
    end

    # ── Total column soil moisture (mrso, kg/m²) ──────────────────────────
    if haskey(cd, "swc") && !isempty(cd["swc"]) && !isempty(z_soil)
        swc_col = Float64.(cd["swc"])
        mrso_sim = column_integral(swc_col, z_soil; scale = rho_water)
        out["mrso"] = align_to_full(mrso_sim, dates_sim)
        # Top 10 cm
        mrsos_sim = column_integral(swc_col, z_soil; scale = rho_water, top_m = 0.10)
        out["mrsos"] = align_to_full(mrsos_sim, dates_sim)
    else
        out["mrso"]  = nothing
        out["mrsos"] = nothing
    end

    # ── Posterior uncertainty (only in posterior file) ─────────────────────
    if post_unc
        out["nep_post_unc"]  = nothing   # filled from analyze_posterior_ensemble if available
        out["hfls_post_unc"] = nothing
        out["hfss_post_unc"] = nothing
    end

    return out
end

@info "Building variable dictionaries…"
prior_vars = build_variable_dict(prior_d, false)
post_vars  = build_variable_dict(post_d,  true)

# ── NetCDF writer ─────────────────────────────────────────────────────────────
function write_callmip_nc(vardata::Dict, stage::String)
    # stage: "Prior" or "Posterior"
    fname = "ClimaLand._Phase1a_Scen1_$(SITE_ID)_Cal_$(stage).nc"
    fpath = joinpath(nc_output_dir, fname)

    ds = NCDatasets.NCDataset(fpath, "c")

    # ── Global attributes ────────────────────────────────────────────────────
    ds.attrib["Model"]                = "ClimaLand."
    ds.attrib["CalLMIP_Phase"]        = "Phase1a"
    ds.attrib["Calibration_Scenario"] = "Scen1"
    ds.attrib["Calibration_stage"]    = stage
    ds.attrib["Cal_Val"]              = "Calibration"
    ds.attrib["Conventions"]          = "CF-1.8"
    ds.attrib["site"]                 = SITE_ID
    ds.attrib["created"]              = string(now())

    # ── Dimensions ───────────────────────────────────────────────────────────
    NCDatasets.defDim(ds, "time", N_DAYS)
    NCDatasets.defDim(ds, "lat",  1)
    NCDatasets.defDim(ds, "lon",  1)

    # ── Latitude coordinate variable (required by ME-org) ────────────────────
    lat_v = NCDatasets.defVar(ds, "lat", Float64, ("lat",))
    lat_v[:] = [SITE_LAT]
    lat_v.attrib["units"]     = "degrees_north"
    lat_v.attrib["long_name"] = "latitude"

    # ── Longitude coordinate variable ────────────────────────────────────────
    lon_v = NCDatasets.defVar(ds, "lon", Float64, ("lon",))
    lon_v[:] = [SITE_LON]
    lon_v.attrib["units"]     = "degrees_east"
    lon_v.attrib["long_name"] = "longitude"

    # ── Time coordinate ───────────────────────────────────────────────────────
    # Units: "days since 1997-01-01 00:00:00", increments 1,2,...,6574
    # NOT fractional days — ME-org requires integer day steps.
    time_v = NCDatasets.defVar(ds, "time", Float64, ("time",))
    time_v[:] = collect(1.0:Float64(N_DAYS))
    time_v.attrib["units"]     = "days since 1997-01-01 00:00:00"
    time_v.attrib["calendar"]  = "standard"
    time_v.attrib["long_name"] = "time"

    # ── Data variables: all (time, lat, lon) ──────────────────────────────────
    all_varnames = stage == "Posterior" ?
        [CMIP_VARS; ["nep_post_unc", "hfls_post_unc", "hfss_post_unc"]] : CMIP_VARS

    for varname in all_varnames
        meta  = get(CMIP_META, varname, nothing)
        units = isnothing(meta) ? "unknown" : meta.units
        long  = isnothing(meta) ? varname   : meta.long_name

        vals  = get(vardata, varname, nothing)
        v = NCDatasets.defVar(ds, varname, Float64, ("time", "lat", "lon"))

        if isnothing(vals) || all(isnan.(vals))
            v[:, 1, 1] = fill(NaN, N_DAYS)
            v.attrib["note"] = "not available in current ClimaLand.jl version"
        else
            @assert length(vals) == N_DAYS "$(varname): expected $N_DAYS days, got $(length(vals))"
            v[:, 1, 1] = vals
        end
        v.attrib["units"]     = units
        v.attrib["long_name"] = long
        v.attrib["cell_methods"] = "time: mean"
    end

    close(ds)
    sz = filesize(fpath)
    @info "Written: $fname  ($(sz ÷ 1024) kB)"
    return fpath
end

# ── Write the 2 files ─────────────────────────────────────────────────────────
@info "Writing CalLMIP NetCDF files to $nc_output_dir…"
prior_nc = write_callmip_nc(prior_vars, "Prior")
post_nc  = write_callmip_nc(post_vars,  "Posterior")

# ── Validate both files ───────────────────────────────────────────────────────
function validate_callmip_nc(fname)
    ds = NCDatasets.NCDataset(fname)

    @assert ds.dim["time"] == N_DAYS    "Expected $N_DAYS days, got $(ds.dim["time"])"
    @assert ds.dim["lat"]  == 1          "Expected lat dim = 1"
    @assert ds.dim["lon"]  == 1          "Expected lon dim = 1"

    @assert haskey(ds, "lat")            "lat must be a variable"
    @assert haskey(ds, "lon")            "lon must be a variable"
    @assert ds["lat"].attrib["units"]     == "degrees_north"
    @assert ds["lat"].attrib["long_name"] == "latitude"
    @assert ds["lon"].attrib["units"]     == "degrees_east"
    @assert ds["lon"].attrib["long_name"] == "longitude"
    @assert ds["lat"][:] isa Vector{<:Real}
    @assert ds["lon"][:] isa Vector{<:Real}

    @assert haskey(ds, "time")            "time variable missing"
    @assert ds["time"].attrib["units"]    == "days since 1997-01-01 00:00:00"
    tv = Float64.(ds["time"][:])
    @assert tv ≈ collect(1.0:Float64(N_DAYS)) "Time axis must be integer days 1:$N_DAYS"

    for var in CMIP_VARS
        haskey(ds, var) || begin
            close(ds)
            error("Required variable $var missing from $(basename(fname))")
        end
        data = Float64.(ds[var][:, 1, 1])
        n_valid = sum(!isnan, data)
        if n_valid == 0
            @warn "$var — all NaN in $(basename(fname))"
        else
            n_valid / N_DAYS < 0.9 && @warn "$var — $(N_DAYS - n_valid) NaN days ($(round(100*(N_DAYS-n_valid)/N_DAYS; digits=1))%)"
        end
    end

    for attr in ["Model", "CalLMIP_Phase", "Calibration_Scenario",
                 "Calibration_stage", "Cal_Val"]
        haskey(ds.attrib, attr) || begin
            close(ds)
            error("Missing global attribute: $attr")
        end
    end

    close(ds)
    println("✓ $(basename(fname)) passed all validation checks")
end

println("\nValidating output files:")
validate_callmip_nc(prior_nc)
validate_callmip_nc(post_nc)

println("\nDone!")
println("  Prior:     $(basename(prior_nc))")
println("  Posterior: $(basename(post_nc))")
println("  Output dir: $nc_output_dir")
