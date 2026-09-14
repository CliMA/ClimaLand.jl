"""
Forward-model engine for the DK-Sor CalLMIP pipeline.

Builds the integrated ClimaLand `LandModel` (canopy + snow + soil + soil-CO2) at the
DK-Sor site and runs it for one or more years. The main entry point `run_prior_year`
returns either monthly (NEE, LHF, SHF) fluxes for the calibration objective or, with
`callmip=true`, the full set of daily CalLMIP diagnostics; it also supports multi-year
continuous runs, spin-up, and parameter/soil/radiation overrides. This file is included
by the calibration, posterior, and output scripts so they share one model definition.
Helpers to load FLUXNET observations and compute skill metrics are included.
"""

import ClimaComms
ClimaComms.@import_required_backends

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.TimeManager: date
import EnsembleKalmanProcesses as EKP     # needed before loading observations.jld2
using NCDatasets
using Dates
using Statistics
using LinearAlgebra
using Printf
using JLD2
using CairoMakie

# ── Configuration ──────────────────────────────────────────────────────────────
const FT          = Float64
const DT          = Float64(900)
const SITE_ID_VAL = :DK_Sor
const CLIMALAND_DIR = pkgdir(ClimaLand)
# Forcing + flux observations come from the CalLMIP ClimaArtifacts (see src/Artifacts.jl):
#   callmip_phase1_forcing → per-site *_Met.nc  (PLUMBER2/FLUXNET2015, from ME-org)
#   callmip_phase1         → per-site *_Flux.nc (callmip-org/Phase1; "1a" = DK-Sor test)
const MET_NC_PATH   = ClimaLand.Artifacts.callmip_phase1_forcing_path("DK-Sor")
const FLUX_NC_PATH  = ClimaLand.Artifacts.callmip_phase1_flux_path("DK-Sor"; phase = "1a")
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_prior_check")

# Years to check. Override with CALLMIP_CHECK_YEARS="1997,1998,..." to test the
# full span (e.g. the prior-crash years 2003 & 2011 that calibration flagged).
const CHECK_YEARS = haskey(ENV, "CALLMIP_CHECK_YEARS") ?
    parse.(Int, split(ENV["CALLMIP_CHECK_YEARS"], ",")) : [2005, 2008, 2010]

# CalLMIP daily diagnostic variables (run_prior_year(...; callmip=true)).
# soilrn/soillhf/soilshf (canonical soil-surface fluxes, exposed for LandModel in
# src/diagnostics) → hfg = −soilrn−soillhf−soilshf, evspsblsoi = soillhf/Lv.
const CALLMIP_SURFACE_VARS = ["nee","lhf","shf","gpp","er","hr","trans","ct","lai","cveg",
                              "soilrn","soillhf","soilshf"]
const CALLMIP_COLUMN_VARS  = ["swc","tsoil","soc"]

const _SITE_LOC    = FluxnetSimulations.get_location(FT, Val(SITE_ID_VAL))
const _SITE_HEIGHT = FluxnetSimulations.get_fluxtower_height(FT, Val(SITE_ID_VAL))
const _SITE_PARAMS = FluxnetSimulations.get_parameters(FT, Val(SITE_ID_VAL))

const lat         = _SITE_LOC.lat
const long        = _SITE_LOC.long
const time_offset = _SITE_LOC.time_offset
const atmos_h     = _SITE_HEIGHT.atmos_h

"""
Run one year, return (nee_monthly, lhf_monthly, shf_monthly).
`param_overrides` is an optional Dict(name => value); empty ⇒ default (prior)
parameters. Pass the posterior mean to evaluate the calibrated model.
"""
function run_prior_year(year; param_overrides::Dict = Dict{String,Float64}(),
                        dt = DT, callmip::Bool = false, with_co2::Bool = true,
                        return_sim::Bool = false, domain_override = nothing,
                        free_drainage::Bool = false, soil_overrides = (;),
                        rad_overrides = (;),
                        T_init_override = nothing, met_path_override = nothing,
                        spinup_days::Int = 0, stop_year::Int = year,
                        init_state = nothing, return_state::Bool = false)
    met_path = isnothing(met_path_override) ? MET_NC_PATH : met_path_override
    # spinup_days>0 starts the run earlier (e.g. a 60-day spin-up).
    # stop_year>year runs CONTINUOUSLY through (stop_year, 12, 31) in a single
    # integration (preserves inter-annual soil/carbon memory); default = one year.
    # NB: the :daily reduction can't emit the final day's mean (its window closes at
    # stop_date 00:00, beyond the forcing extent), so the last calendar day is fill.
    start_date = DateTime(year, 1, 1) - Day(spinup_days)
    stop_date  = DateTime(stop_year + 1, 1, 1)

    # Non-TOML params (site/radiation, e.g. α_NIR_leaf) are pulled out of
    # param_overrides and applied to rad_overrides below, not written to the TOML.
    param_overrides = Dict(param_overrides)
    for nm in ("α_NIR_leaf", "α_PAR_leaf", "ϵ_canopy")
        if haskey(param_overrides, nm)
            rad_overrides = merge((; Symbol(nm) => pop!(param_overrides, nm)), rad_overrides)
        end
    end

    if isempty(param_overrides)
        local_toml = LP.create_toml_dict(FT)          # default (prior) parameters
    else
        tmpf = tempname() * ".toml"
        open(tmpf, "w") do io
            for (nm, v) in param_overrides
                println(io, "[\"$nm\"]"); println(io, "value = $(Float64(v))")
                println(io, "type  = \"float\""); println(io)
            end
        end
        local_toml = LP.create_toml_dict(FT; override_files = [tmpf])
        rm(tmpf; force = true)
    end

    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, Val(SITE_ID_VAL))
    # Optional domain override (e.g. test the global model's 50 m / (10,0.5) /
    # 15-layer column, whose deep base stays below the seasonal freezing depth).
    if !isnothing(domain_override)
        haskey(domain_override, :zmin)     && (zmin = domain_override.zmin)
        haskey(domain_override, :nelements)&& (nelements = domain_override.nelements)
        haskey(domain_override, :dz_tuple) && (dz_tuple = FT.(domain_override.dz_tuple))
    end
    land_domain = Column(;
        zlim      = (FT(zmin), FT(zmax)),
        nelements, dz_tuple,
        longlat   = (long, lat),
    )
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    forcing = FluxnetSimulations.prescribed_forcing_netcdf(
        met_path, lat, long, time_offset, atmos_h,
        start_date, local_toml, FT,
    )

    LAI, maxLAI = NCDataset(met_path, "r") do ds
        time_vals = ds["time"][:]
        lai_data  = Float64.(coalesce.(ds["LAI_alternative"][1, 1, :], NaN))
        lai_secs  = Float64[
            Second(t - Hour(time_offset) - start_date).value
            for t in time_vals
        ]
        valid = .!isnan.(lai_data)
        TimeVaryingInput(lai_secs[valid], lai_data[valid]), maximum(lai_data[valid])
    end

    (; soil_ν, θ_r, soil_K_sat, soil_vg_n, soil_vg_α,
       ν_ss_om, ν_ss_quartz, ν_ss_gravel,
       z_0m_soil, z_0b_soil, soil_α_PAR, soil_α_NIR, soil_ϵ,
       Ω, χl, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
       ϵ_canopy, ac_canopy, Drel, g0,
       SAI, f_root_to_shoot, rooting_depth, h_canopy,
       conductivity_model, retention_model, plant_ν, plant_S_s) = _SITE_PARAMS

    # Optional soil-parameter overrides (e.g. test alternative values:
    # soil_K_sat=1e-5, ν_ss_om=0.03 — better-drained, lower-organic → higher κ).
    haskey(soil_overrides, :soil_K_sat) && (soil_K_sat = FT(soil_overrides.soil_K_sat))
    haskey(soil_overrides, :ν_ss_om)    && (ν_ss_om    = FT(soil_overrides.ν_ss_om))
    haskey(soil_overrides, :ν_ss_quartz)&& (ν_ss_quartz= FT(soil_overrides.ν_ss_quartz))
    haskey(soil_overrides, :ν_ss_gravel)&& (ν_ss_gravel= FT(soil_overrides.ν_ss_gravel))

    # Optional radiation/leaf-optics overrides (site params, not in the TOML) to
    # probe net-radiation control of SHF (α_NIR dominates SW_n absorbed energy).
    haskey(rad_overrides, :α_NIR_leaf) && (α_NIR_leaf = FT(rad_overrides.α_NIR_leaf))
    haskey(rad_overrides, :τ_NIR_leaf) && (τ_NIR_leaf = FT(rad_overrides.τ_NIR_leaf))
    haskey(rad_overrides, :α_PAR_leaf) && (α_PAR_leaf = FT(rad_overrides.α_PAR_leaf))
    haskey(rad_overrides, :ϵ_canopy)   && (ϵ_canopy   = FT(rad_overrides.ϵ_canopy))

    RAI = maxLAI * f_root_to_shoot

    prognostic_land_components =
        with_co2 ? (:canopy, :snow, :soil, :soilco2) : (:canopy, :snow, :soil)

    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR, NIR_albedo = soil_α_NIR)
    retention_parameters = (;
        ν      = soil_ν,
        θ_r,
        K_sat  = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

    # Bottom BC: default is the convenience-ctor zero-flux (water+heat); with
    # free_drainage=true use EnergyWaterFreeDrainage (water drains out the base,
    # the physical case for the permeable Sorø subsurface).
    bottom_bc_kw = free_drainage ?
        (; bottom_bc = Soil.EnergyWaterFreeDrainage()) : (;)
    soil = Soil.EnergyHydrology{FT}(
        land_domain, forcing, local_toml;
        prognostic_land_components,
        additional_sources   = (ClimaLand.RootExtraction{FT}(),),
        albedo               = soil_albedo,
        runoff               = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters, composition_parameters,
        S_s = FT(1e-3), z_0m = z_0m_soil, z_0b = z_0b_soil,
        bottom_bc_kw...,
    )

    soilco2 = if with_co2
        co2_drivers = Soil.Biogeochemistry.SoilDrivers(
            Soil.Biogeochemistry.PrognosticMet(soil.parameters), forcing.atmos,
        )
        Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, co2_drivers, local_toml)
    else
        nothing
    end

    radiation_parameters = (;
        Ω, G_Function = CLMGFunction(χl),
        α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
    )
    radiative_transfer = Canopy.TwoStreamModel{FT}(
        canopy_domain, local_toml; radiation_parameters, ϵ_canopy,
    )

    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI, SAI, RAI, rooting_depth, height = h_canopy,
    )
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        (; atmos = forcing.atmos, radiation = forcing.radiation,
           ground = ClimaLand.PrognosticGroundConditions{FT}()),
        LAI, local_toml;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis       = Canopy.PModel{FT}(canopy_domain, local_toml),
        conductance          = Canopy.PModelConductance{FT}(local_toml; Drel),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, local_toml; soil_params = (; ν = soil_ν, θ_r)),
        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            canopy_domain, local_toml;
            ν = plant_ν, S_s = plant_S_s, conductivity_model, retention_model),
        energy  = Canopy.BigLeafEnergyModel{FT}(local_toml; ac_canopy),
        biomass,
    )
    snow = Snow.SnowModel(
        FT, canopy_domain, forcing, local_toml, dt;
        prognostic_land_components,
    )
    land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)

    # Default IC temperature = met Tair at the start of the year (Jan → near
    # freezing), applied uniformly to the whole column. T_init_override lets us
    # test a realistic warm deep-soil temperature (annual mean ~281 K) — the
    # real deep soil is NOT at January's surface temperature.
    T_init_K = if !isnothing(T_init_override)
        Float64(T_init_override)
    else
        NCDataset(met_path, "r") do ds
            t_dates = DateTime.(ds["time"][:])
            idx_yr  = findfirst(t -> Dates.year(t) == year, t_dates)
            isnothing(idx_yr) ? 283.15 : Float64(coalesce(ds["Tair"][1, 1, idx_yr], 283.15))
        end
    end

    function set_ic!(Y, p, t, model)
        FT_l = eltype(Y.soil.ρe_int)
        if isnothing(init_state)
            # Cold start from arbitrary IC (soil near saturation, T = met Tair).
            ν_l  = FT_l(soil_ν)
            θr_l = FT_l(θ_r)
            Y.soil.ϑ_l .= θr_l + (ν_l - θr_l) * FT_l(0.95)
            Y.soil.θ_i .= FT_l(0)
            ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l, Y.soil.θ_i,
                soil.parameters.ρc_ds, soil.parameters.earth_param_set)
            Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
                Y.soil.θ_i, ρc_s, FT_l(T_init_K),
                soil.parameters.earth_param_set)
            Y.snow.S .= FT_l(0); Y.snow.S_l .= FT_l(0); Y.snow.U .= FT_l(0)
            Y.canopy.energy.T .= FT_l(T_init_K)
            Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν
        else
            # Restart the WATER + ENERGY (temperature) prognostic states from a
            # spun-up state (parent() copies are robust to space-identity). Carbon
            # is NOT carried — it is prescribed below (SOC/CO2/O2 from files).
            parent(Y.soil.ϑ_l)              .= parent(init_state.soil.ϑ_l)
            parent(Y.soil.θ_i)              .= parent(init_state.soil.θ_i)
            parent(Y.soil.ρe_int)           .= parent(init_state.soil.ρe_int)
            parent(Y.snow.S)                .= parent(init_state.snow.S)
            parent(Y.snow.S_l)              .= parent(init_state.snow.S_l)
            parent(Y.snow.U)                .= parent(init_state.snow.U)
            parent(Y.canopy.energy.T)       .= parent(init_state.canopy.energy.T)
            parent(Y.canopy.hydraulics.ϑ_l) .= parent(init_state.canopy.hydraulics.ϑ_l)
        end
        # Carbon pools are PRESCRIBED (from ClimaLand defaults), never spun up.
        if with_co2
            Y.soilco2.CO2 .= FT_l(6e-5)
            Y.soilco2.O2  .= FT_l(0.21)
            SOC_top = FT_l(15.0); SOC_bot = FT_l(0.5)
            τ_soc   = FT_l(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end

    ow = ClimaDiagnostics.Writers.DictWriter()
    if callmip
        diags = vcat(
            ClimaLand.default_diagnostics(land, start_date, "";
                output_writer = ow, output_vars = CALLMIP_SURFACE_VARS,
                reduction_period = :daily),
            ClimaLand.default_diagnostics(land, start_date, "";
                output_writer = ow, output_vars = CALLMIP_COLUMN_VARS,
                reduction_period = :daily),
        )
    else
        diags = ClimaLand.default_diagnostics(land, start_date, "";
            output_writer = ow, output_vars = ["lhf", "shf", "nee"],
            reduction_period = :daily)
    end
    simulation = LandSimulation(
        start_date, stop_date, dt, land;
        set_ic!, updateat = Second(dt), diagnostics = diags,
    )

    # Debug hook: return the (un-solved) simulation so a driver can step it
    # manually and inspect the soil state at the divergence onset.
    return_sim && return simulation

    solve!(simulation)

    # Spin-up hook: return the final prognostic state so a driver can restart the
    # output run from a spun-up (equilibrated water+temp) state.
    return_state && return simulation._integrator.u

    if callmip
        # Daily CalLMIP output: surface vars → Vector, column vars → (n_z × n_days).
        surface_data = Dict{String, Vector{Float64}}()
        ref_dates    = Date[]
        for var in CALLMIP_SURFACE_VARS
            try
                times, data =
                    ClimaLand.Diagnostics.diagnostic_as_vectors(ow, "$(var)_1d_average")
                surface_data[var] = Float64.(data)
                # ClimaDiagnostics timestamps each :daily mean at the END of its
                # averaging window, so the value written at 00:00 of day D+1 is the
                # mean over day D. Subtract 1 day so each daily mean is labelled by
                # the calendar day it represents (verified empirically: the model
                # otherwise leads the FLUXNET obs by exactly 1 day). This also fills
                # 1997-01-01 and drops the spurious trailing 2015-01-01.
                isempty(ref_dates) && (ref_dates = Date.(date.(times)) .- Day(1))
            catch
                @warn "  [$year] surface '$var' not found"
            end
        end
        column_data = Dict{String, Matrix{Float64}}()
        for var in CALLMIP_COLUMN_VARS
            key = "$(var)_1d_average"
            if haskey(ow.dict, key)
                raw = ow.dict[key]
                tt  = sort(collect(keys(raw)))
                column_data[var] =
                    reduce(hcat, [Float64.(vec(parent(raw[t]))) for t in tt])
            else
                @warn "  [$year] column '$var' not found"
            end
        end
        z_soil = try
            Float64.(vec(parent(ClimaCore.Fields.coordinate_field(
                axes(land.soil.domain.fields.z)).z)))
        catch
            Float64[]
        end
        return surface_data, column_data, z_soil, ref_dates
    end

    function get_monthly(diag_name, scale = 1.0)
        times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(ow, diag_name)
        months = Dates.month.(date.(times))
        return [mean(data[months .== m]) * scale for m in 1:12]
    end

    nee_m = get_monthly("nee_1d_average", 12.0 * 86400.0)  # mol/m²/s → gC/m²/d
    lhf_m = get_monthly("lhf_1d_average")
    shf_m = get_monthly("shf_1d_average")
    return nee_m, lhf_m, shf_m
end

"""Load monthly obs for a given year from the FLUXNET daily NC file."""
function load_obs_year(year)
    NCDataset(FLUX_NC_PATH, "r") do ds
        dates   = DateTime.(ds["time"][:])
        nee_all = Float64.(coalesce.(ds["NEE_daily"][:], NaN))
        qle_all = Float64.(coalesce.(ds["Qle_daily"][:], NaN))
        qh_all  = Float64.(coalesce.(ds["Qh_daily"][:],  NaN))
        for arr in (nee_all, qle_all, qh_all)
            arr[arr .>= 1.0e19] .= NaN
        end

        nee_m = Float64[]; lhf_m = Float64[]; shf_m = Float64[]
        for mon in 1:12
            mask = (Dates.year.(dates) .== year) .& (Dates.month.(dates) .== mon)
            nee_v = nee_all[mask]; qle_v = qle_all[mask]; qh_v = qh_all[mask]
            nee_f = filter(isfinite, nee_v); push!(nee_m, isempty(nee_f) ? NaN : mean(nee_f))
            lhf_f = filter(isfinite, qle_v); push!(lhf_m, isempty(lhf_f) ? NaN : mean(lhf_f))
            shf_f = filter(isfinite, qh_v);  push!(shf_m, isempty(shf_f) ? NaN : mean(shf_f))
        end
        return nee_m, lhf_m, shf_m
    end
end

function rmse_r2(mod, obs)
    mask = isfinite.(mod) .& isfinite.(obs)
    sum(mask) < 2 && return NaN, NaN
    m = mod[mask]; o = obs[mask]
    rmse = sqrt(mean((m .- o).^2))
    r2   = 1.0 - sum((m .- o).^2) / sum((o .- mean(o)).^2)
    return rmse, r2
end

# ── Run simulations ────────────────────────────────────────────────────────────
# Only when run directly (not when `include`d by compute_prior_post_stats.jl,
# which reuses run_prior_year / load_obs_year / rmse_r2).
if abspath(PROGRAM_FILE) == @__FILE__
mkpath(OUTDIR)

results_mod = Dict{Int, NamedTuple}()
results_obs = Dict{Int, NamedTuple}()

ok_years = Int[]
for yr in CHECK_YEARS
    @info "Running year $yr with default parameters..."
    try
        nee_m, lhf_m, shf_m = run_prior_year(yr)
        # Sanity: physical fluxes, no NaN/blowup.
        finite = all(isfinite, nee_m) && all(isfinite, lhf_m) && all(isfinite, shf_m)
        @info @sprintf("  %s year %d  NEE[%.2f,%.2f] LHF[%.1f,%.1f] SHF[%.1f,%.1f]",
            finite ? "OK  " : "BAD ", yr,
            minimum(nee_m), maximum(nee_m),
            minimum(lhf_m), maximum(lhf_m),
            minimum(shf_m), maximum(shf_m))
        if finite
            results_mod[yr] = (; nee = nee_m, lhf = lhf_m, shf = shf_m)
            nee_o, lhf_o, shf_o = load_obs_year(yr)
            results_obs[yr] = (; nee = nee_o, lhf = lhf_o, shf = shf_o)
            push!(ok_years, yr)
        end
    catch e
        @warn "  FAIL year $yr: $(typeof(e)): $(sprint(showerror, e))"
    end
end

@info "PRIOR FORWARD-MODEL CHECK: $(length(ok_years))/$(length(CHECK_YEARS)) years OK: $ok_years"
failed = setdiff(CHECK_YEARS, ok_years)
isempty(failed) || @warn "FAILED years: $failed"

# Only plot if every requested year ran (keeps the multi-year metrics honest).
if length(ok_years) != length(CHECK_YEARS)
    @info "Skipping plot (not all years succeeded); foundational check is the result."
    exit(0)
end

# ── Aggregate metrics across all 3 years ─────────────────────────────────────
function aggregate_metric(var_key)
    mod_all = vcat([results_mod[yr][var_key] for yr in CHECK_YEARS]...)
    obs_all = vcat([results_obs[yr][var_key] for yr in CHECK_YEARS]...)
    return rmse_r2(mod_all, obs_all)
end

nee_rmse, nee_r2 = aggregate_metric(:nee)
lhf_rmse, lhf_r2 = aggregate_metric(:lhf)
shf_rmse, shf_r2 = aggregate_metric(:shf)

@info "Prior performance (2005/2008/2010 combined):"
@info "  NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))"
@info "  LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))"
@info "  SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))"

# ── Plot ───────────────────────────────────────────────────────────────────────
months = 1:12
var_labels = ["NEE (gC m⁻² d⁻¹)", "LHF (W m⁻²)", "SHF (W m⁻²)"]
var_keys   = [:nee, :lhf, :shf]
metrics    = [(nee_rmse, nee_r2), (lhf_rmse, lhf_r2), (shf_rmse, shf_r2)]

fig = Figure(size = (1000, 800))

for (row, (vk, ylabel, (rmse, r2))) in enumerate(zip(var_keys, var_labels, metrics))
    for (col, yr) in enumerate(CHECK_YEARS)
        ax = Axis(fig[row, col];
            xlabel = col == 2 ? "Month" : "",
            ylabel = col == 1 ? ylabel : "",
            title  = row == 1 ? string(yr) : "",
            xticks = (1:12, string.(1:12)),
        )

        obs_v = results_obs[yr][vk]
        mod_v = results_mod[yr][vk]

        # Obs: black dots for valid months
        valid_months = findall(isfinite, obs_v)
        scatter!(ax, valid_months, obs_v[valid_months];
            color = :black, markersize = 8, label = "FLUXNET")

        # Model: colored line
        lines!(ax, months, mod_v; color = Makie.wong_colors()[row], linewidth = 2,
            label = "ClimaLand prior")

        if col == 1 && row == 1
            axislegend(ax; position = :rt, labelsize = 10)
        end
    end
end

# Overall title with metrics
Label(fig[0, :],
    "DK-Sor prior (default params) vs FLUXNET — 2005/2008/2010\n" *
    "NEE: RMSE=$(round(nee_rmse; digits=2)) gC/m²/d, R²=$(round(nee_r2; digits=3))  |  " *
    "LHF: RMSE=$(round(lhf_rmse; digits=2)) W/m², R²=$(round(lhf_r2; digits=3))  |  " *
    "SHF: RMSE=$(round(shf_rmse; digits=2)) W/m², R²=$(round(shf_r2; digits=3))",
    fontsize = 13, font = :bold,
)

save(joinpath(OUTDIR, "prior_vs_obs.png"), fig; px_per_unit = 2)
@info "Plot saved → $(joinpath(OUTDIR, "prior_vs_obs.png"))"

end  # if abspath(PROGRAM_FILE) == @__FILE__
