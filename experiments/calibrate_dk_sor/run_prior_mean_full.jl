"""
Run the calibration forward model at the prior mean parameters for the full
period 1997–2014 (spinup 1997, analysis 1998–2013), as required by ME-org
(CalLMIP Phase 1, 6574-day period).

Uses prior mean values passed directly to constructors (avoids TOML unicode issues).
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date
import ClimaParams as CP
using Insolation

using Dates
using NCDatasets
using Statistics
using CairoMakie
CairoMakie.activate!()

# Needed to access prescribed_forcing_netcdf internals for recycled spinup
import ClimaLand.Parameters as LP

const FT = Float64
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const SITE_ID = "DK-Sor"
const DT = Float64(900)

# Pass --smoke (or SMOKE=true env var) for a quick validation run:
# runs full 1996 spinup + Jan 1997 analysis so NEE/H/LE are non-NaN
const SMOKE_TEST = "--smoke" in ARGS || get(ENV, "SMOKE", "false") == "true"

# ── Prior mean values (from run_calibration.jl) ───────────────────────────────
# Centers match ClimaLand toml/default_parameters.toml exactly.
const PRIOR_TOML = joinpath(@__DIR__, "prior_mean_parameters.toml")
open(PRIOR_TOML, "w") do io
    write(io, """
[moisture_stress_c]
value = 0.27
type = "float"
used_in = ["getindex"]

["pmodel_cstar"]
value = 0.43
type = "float"
used_in = ["getindex"]

["pmodel_β"]
value = 20.0
type = "float"
used_in = ["getindex"]

["leaf_Cd"]
value = 0.07
type = "float"
used_in = ["getindex"]

["canopy_z_0m_coeff"]
value = 0.10
type = "float"
used_in = ["getindex"]

["canopy_z_0b_coeff"]
value = 0.0007
type = "float"
used_in = ["getindex"]

["canopy_d_coeff"]
value = 0.65
type = "float"
used_in = ["getindex"]

["canopy_K_lw"]
value = 0.85
type = "float"
used_in = ["getindex"]

["canopy_emissivity"]
value = 0.98
type = "float"
used_in = ["getindex"]

["root_leaf_nitrogen_ratio"]
value = 1.0
type = "float"
used_in = ["getindex"]

["stem_leaf_nitrogen_ratio"]
value = 0.1
type = "float"
used_in = ["getindex"]

["soilCO2_activation_energy"]
value = 61000.0
type = "float"
used_in = ["Land"]

["soilCO2_pre_exponential_factor"]
value = 23835.0
type = "float"
used_in = ["Land"]

["michaelis_constant"]
value = 0.005
type = "float"
used_in = ["Land"]

["O2_michaelis_constant"]
value = 0.004
type = "float"
used_in = ["Land"]

["ac_canopy"]
value = 2500.0
type = "float"
used_in = ["getindex"]
""")
end
println("Wrote prior mean TOML: $PRIOR_TOML")

# ── Setup ─────────────────────────────────────────────────────────────────────
site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)

# Spinup: 1996-01-01 to 1997-01-01 using recycled 1997 forcing (1996 is a leap year = 366 days)
# Analysis: 1997-01-01 to 2015-01-01 (= 1997-2014 inclusive = 6574 days for ME-org)
const SPINUP_SHIFT_DAYS = 366  # 1996 is a leap year
sim_start   = DateTime(1996, 1, 1)
sim_stop    = SMOKE_TEST ? DateTime(1997, 2, 1) : DateTime(2015, 1, 1)
spinup_date = DateTime(1997, 1, 1)
SMOKE_TEST && println("*** SMOKE TEST MODE: 1996 spinup + Jan 1997 analysis → $(sim_stop) ***")
println("Simulating $sim_start → $sim_stop (spinup 1996 from recycled 1997 data, analysis from $spinup_date)")

toml_dict = LP.create_toml_dict(FT; override_files = [PRIOR_TOML])

(; dz_tuple, nelements, zmin, zmax) = FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long)          = FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h)                         = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

land_domain   = Column(; zlim = (zmin, zmax), nelements, dz_tuple, longlat = (long, lat))
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

met_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_1997-2014_FLUXNET2015_Met.nc")

# ── Build recycled forcing: prepend shifted 1997 data for the 1996 spinup year ──
# The natural data starts at 1997-01-01. We shift those 1997 timestamps back by
# 366 days (leap year 1996) so the spinup period 1996-01-01→1997-01-01 uses real
# 1997 meteorology. Then the actual 1997-2014 data follows seamlessly.
#
# prescribed_forcing_netcdf uses start_date to compute seconds_since_start.
# We pass sim_start = 1997-01-01 so it builds the natural [0, …] time axis,
# then we manually shift all timestamps to [+366d, …] and prepend the recycled
# spinup block at [0, …+366d).  Finally we pass the combined arrays directly
# to TimeVaryingInput and build PrescribedAtmosphere ourselves.

# Step 1: load raw met data (seconds relative to 1997-01-01)
period_offset = Hour(time_offset)
_ds = NCDataset(met_nc_path, "r")
_times = _ds["time"][:]
function _read_var(varname)
    v = _ds[varname]
    ndims(v) == 3 ? Float64.(coalesce.(v[1,1,:], NaN)) :
    ndims(v) == 1 ? Float64.(coalesce.(v[:], NaN)) :
    error("Unexpected dims for $varname")
end
_T    = _read_var("Tair")
_SW   = _read_var("SWdown")
_LW   = _read_var("LWdown")
_q    = _read_var("Qair")
_P    = _read_var("Psurf")
_prec = _read_var("Precip")
_wind = _read_var("Wind")
_co2  = _read_var("CO2air")
_lai_raw  = Float64.(coalesce.(_ds["LAI"][1,1,:], NaN))
close(_ds)

# seconds since 1997-01-01 (the actual data reference)
_ref_date = DateTime(1997, 1, 1)
_t97 = Float64[Second(t - period_offset - _ref_date).value for t in _times]

# Step 2: extract 1997 data only (t in [0, 366 days))
_shift_s = Float64(SPINUP_SHIFT_DAYS * 86400)
_mask97  = (_t97 .>= 0.0) .& (_t97 .< _shift_s)
# _t97[_mask97] is in [0, 366d) relative to 1997-01-01.
# sim_start is 1996-01-01 = 366d before 1997-01-01, so _t97 is already
# in [0, 366d) relative to sim_start — use as-is for the spinup block.
_t_spinup = _t97[_mask97]  # → [0, 366d) relative to sim_start

# Step 3: full time axis for combined dataset (spinup + actual), relative to sim_start = 1996-01-01
# sim_start is 366 days before 1997-01-01, so offset all _t97 by +366d
_t_full_s = _shift_s  # seconds from 1996-01-01 to 1997-01-01
_t_actual = _t97 .+ _t_full_s  # actual 1997-2014 data shifted to sim_start basis (→ [366d, 18y))

# Combine: [recycled 1997 for spinup] ++ [actual 1997-2014]
# NOTE: _t97 values for UTC-offset timestamps at the start of the file can be
# negative (e.g., Denmark UTC+1 shifts 1996-12-31 23:00 UTC into _t97 = -3600s).
# Those negative-_t97 points land at _t_actual = shift_s + negative = < shift_s,
# overlapping the spinup range [0, shift_s).  Filter them out with _t97 >= 0.
function _make_tvi(t_actual_vals, data, t_spinup_vals)
    valid_sp  = .!isnan.(data[_mask97])
    valid_act = (.!isnan.(data)) .& (_t97 .>= 0.0)   # exclude pre-1997 UTC-offset points
    t_comb = vcat(t_spinup_vals[valid_sp], t_actual_vals[valid_act])
    v_comb = vcat(data[_mask97][valid_sp], data[valid_act])
    return TimeVaryingInput(t_comb, v_comb)
end

_make_co2(t_actual_vals, data, t_spinup_vals) = begin
    valid_sp  = .!isnan.(data[_mask97])
    valid_act = (.!isnan.(data)) .& (_t97 .>= 0.0)   # exclude pre-1997 UTC-offset points
    t_comb = vcat(t_spinup_vals[valid_sp], t_actual_vals[valid_act])
    v_comb = vcat(data[_mask97][valid_sp] .* FT(1e-6), data[valid_act] .* FT(1e-6))
    TimeVaryingInput(t_comb, v_comb)
end

atmos_T    = _make_tvi(_t_actual, _T,    _t_spinup)
atmos_P    = _make_tvi(_t_actual, _P,    _t_spinup)
atmos_u    = _make_tvi(_t_actual, _wind, _t_spinup)
atmos_q    = _make_tvi(_t_actual, _q,    _t_spinup)
LW_d       = _make_tvi(_t_actual, _LW,   _t_spinup)
SW_d       = _make_tvi(_t_actual, _SW,   _t_spinup)
c_co2      = _make_co2(_t_actual, _co2,  _t_spinup)

# Separate liquid/snow precip (use phase partitioning from temperature: T < 273.15 → snow)
# Also exclude pre-1997 UTC-offset points (same reason as in _make_tvi above)
_valid_prec     = .!isnan.(_prec) .& .!isnan.(_T) .& (_t97 .>= 0.0)
_valid_prec97   = _valid_prec .& _mask97          # spinup: valid prec in 1997 window

# Full (actual 1997-2014) precip, time-shifted to sim_start basis
_rain_vals  = copy(_prec[_valid_prec])
_snow_vals  = zero(_rain_vals)
mask_snow   = _T[_valid_prec] .< FT(273.15)
_snow_vals[mask_snow] .= _rain_vals[mask_snow]
_rain_vals[mask_snow] .= 0.0

# Spinup (recycled 1997) precip, time-shifted back by SPINUP_SHIFT_DAYS
_rain97 = copy(_prec[_valid_prec97])
_snow97 = zero(_rain97)
mask_snow97 = _T[_valid_prec97] .< FT(273.15)
_snow97[mask_snow97] .= _rain97[mask_snow97]; _rain97[mask_snow97] .= 0.0

# Spinup timestamps: _t_spinup is relative to sim_start, indexed over all of _mask97;
# we need only the subset where precip is also valid
_t_spinup_prec = _t_spinup[_valid_prec97[_mask97]]

atmos_P_liq  = TimeVaryingInput(
    vcat(_t_spinup_prec, _t_actual[_valid_prec]),
    vcat(_rain97, _rain_vals))
atmos_P_snow = TimeVaryingInput(
    vcat(_t_spinup_prec, _t_actual[_valid_prec]),
    vcat(_snow97, _snow_vals))

atmos = ClimaLand.PrescribedAtmosphere(
    atmos_P_liq, atmos_P_snow,
    atmos_T, atmos_u, atmos_q, atmos_P,
    sim_start, FT(atmos_h), toml_dict;
    gustiness = FT(1),
    c_co2 = c_co2,
)
radiation = ClimaLand.PrescribedRadiativeFluxes(FT, SW_d, LW_d, sim_start)

# LAI with recycled 1997 spinup
_valid_lai97  = .!isnan.(_lai_raw) .& _mask97
_valid_lai_all = .!isnan.(_lai_raw) .& (_t97 .>= 0.0)  # exclude pre-1997 UTC-offset points
lai_seconds = vcat(
    _t_spinup[_valid_lai97[_mask97]],
    _t_actual[_valid_lai_all],
)
lai_vals = vcat(_lai_raw[_valid_lai97], _lai_raw[_valid_lai_all])
LAI = TimeVaryingInput(lai_seconds, lai_vals)

# ── Explicit site-specific parameters ────────────────────────────────────────
χl          = FT(0.25)
α_PAR_leaf  = FT(0.1)
α_NIR_leaf  = FT(0.45)
τ_PAR_leaf  = FT(0.05)
τ_NIR_leaf  = FT(0.25)
Ω           = FT(1)
rooting_depth = FT(0.3)
ν           = FT(0.45)
θ_r         = FT(0.07)
K_sat       = FT(1e-5)
vg_n        = FT(1.6)
vg_α        = FT(1.6)

hydrology_cm       = ClimaLand.Soil.vanGenuchten{FT}(; α = vg_α, n = vg_n)
retention_parameters  = (; ν, hydrology_cm, θ_r, K_sat)
composition_parameters = (; ν_ss_om = FT(0.03), ν_ss_quartz = FT(0.47), ν_ss_gravel = FT(0.12))

# ── Build model (mirrors run_dk_sor_default.jl exactly) ──────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing_nt = (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}())

biomass = Canopy.PrescribedBiomassModel{FT}(
    land_domain, LAI, toml_dict;
    rooting_depth,
    height = FT(25), SAI = FT(1.5), RAI = FT(1.5),
)
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf,
)
radiative_transfer = Canopy.TwoStreamModel{FT}(canopy_domain, toml_dict; radiation_parameters)
hydraulics = Canopy.PlantHydraulicsModel{FT}(canopy_domain, toml_dict)
canopy = Canopy.CanopyModel{FT}(
    canopy_domain, forcing_nt, LAI, toml_dict;
    prognostic_land_components,
    photosynthesis = Canopy.PModel{FT}(canopy_domain, toml_dict),
    conductance    = Canopy.PModelConductance{FT}(toml_dict),
    soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict; soil_params = (; ν, θ_r)),
    biomass, radiative_transfer, hydraulics,
)

soil = ClimaLand.Soil.EnergyHydrology{FT}(
    land_domain, (; atmos, radiation), toml_dict;
    prognostic_land_components, retention_parameters, composition_parameters,
    S_s = FT(1e-3),
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
)

land = LandModel{FT}(
    (; atmos, radiation), LAI, toml_dict, land_domain, DT;
    prognostic_land_components, canopy, soil,
)

# ── Initial conditions (same as model_interface.jl) ──────────────────────────
function set_ic!(Y, p, t, model)
    earth_param_set = ClimaLand.get_earth_param_set(model.soil)
    evaluate!(p.drivers.T, atmos.T, t)
    (; θ_r, ν, ρc_ds) = model.soil.parameters
    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) .* FT(0.95)
    Y.soil.θ_i .= FT(0)
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, ρc_ds, earth_param_set)
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc_s, p.drivers.T, earth_param_set)
    Y.snow.S .= FT(0); Y.snow.S_l .= FT(0); Y.snow.U .= FT(0)
    if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
        Y.canopy.energy.T .= p.drivers.T
    end
    Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν
    if !isnothing(model.soilco2)
        Y.soilco2.CO2 .= FT(0.000412); Y.soilco2.O2_f .= FT(0.21)
        SOC_top = FT(15.0); SOC_bot = FT(0.5)
        τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end
end

output_writer = ClimaDiagnostics.Writers.DictWriter()
diags = ClimaLand.default_diagnostics(
    land, sim_start; output_writer, output_vars = :short, reduction_period = :daily)

simulation = LandSimulation(sim_start, sim_stop, DT, land;
    set_ic! = set_ic!, updateat = Second(DT), diagnostics = diags)

println("Running calibration prior-mean model...")
@time solve!(simulation)
println("Done.")

# ── Save raw diagnostics to disk immediately (so crashes in post-processing don't lose data)
using JLD2
const DIAG_CACHE = joinpath(@__DIR__, "calibrate_dk_sor_output",
    SMOKE_TEST ? "prior_mean_full_diag_smoke.jld2" : "prior_mean_full_diag.jld2")
mkpath(dirname(DIAG_CACHE))

function _extract_raw(sim, name)
    writer = nothing
    for d in sim.diagnostics
        haskey(d.output_writer.dict, name) && (writer = d.output_writer; break)
    end
    isnothing(writer) && return nothing, nothing
    times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, name)
    return times, Float64.(data)
end

raw = Dict{String,Any}()
for var in ("nee_1d_average","lhf_1d_average","shf_1d_average","gpp_1d_average","er_1d_average")
    t, d = _extract_raw(simulation, var)
    raw["times"] = t  # same for all vars
    raw[var] = d
end
raw["sim_start"] = sim_start
JLD2.save(DIAG_CACHE, raw)
println("Diagnostics saved to: $DIAG_CACHE")

# ── Extract diagnostics ───────────────────────────────────────────────────────
function get_diag(sim, name)
    writer = nothing
    for d in sim.diagnostics
        haskey(d.output_writer.dict, name) && (writer = d.output_writer; break)
    end
    isnothing(writer) && error("$name not found")
    times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, name)
    # times may be Vector{DateTime} or Vector of seconds-since-sim_start
    dates = Date.(date.(times))
    return dates, Float64.(data)
end

nee_dates, nee_vals = get_diag(simulation, "nee_1d_average")
lhf_dates, lhf_vals = get_diag(simulation, "lhf_1d_average")
shf_dates, shf_vals = get_diag(simulation, "shf_1d_average")
gpp_dates, gpp_vals = get_diag(simulation, "gpp_1d_average")
er_dates, er_vals   = get_diag(simulation, "er_1d_average")

# convert NEE: mol CO2/m2/s → gC/m2/d
nee_gC = nee_vals .* 12.0 .* 86400.0
gpp_gC = gpp_vals .* 12.0 .* 86400.0
er_gC  = er_vals  .* 12.0 .* 86400.0

# filter to full analysis period 1997-2014 (post-spinup)
mask_analysis = (nee_dates .>= Date(1997,1,1)) .& (nee_dates .< Date(2015,1,1))
nee_dates2 = nee_dates[mask_analysis]
nee_gC2    = nee_gC[mask_analysis]
gpp_gC2    = gpp_gC[mask_analysis]
er_gC2     = er_gC[mask_analysis]
lhf2       = lhf_vals[mask_analysis]
shf2       = shf_vals[mask_analysis]

# Load observations for 1997-2013 (obs file only goes through 2013)
flux_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
flux_ds = NCDataset(flux_nc_path, "r")
flux_times = Date.(flux_ds["time"][:])
nee_obs_raw = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
qle_obs_raw = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
qh_obs_raw  = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))
close(flux_ds)

mask_obs = (flux_times .>= Date(1997,1,1)) .& (flux_times .< Date(2014,1,1))
obs_dates = flux_times[mask_obs]
nee_obs   = nee_obs_raw[mask_obs]
qle_obs   = qle_obs_raw[mask_obs]
qh_obs    = qh_obs_raw[mask_obs]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = Figure(size=(1200, 1000))

ax1 = Axis(fig[1,1]; ylabel="NEE (gC/m²/d)", title="DK-Sor 1997–2014 — Prior Mean vs Obs")
lines!(ax1, nee_dates2, nee_gC2; color=:blue, linewidth=1.5, label="Prior mean (calib model)")
lines!(ax1, obs_dates, nee_obs; color=:green, linewidth=1.5, label="Obs")
axislegend(ax1; position=:rt, framevisible=false)

ax2 = Axis(fig[2,1]; ylabel="GPP (gC/m²/d)")
lines!(ax2, nee_dates2, gpp_gC2; color=:blue, linewidth=1.5, label="GPP prior mean")
lines!(ax2, nee_dates2, er_gC2;  color=:red,  linewidth=1.5, label="ER prior mean")
axislegend(ax2; position=:rt, framevisible=false)

ax3 = Axis(fig[3,1]; ylabel="Qle (W/m²)")
lines!(ax3, nee_dates2, lhf2; color=:blue, linewidth=1.5, label="Prior mean")
lines!(ax3, obs_dates, qle_obs; color=:green, linewidth=1.5, label="Obs")
axislegend(ax3; position=:rt, framevisible=false)

ax4 = Axis(fig[4,1]; xlabel="Date", ylabel="Qh (W/m²)")
lines!(ax4, nee_dates2, shf2; color=:blue, linewidth=1.5, label="Prior mean")
lines!(ax4, obs_dates, qh_obs; color=:green, linewidth=1.5, label="Obs")
axislegend(ax4; position=:rt, framevisible=false)

outpath = joinpath(@__DIR__, "calibrate_dk_sor_output", "prior_mean_full_1997_2014.png")
mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")
_stats(v) = isempty(filter(!isnan,v)) ? "all NaN" : "min=$(round(minimum(filter(!isnan,v)), digits=2)), max=$(round(maximum(filter(!isnan,v)), digits=2)), mean=$(round(mean(filter(!isnan,v)), digits=2))"
println("\nNEE stats (1997-2014): $(_stats(nee_gC2)) gC/m²/d")
println("GPP stats (1997-2014): $(_stats(gpp_gC2))")
println("ER stats  (1997-2014): $(_stats(er_gC2))")
