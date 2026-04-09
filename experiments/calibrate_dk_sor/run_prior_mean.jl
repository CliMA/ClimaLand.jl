"""
Run the calibration forward model at the prior mean parameters for 2004,
to compare against experiments/integrated/fluxnet/run_dk_sor_default.jl.

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
using ClimaLand.Canopy.PlantHydraulics
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

const FT = Float64
const climaland_dir = abspath(joinpath(@__DIR__, "..", ".."))
const SITE_ID = "DK-Sor"
const DT = Float64(900)

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

sim_start = DateTime(2003, 1, 1)
sim_stop  = DateTime(2005, 1, 1)
spinup_date = DateTime(2004, 1, 1)
println("Simulating $sim_start → $sim_stop (spinup until $spinup_date)")

toml_dict = LP.create_toml_dict(FT; override_files = [PRIOR_TOML])

(; dz_tuple, nelements, zmin, zmax) = FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long)          = FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h)                         = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

land_domain   = Column(; zlim = (zmin, zmax), nelements, dz_tuple, longlat = (long, lat))
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

met_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path, lat, long, time_offset, atmos_h, sim_start, toml_dict, FT)

# FLUXNET LAI — same as run_dk_sor_default.jl
surface_space = canopy_domain.space.surface
met_ds = NCDataset(met_nc_path, "r")
lai_data = Float64.(coalesce.(met_ds["LAI"][1, 1, :], NaN))
lai_times = met_ds["time"][:]
close(met_ds)
lai_seconds = [Float64(Second(t - Hour(time_offset) - sim_start).value) for t in lai_times]
valid_lai = .!isnan.(lai_data)
LAI = TimeVaryingInput(lai_seconds[valid_lai], lai_data[valid_lai])

# ── Explicit site-specific parameters — same as run_dk_sor_default.jl ────────
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
hydraulics = Canopy.PlantHydraulicsModel{FT}(canopy_domain, toml_dict; n_stem=1, h_stem=FT(24), h_leaf=FT(1))
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
    for i in 1:(model.canopy.hydraulics.n_stem + model.canopy.hydraulics.n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .= model.canopy.hydraulics.parameters.ν
    end
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

# ── Extract diagnostics ───────────────────────────────────────────────────────
function get_diag(sim, name)
    writer = nothing
    for d in sim.diagnostics
        haskey(d.output_writer.dict, name) && (writer = d.output_writer; break)
    end
    isnothing(writer) && error("$name not found")
    times, data = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, name)
    dates = Date.(times isa Vector{DateTime} ? times : date.(times))
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

# filter to 2004 post-spinup
mask2004 = (nee_dates .>= Date(2004,1,1)) .& (nee_dates .< Date(2005,1,1))
nee_dates2 = nee_dates[mask2004]
nee_gC2    = nee_gC[mask2004]
gpp_gC2    = gpp_gC[mask2004]
er_gC2     = er_gC[mask2004]
lhf2       = lhf_vals[mask2004]
shf2       = shf_vals[mask2004]

# Load observations for 2004
flux_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
flux_ds = NCDataset(flux_nc_path, "r")
flux_times = Date.(flux_ds["time"][:])
nee_obs_raw = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
qle_obs_raw = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
qh_obs_raw  = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))
close(flux_ds)

mask_obs = (flux_times .>= Date(2004,1,1)) .& (flux_times .< Date(2005,1,1))
obs_dates = flux_times[mask_obs]
nee_obs   = nee_obs_raw[mask_obs]
qle_obs   = qle_obs_raw[mask_obs]
qh_obs    = qh_obs_raw[mask_obs]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = Figure(size=(1200, 1000))

ax1 = Axis(fig[1,1]; ylabel="NEE (gC/m²/d)", title="DK-Sor 2004 — Prior Mean vs Default vs Obs")
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

outpath = joinpath(@__DIR__, "calibrate_dk_sor_output", "prior_mean_2004.png")
mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")
println("\nNEE stats (2004): min=$(round(minimum(nee_gC2), digits=2)), max=$(round(maximum(nee_gC2), digits=2)), mean=$(round(mean(nee_gC2), digits=2)) gC/m²/d")
println("GPP stats (2004): min=$(round(minimum(gpp_gC2), digits=2)), max=$(round(maximum(gpp_gC2), digits=2))")
println("ER stats  (2004): min=$(round(minimum(er_gC2), digits=2)), max=$(round(maximum(er_gC2), digits=2))")
