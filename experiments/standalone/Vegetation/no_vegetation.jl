import SciMLBase
using CairoMakie
using Statistics
using Dates
using Insolation

# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
using StaticArrays
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ClimaLand.Simulations: LandSimulation, solve!
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

const FT = Float32;
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict);

time_offset = -6
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
land_domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
atmos_h = FT(32)
site_ID = "US-MOz"
start_date = DateTime(2010, 1, 2)
N_days = 10
stop_date = start_date + Day(N_days)
dt = 225.0;

# Prescribed forcing from Fluxnet data
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
)
ground = PrescribedGroundConditions{FT}(;
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    ϵ = FT(0.99),
);
forcing = (; atmos, radiation, ground);

function fakeLAIfunction(t)
    if FT(t) < 30 * 24 * 3600
        0.0
    elseif FT(t) < (364 - 30) * 24 * 3600.0
        max(2.0 * sin(2 * π / (730 * 24 * 3600) * (FT(t) - 30 * 24 * 3600)), 0)
    else
        0.0
    end
end

LAI = TimeVaryingInput(fakeLAIfunction)
SAI = RAI = FT(0)
hydraulics = PlantHydraulicsModel{FT}(land_domain, LAI, toml_dict; SAI, RAI);
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    land_domain,
    forcing,
    LAI,
    toml_dict;
    hydraulics,
)

function set_ic!(Y, p, t0, model)
    atmos = model.boundary_conditions.atmos
    ψ_leaf_0 = FT(-2e5 / 9800)
    (; retention_model, ν, S_s) = model.hydraulics.parameters
    S_l_ini = inverse_water_retention_curve(retention_model, ψ_leaf_0, ν, S_s)
    Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini)
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
end

n = 16
saveat = Second(n * dt)
saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, saveat);
sv = saving_cb.affect!.saved_values;

simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    canopy;
    set_ic! = set_ic!,
    user_callbacks = (saving_cb,),
    updateat = Second(1800),
    solver_kwargs = (; saveat),
    diagnostics = nothing,
)
sol = solve!(simulation)

savedir =
    generate_output_path("experiments/standalone/Vegetation/no_vegetation");
T = [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]
T_atmos = [parent(sv.saveval[k].drivers.T)[1] for k in 1:length(sol.t)]
ϑ = [parent(sol.u[k].canopy.hydraulics.ϑ_l.:1)[1] for k in 1:length(sol.t)]
GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] * 1e6 for
    k in 1:length(sol.t)
]
resp = [
    parent(sv.saveval[k].canopy.autotrophic_respiration.Ra)[1] * 1e6 for
    k in 1:length(sol.t)
]
SW_n = [
    parent(sv.saveval[k].canopy.radiative_transfer.SW_n)[1] for
    k in 1:length(sol.t)
]
LW_n = [
    parent(sv.saveval[k].canopy.radiative_transfer.LW_n)[1] for
    k in 1:length(sol.t)
]
SHF = [
    parent(sv.saveval[k].canopy.turbulent_fluxes.shf)[1] for
    k in 1:length(sol.t)
]
LHF = [
    parent(sv.saveval[k].canopy.turbulent_fluxes.lhf)[1] for
    k in 1:length(sol.t)
]
RE = [
    parent(sv.saveval[k].canopy.energy.fa_energy_roots)[1] for
    k in 1:length(sol.t)
]
R = [
    parent(sv.saveval[k].canopy.hydraulics.fa_roots)[1] for k in 1:length(sol.t)
]
Tr = [
    parent(sv.saveval[k].canopy.turbulent_fluxes.transpiration)[1] for
    k in 1:length(sol.t)
]
times = FT.(sol.t) ./ 24 ./ 3600
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Temperature (K)")
lines!(ax, times, T, label = "Canopy")
lines!(ax, times, T_atmos, label = "Atmos")
axislegend(ax)
ax = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "Volumetric Water")
lines!(ax, times, ϑ, label = "Canopy")
axislegend(ax)
ax = Axis(fig[3, 1], xlabel = "Time (days)", ylabel = "LAI")
lines!(ax, times, fakeLAIfunction.(sol.t), label = "Canopy")
axislegend(ax)
save(joinpath(savedir, "no_veg_state.png"), fig)
fig2 = Figure()
ax = Axis(fig2[1, 1], xlabel = "Time (days)", ylabel = "Energy Fluxes")
lines!(ax, times, SW_n, label = "SW")
lines!(ax, times, LW_n, label = "LW")
lines!(ax, times, SHF, label = "SHF")
lines!(ax, times, LHF, label = "LHF")
lines!(ax, times, RE, label = "RE")
axislegend(ax)
ax = Axis(fig2[2, 1], xlabel = "Time (days)", ylabel = "Water Fluxes")
lines!(ax, times, Tr, label = "Transpiration")
lines!(ax, times, R, label = "R")

axislegend(ax)
ax = Axis(fig2[3, 1], xlabel = "Time (days)", ylabel = "Carbon Fluxes")
lines!(ax, times, GPP, label = "GPP")
lines!(ax, times, resp, label = "Respiration")
axislegend(ax)
save(joinpath(savedir, "no_veg_fluxes.png"), fig2)
