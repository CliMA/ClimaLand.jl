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
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaUtilities.OutputPathGenerator: generate_output_path
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
N_days = 60
stop_date = start_date + Day(N_days)
dt = FT(225);
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

function fakeLAIfunction2(t)
    if FT(t) < 10 * 24 * 3600
        0.0
    elseif FT(t) < (60 - 10) * 24 * 3600.0
        max(2.0 * sin(2 * π / (730 * 24 * 3600) * (FT(t) - 10 * 24 * 3600)), 0)
    else
        0.0
    end
end
LAI = TimeVaryingInput(fakeLAIfunction2)
f_root_to_shoot = FT(3.5)
SAI = FT(1.0)
RAI = FT(3 * f_root_to_shoot)
ν = FT(0.7)
S_s = FT(1e-2 * 0.0098)
rooting_depth = FT(1.0);

K_sat_plant = FT(1.8e-6)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
a = FT(0.05 * 0.0098)
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)

retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ν = ν,
    S_s = S_s,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
);

n_stem = Int64(1)
n_leaf = Int64(1)
h_leaf = FT(9.5)
h_stem = FT(9)
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [FT(0), h_stem, h_stem + h_leaf]
plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_surfaces = compartment_surfaces,
    compartment_midpoints = compartment_midpoints,
);
height = h_leaf + h_stem
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    land_domain,
    forcing,
    LAI,
    toml_dict;
    hydraulics = plant_hydraulics,
    biomass,
);


function set_ic!(Y, p, t0, model)
    atmos = model.boundary_conditions.atmos
    (; retention_model, ν, S_s) = model.hydraulics.parameters
    ψ_leaf_0 = FT(-2e5 / 9800)
    ψ_stem_0 = FT(-1e5 / 9800)

    S_l_ini =
        inverse_water_retention_curve.(
            retention_model,
            [ψ_stem_0, ψ_leaf_0],
            ν,
            S_s,
        )

    Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini[1])
    Y.canopy.hydraulics.ϑ_l.:2 .= augmented_liquid_fraction.(ν, S_l_ini[2])
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
savedir = generate_output_path(
    "experiments/standalone/Vegetation/varying_lai_with_stem",
);
T = [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]
T_atmos = [parent(sv.saveval[k].drivers.T)[1] for k in 1:length(sol.t)]
ϑ_l = [parent(sol.u[k].canopy.hydraulics.ϑ_l.:2)[1] for k in 1:length(sol.t)]
ϑ_s = [parent(sol.u[k].canopy.hydraulics.ϑ_l.:1)[1] for k in 1:length(sol.t)]
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
R_stem_leaf =
    [parent(sv.saveval[k].canopy.hydraulics.fa.:1)[1] for k in 1:length(sol.t)]
Tr = [
    parent(sv.saveval[k].canopy.turbulent_fluxes.transpiration)[1] for
    k in 1:length(sol.t)
]
times = float.(sol.t) ./ 24 ./ 3600
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Temperature (K)")
lines!(ax, times, T, label = "Canopy")
lines!(ax, times, T_atmos, label = "Atmos")
axislegend(ax)
ax = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "Volumetric Water")
lines!(ax, times, ϑ_l, label = "Leaf")
lines!(ax, times, ϑ_s, label = "Stem")
axislegend(ax)
ax = Axis(fig[3, 1], xlabel = "Time (days)", ylabel = "LAI")
lines!(ax, times, fakeLAIfunction2.(sol.t), label = "Canopy")
axislegend(ax)
save(joinpath(savedir, "varying_lai_with_stem_state.png"), fig)
fig2 = Figure()
ax = Axis(fig2[1, 1], xlabel = "Time (days)", ylabel = "Energy Fluxes")
lines!(ax, times, SW_n, label = "SW_n")
lines!(ax, times, LW_n, label = "LW_n")
lines!(ax, times, SHF, label = "SHF")
lines!(ax, times, LHF, label = "LHF")
lines!(ax, times, RE, label = "RE")
axislegend(ax)
ax = Axis(fig2[2, 1], xlabel = "Time (days)", ylabel = "Water Fluxes")
ylims!(ax, (0, 1e-7))
lines!(ax, times, Tr, label = "Transpiration")
lines!(ax, times, R, label = "R")
lines!(ax, times, R_stem_leaf, label = "R_stem_leaf")

axislegend(ax)
ax = Axis(fig2[3, 1], xlabel = "Time (days)", ylabel = "Carbon Fluxes")
lines!(ax, times, GPP, label = "GPP")
lines!(ax, times, resp, label = "Respiration")
axislegend(ax)
save(joinpath(savedir, "varying_lai_with_stem_fluxes.png"), fig2)
