import SciMLBase
using CairoMakie
using Statistics
using Dates
using Insolation

# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using StaticArrays
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float32;
earth_param_set = LP.LandParameters(FT);
f_root_to_shoot = FT(3.5)
SAI = FT(0.0)
plant_ν = FT(2.46e-4) # kg/m^2
n_stem = Int64(0)
n_leaf = Int64(1)
h_leaf = FT(9.5)
compartment_midpoints = [h_leaf / 2]
compartment_surfaces = [FT(0), h_leaf]
land_domain = Point(; z_sfc = FT(0.0))
include(
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/data_tools.jl"),
);
time_offset = 7
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
atmos_h = FT(32)
site_ID = "US-MOz"
data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv"

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
);

z0_m = FT(2)
z0_b = FT(0.2)

shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
);
ψ_soil0 = FT(0.0)

soil_driver = PrescribedGroundConditions(
    FT;
    root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ = t -> ψ_soil0,
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    T = t -> 298.0,
    ϵ = FT(0.99),
);

rt_params = TwoStreamParameters(
    FT;
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.45),
    τ_PAR_leaf = FT(0.05),
    τ_NIR_leaf = FT(0.25),
    Ω = FT(0.69),
    λ_γ_PAR = FT(5e-7),
)

rt_model = TwoStreamModel{FT}(rt_params);

cond_params = MedlynConductanceParameters(FT; g1 = FT(141.0))

stomatal_model = MedlynConductanceModel{FT}(cond_params);

is_c3 = FT(1) # set the photosynthesis mechanism to C3
photo_params = FarquharParameters(FT, is_c3; Vcmax25 = FT(5e-5))

photosynthesis_model = FarquharModel{FT}(photo_params);

AR_params = AutotrophicRespirationParameters(FT)
AR_model = AutotrophicRespirationModel{FT}(AR_params);

function fakeLAIfunction(t)
    if t < 30 * 24 * 3600
        0.0
    elseif t < (364 - 30) * 24 * 3600.0
        max(2.0 * sin(2 * π / (730 * 24 * 3600) * (t - 30 * 24 * 3600)), 0)
    else
        0.0
    end
end


f_root_to_shoot = FT(3.5)
SAI = FT(0)
RAI = FT(2 * f_root_to_shoot)
ai_parameterization =
    PrescribedSiteAreaIndex{FT}(TimeVaryingInput(fakeLAIfunction), SAI, RAI)
rooting_depth = FT(1.0);

K_sat_plant = FT(1.8e-8)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
a = FT(0.05 * 0.0098)

conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)

retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);

ν = FT(0.7)
S_s = FT(1e-2 * 0.0098)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = ν,
    S_s = S_s,
    rooting_depth = rooting_depth,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
);


plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_surfaces = compartment_surfaces,
    compartment_midpoints = compartment_midpoints,
);

energy_model = ClimaLand.Canopy.BigLeafEnergyModel{FT}(
    BigLeafEnergyParameters{FT}(FT(1e4)),
)

canopy = ClimaLand.Canopy.CanopyModel{FT}(;
    parameters = shared_params,
    domain = land_domain,
    autotrophic_respiration = AR_model,
    radiative_transfer = rt_model,
    photosynthesis = photosynthesis_model,
    conductance = stomatal_model,
    energy = energy_model,
    hydraulics = plant_hydraulics,
    boundary_conditions = Canopy.AtmosDrivenCanopyBC(
        atmos,
        radiation,
        soil_driver,
    ),
);


Y, p, coords = ClimaLand.initialize(canopy)
exp_tendency! = make_exp_tendency(canopy);


ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini = inverse_water_retention_curve(retention_model, ψ_leaf_0, ν, S_s)

Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini)



t0 = 0.0
N_days = 364
tf = t0 + 3600 * 24 * N_days
dt = 225.0;
evaluate!(Y.canopy.energy.T, atmos.T, t0)
set_initial_cache! = make_set_initial_cache(canopy)
set_initial_cache!(p, Y, t0);


n = 16
saveat = Array(t0:(n * dt):tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat);

updateat = Array(t0:1800:tf)
drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb);

timestepper = CTS.RK4();
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
);

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat);

savedir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Vegetation");
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
    parent(sv.saveval[k].canopy.energy.turbulent_fluxes.shf)[1] for
    k in 1:length(sol.t)
]
LHF = [
    parent(sv.saveval[k].canopy.energy.turbulent_fluxes.lhf)[1] for
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
    parent(sv.saveval[k].canopy.energy.turbulent_fluxes.transpiration)[1]
    for k in 1:length(sol.t)
]
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Temperature (K)")
lines!(ax, sol.t ./ 24 ./ 3600, T, label = "Canopy")
lines!(ax, sol.t ./ 24 ./ 3600, T_atmos, label = "Atmos")
axislegend(ax)
ax = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "Volumetric Water")
lines!(ax, sol.t ./ 24 ./ 3600, ϑ, label = "Canopy")
axislegend(ax)
ax = Axis(fig[3, 1], xlabel = "Time (days)", ylabel = "LAI")
lines!(ax, sol.t ./ 24 ./ 3600, fakeLAIfunction.(sol.t), label = "Canopy")
axislegend(ax)
save(joinpath(savedir, "varying_lai_state.png"), fig)
fig2 = Figure()
ax = Axis(fig2[1, 1], xlabel = "Time (days)", ylabel = "Energy Fluxes")
lines!(ax, sol.t ./ 24 ./ 3600, SW_n, label = "SW_n")
lines!(ax, sol.t ./ 24 ./ 3600, LW_n, label = "LW_n")
lines!(ax, sol.t ./ 24 ./ 3600, SHF, label = "SHF")
lines!(ax, sol.t ./ 24 ./ 3600, LHF, label = "LHF")
lines!(ax, sol.t ./ 24 ./ 3600, RE, label = "RE")
axislegend(ax)
ax = Axis(fig2[2, 1], xlabel = "Time (days)", ylabel = "Water Fluxes")
lines!(ax, sol.t ./ 24 ./ 3600, Tr, label = "Transpiration")
lines!(ax, sol.t ./ 24 ./ 3600, R, label = "R")

axislegend(ax)
ax = Axis(fig2[3, 1], xlabel = "Time (days)", ylabel = "Carbon Fluxes")
lines!(ax, sol.t ./ 24 ./ 3600, GPP, label = "GPP")
lines!(ax, sol.t ./ 24 ./ 3600, resp, label = "Respiration")
axislegend(ax)
save(joinpath(savedir, "varying_lai_fluxes.png"), fig2)
