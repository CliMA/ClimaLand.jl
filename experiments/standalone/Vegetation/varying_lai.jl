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
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaUtilities.OutputPathGenerator: generate_output_path
const FT = Float32;
earth_param_set = LP.LandParameters(FT);

time_offset = 7
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
land_domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
atmos_h = FT(32)
site_ID = "US-MOz"
start_date = DateTime(2010) + Hour(time_offset)
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
    if t < 30 * 24 * 3600
        0.0
    elseif t < (364 - 30) * 24 * 3600.0
        max(2.0 * sin(2 * π / (730 * 24 * 3600) * (t - 30 * 24 * 3600)), 0)
    else
        0.0
    end
end
LAI = TimeVaryingInput(fakeLAIfunction)

# Overwrite some plant hydraulics defaults
n_leaf = 1
h_leaf = FT(9.5)
f_root_to_shoot = FT(3.5)
RAI = FT(2 * f_root_to_shoot)
ν = FT(0.7)
hydraulics = PlantHydraulicsModel{FT}(land_domain, LAI; n_leaf, h_leaf, RAI, ν)

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    land_domain,
    forcing,
    LAI,
    earth_param_set;
    hydraulics,
);

Y, p, coords = ClimaLand.initialize(canopy)
exp_tendency! = make_exp_tendency(canopy);
imp_tendency! = make_imp_tendency(canopy)
jacobian! = make_jacobian(canopy);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);


(; retention_model, ν, S_s) = canopy.hydraulics.parameters;
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

# Set up timestepper
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat);

savedir = generate_output_path("experiments/standalone/Vegetation/varying_lai");
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
