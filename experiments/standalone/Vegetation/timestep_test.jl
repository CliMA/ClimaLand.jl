# # Timestep test for standalone CanopyModel

# This experiment runs the standalone canopy model for 20 days at a single point
# and outputs useful plots related to the timestepping and stability of
# the simulation using different timesteps. This was initially implemented to
# test the behavior of stepping canopy temperature (T) implicitly.

# The simulation is run multiple times: first, with a small reference timestep,
# then with increasing timestep sizes. The goal here is to see which timesteps
# are small enough to produce a stable simulation, and which timesteps are too
# large to produce accurate results.

# The runs with larger timesteps are compared with the reference timestep run
# using multiple metrics, including:
# - mean, 99%th, and 95th% error over the course of the simulation
# - T at the end of the simulation (used to test convergence with each dt)
# - T throughout the entire simulation

# Note that the simulation is run for 20 days + 80 seconds. This somewhat
# unusual length is necessary to test the convergence behavior of the solver.
# Convergence behavior must be assessed when T is actively changing, rather
# than when it is in equilibrium. T changes in the time period after a driver
# update, which occurs every 3 hours (so exactly at 20 days). Running for
# 20 days + 80 seconds allows us to look at results over a longer simulation,
# as well as convergence behavior. See the discussion in
# https://github.com/CliMA/ClimaLand.jl/pull/675 for more information.

# Simulation Setup
# Space: single point
# Simulation duration: 20 days + 80 seconds (= 480.02 hours)
# Timestep: variable, ranging from 0.6s to 3600s
# Timestepper: ARS111
# Fixed number of iterations: 6
# Jacobian update: Every Newton iteration
# Atmospheric and radiation driver updates: Every 3 hours

import SciMLBase
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using CairoMakie
using Statistics
import StatsBase: percentile
using Dates
using Insolation
using StaticArrays

# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float64;
earth_param_set = LP.LandParameters(FT);
f_root_to_shoot = FT(3.5)
plant_ν = FT(2.46e-4) # kg/m^2
n_stem = Int64(1)
n_leaf = Int64(1)
h_leaf = FT(9.5)
h_stem = FT(9)
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [FT(0), h_stem, h_stem + h_leaf]
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

is_c3 = FT(1)
photo_params = FarquharParameters(FT, is_c3; Vcmax25 = FT(5e-5))
photosynthesis_model = FarquharModel{FT}(photo_params);

AR_params = AutotrophicRespirationParameters(FT)
AR_model = AutotrophicRespirationModel{FT}(AR_params);

f_root_to_shoot = FT(3.5)
SAI = FT(1.0)
RAI = FT(3f_root_to_shoot)
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

K_sat_plant = FT(1.8e-6)
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
    rooting_depth = FT(0.5),
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
ac_canopy = FT(1e3)
energy_model = ClimaLand.Canopy.BigLeafEnergyModel{FT}(
    BigLeafEnergyParameters{FT}(ac_canopy),
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

exp_tendency! = make_exp_tendency(canopy)
imp_tendency! = make_imp_tendency(canopy)
jacobian! = make_jacobian(canopy)
drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(drivers)

ψ_leaf_0 = FT(-2e5 / 9800)
ψ_stem_0 = FT(-1e5 / 9800)
S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        ν,
        S_s,
    )

seconds_per_day = 3600 * 24.0
t0 = 150seconds_per_day
N_days = 20.0
tf = t0 + N_days * seconds_per_day + 80
sim_time = round((tf - t0) / 3600, digits = 2) # simulation length in hours
set_initial_cache! = make_set_initial_cache(canopy)

timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

ref_dt = 6.0
dts = [ref_dt, 12.0, 48.0, 225.0, 450.0, 900.0, 1800.0, 3600.0]

ref_T = []
mean_err = []
p95_err = []
p99_err = []
sol_err = []
T_states = []
times = []
ΔT = FT(0)
for dt in dts
    @info dt

    # Initialize model before each simulation
    Y, p, coords = ClimaLand.initialize(canopy)
    jac_kwargs = (;
        jac_prototype = ClimaLand.ImplicitEquationJacobian(Y),
        Wfact = jacobian!,
    )

    # Set initial conditions
    Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini[1])
    Y.canopy.hydraulics.ϑ_l.:2 .= augmented_liquid_fraction.(ν, S_l_ini[2])
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
    set_initial_cache!(p, Y, t0)

    # Create callback for driver updates
    saveat = vcat(Array(t0:(3 * 3600):tf), tf)
    updateat = vcat(Array(t0:(3 * 3600):tf), tf)
    cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    # Solve simulation
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    @time sol =
        SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)

    # Save results for comparison
    if dt == ref_dt
        global ref_T =
            [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]
    else
        T = [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]

        global ΔT = abs.(T .- ref_T)
        push!(mean_err, FT(mean(ΔT)))
        push!(p95_err, FT(percentile(ΔT, 95)))
        push!(p99_err, FT(percentile(ΔT, 99)))
        push!(sol_err, ΔT[end])
        push!(T_states, T)
        push!(times, sol.t)
    end
end

savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments",
    "standalone",
    "Vegetation",
    "timestep_test",
);
!ispath(savedir) && mkpath(savedir)

# Create plot with statistics
# Compare T state computed with small vs large dt
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "Timestep (minutes)",
    ylabel = "Temperature (K)",
    xscale = log10,
    yscale = log10,
    title = "Error in T over $(sim_time / 24) day sim, dts in [$(dts[1]), $(dts[end])]",
)
dts = dts[2:end] ./ 60
lines!(ax, dts, FT.(mean_err), label = "Mean Error")
lines!(ax, dts, FT.(p95_err), label = "95th% Error")
lines!(ax, dts, FT.(p99_err), label = "99th% Error")
axislegend(ax)
save(joinpath(savedir, "errors.png"), fig)

# Create convergence plot
fig2 = Figure()
ax2 = Axis(
    fig2[1, 1],
    xlabel = "log(dt)",
    ylabel = "log(|T[end] - T_ref[end]|)",
    xscale = log10,
    yscale = log10,
    title = "Convergence of T; sim time = $(sim_time / 24) days, dts in [$(dts[1]), $(dts[end])]",
)
scatter!(ax2, dts, FT.(sol_err))
lines!(ax2, dts, dts / 1e3)
save(joinpath(savedir, "convergence.png"), fig2)

# Plot T throughout full simulation for each run
fig3 = Figure()
ax3 = Axis(
    fig3[1, 1],
    xlabel = "time (hr)",
    ylabel = "T (K)",
    title = "T throughout simulation; length = $(sim_time / 24) days, dts in [$(dts[1]), $(dts[end])]",
)
times = times ./ 3600.0 # hours
for i in 1:length(times)
    lines!(ax3, times[i], T_states[i], label = "dt $(dts[i]) s")
end
axislegend(ax3, position = :rt)
save(joinpath(savedir, "states.png"), fig3)
