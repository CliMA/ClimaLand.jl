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
ClimaComms.@import_required_backends
using CairoMakie
using Statistics
import StatsBase: percentile
using Dates
using Insolation
using StaticArrays
import ClimaUtilities.OutputPathGenerator: generate_output_path

# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations

const FT = Float64;
earth_param_set = LP.LandParameters(FT);
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
toml_dict = LP.create_toml_dict(FT, default_params_filepath)

# Site-specific information
time_offset = 7 # difference from UTC in hours
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
land_domain = Point(; z_sfc = FT(0.0), longlat = (long, lat))
atmos_h = FT(32)
site_ID = "US-MOz"
start_date = DateTime(2010) + Hour(time_offset) + Day(150)
N_days = 20.0
stop_date = start_date + Day(N_days) + Second(80)

# Get prescribed atmospheric and radiation forcing
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

# Read in LAI from MODIS data
surface_space = land_domain.space.surface
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# Overwrite energy parameter for stability
energy = BigLeafEnergyModel{FT}(; ac_canopy = FT(1e3))

# Construct canopy model
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    land_domain,
    forcing,
    LAI,
    toml_dict;
    energy,
);

(; retention_model, ν, S_s) = canopy.hydraulics.parameters;
ψ_leaf_0 = FT(-2e5 / 9800)
ψ_stem_0 = FT(-1e5 / 9800)
S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        ν,
        S_s,
    )

timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);
function set_ic!(Y, p, t0, model)
    Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini[1])
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
end
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

    # Create update times for driver and saving callback
    saveat = vcat(Array(start_date:Second(3 * 3600):stop_date), stop_date)
    updateat = vcat(Array(start_date:Second(3 * 3600):stop_date), stop_date)
    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        canopy;
        set_ic! = set_ic!,
        updateat = updateat,
        solver_kwargs = (; saveat = deepcopy(saveat)),
        timestepper = ode_algo,
        user_callbacks = (),
    )
    @time sol = solve!(simulation)
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
        push!(times, float.(sol.t))
    end
end
savedir = generate_output_path(
    joinpath("experiments", "standalone", "Vegetation", "timestep_test"),
);

# Create plot with statistics
sim_time = round(Dates.value(Second(stop_date - start_date)) / 3600, digits = 2) # simulation length in hours

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
axislegend(ax, position = :rb)
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
