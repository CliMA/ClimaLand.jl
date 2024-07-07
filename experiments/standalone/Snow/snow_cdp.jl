import SciMLBase
import ClimaTimeSteppers as CTS
using Plots
import ClimaParams as CP
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaCore
using ClimaLand.Snow
using ClimaLand.Domains
using ClimaComms
import ClimaLand
import ClimaLand.Parameters as LP
climaland_dir = pkgdir(ClimaLand)

FT = Float32
param_set = LP.LandParameters(FT)
context = ClimaComms.context()
# This reads in the data and sets up the drivers, as well as computes the IC from the data
include(joinpath(climaland_dir, "experiments/standalone/Snow/process_cdp.jl"))
savedir = joinpath(climaland_dir, "experiments/standalone/Snow/")
t0 = FT(0.0)
tf = FT(seconds[end])
ndays = (tf - t0) / 3600 / 24
Δt = FT(60 * 60)

domain = Point(; z_sfc = FT(0))

parameters =
    SnowParameters{FT}(Δt; α_snow = α, ρ_snow = ρ, earth_param_set = param_set)
model = ClimaLand.Snow.SnowModel(
    parameters = parameters,
    domain = domain,
    atmosphere = atmos,
    radiation = radiation,
)
Y, p, coords = ClimaLand.initialize(model)

# Set initial conditions
Y.snow.S .= FT(SWE[1]) # first data point
Y.snow.U .=
    ClimaLand.Snow.energy_from_q_l_and_swe(FT(SWE[1]), FT(0), parameters) # with q_l = 0

set_initial_cache! = ClimaLand.make_set_initial_cache(model)
set_initial_cache!(p, Y, t0)
exp_tendency! = ClimaLand.make_exp_tendency(model)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
)
saveat = FT.(collect(t0:Δt:tf))
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat);
updateat = copy(saveat)
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)

sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    adaptive = false,
    saveat = saveat,
    callback = cb,
);

# Plotting
q_l = [parent(sv.saveval[k].snow.q_l)[1] for k in 1:length(sol.t)];
T = [parent(sv.saveval[k].snow.T)[1] for k in 1:length(sol.t)];
evaporation = [
    parent(sv.saveval[k].snow.turbulent_fluxes.vapor_flux)[1] for
    k in 1:length(sol.t)
];
R_n = [parent(sv.saveval[k].snow.R_n)[1] for k in 1:length(sol.t)];
water_runoff =
    [parent(sv.saveval[k].snow.water_runoff)[1] for k in 1:length(sol.t)];
rain = [parent(sv.saveval[k].drivers.P_liq)[1] for k in 1:length(sv.t)];
snow = [parent(sv.saveval[k].drivers.P_snow)[1] for k in 1:length(sv.t)];


S = [parent(sol.u[k].snow.S)[1] for k in 1:length(sol.t)];
U = [parent(sol.u[k].snow.U)[1] for k in 1:length(sol.t)];
t = sol.t;

plot1a = plot()
plot!(
    plot1a,
    t ./ 3600 ./ 24,
    S,
    label = "Model",
    legend = :topright,
    xticks = :none,
    ylabel = "Snow water",
    xlim = (0, ndays),
)
scatter!(plot1a, seconds[snow_data_avail] ./ 24 ./ 3600, SWE, label = "Data")

plot1b = plot(xlim = (0, ndays))
plot!(
    plot1b,
    t ./ 3600 ./ 24,
    q_l,
    xlabel = "Time (days)",
    label = "",
    ylabel = "Liquid mass fraction",
)

plot(plot1a, plot1b, layout = (2, 1), size = (800, 800))
savefig(joinpath(savedir, "snow_water_content.png"))



plot2a = plot(xlim = (0, ndays), ylim = (250, 280))
plot!(
    plot2a,
    t ./ 3600 ./ 24,
    T,
    label = "Model",
    legend = :bottomright,
    xticks = :none,
    ylabel = "Snow Temperature",
)
scatter!(
    plot2a,
    seconds[snow_data_avail] ./ 3600 ./ 24,
    T_snow .+ 273.15,
    label = "Data",
)

plot2b = plot(xlim = (0, ndays))
plot!(
    plot2b,
    t ./ 3600 ./ 24,
    U,
    label = "",
    xlabel = "Time (days)",
    ylabel = "Snow Energy per Area",
)



plot(plot2a, plot2b, layout = (2, 1), size = (800, 800))
savefig(joinpath(savedir, "snow_energy_content.png"))

plot3 = plot(
    xlabel = "Time (days)",
    ylabel = "Cumulative height (m)",
    xlim = (0, ndays),
    title = "Water Fluxes",
)
plot!(plot3, t ./ 3600 ./ 24, cumsum(snow) .* Δt, label = "Snow", color = "red")
plot!(
    plot3,
    t ./ 3600 ./ 24,
    cumsum(rain) .* Δt,
    label = "Rain",
    color = "green",
)
plot!(
    plot3,
    t ./ 3600 ./ 24,
    cumsum(evaporation) .* Δt,
    label = "Evaporation",
    color = "purple",
)
plot!(
    plot3,
    t ./ 3600 ./ 24,
    cumsum(water_runoff) .* Δt,
    label = "Runoff",
    color = "blue",
)
savefig(joinpath(savedir, "water_fluxes.png"))
