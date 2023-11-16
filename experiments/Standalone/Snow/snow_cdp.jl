import SciMLBase
import ClimaTimeSteppers as CTS
using Plots
using ClimaCore
using ClimaLSM.Snow
using ClimaLSM.Domains
import ClimaLSM
import ClimaLSM.Parameters as LSMP
climalsm_dir = pkgdir(ClimaLSM)
include(joinpath(climalsm_dir, "parameters", "create_parameters.jl"))
# Debug type promotion when we use Float32
FT = Float64
param_set = create_lsm_parameters(FT)
# This reads in the data and sets up the drivers, as well as computes the IC from the data
include(joinpath(climalsm_dir, "experiments/Standalone/Snow/process_cdp.jl"))
savedir = joinpath(climalsm_dir, "experiments/Standalone/Snow/")

domain = Point(; z_sfc = FT(0))

model = ClimaLSM.Snow.SnowModel(
    parameters = parameters,
    domain = domain,
    atmosphere = atmos,
    radiation = rad,
)
Y, p, coords = ClimaLSM.initialize(model)

# Set initial conditions
Y.snow.S .= S_0
Y.snow.U .= U_0

update_aux! = ClimaLSM.make_update_aux(model)
t0 = FT(0.0)
ndays = 180
tf = FT(ndays * 3600 * 24)
#tf = FT(seconds[end])
update_aux!(p, Y, t0)
exp_tendency! = ClimaLSM.make_exp_tendency(model)
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLSM.dss!),
    Y,
    (t0, tf),
    p,
)
saveat = FT.(collect(t0:(24 * 3600):tf))
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat);

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
evaporation =
    [parent(sv.saveval[k].snow.evaporation)[1] for k in 1:length(sol.t)];
turbulent_energy_flux = [
    parent(sv.saveval[k].snow.turbulent_energy_flux)[1] for k in 1:length(sol.t)
];
R_n = [parent(sv.saveval[k].snow.R_n)[1] for k in 1:length(sol.t)];
energy_runoff =
    [parent(sv.saveval[k].snow.energy_runoff)[1] for k in 1:length(sol.t)];
water_runoff =
    [parent(sv.saveval[k].snow.water_runoff)[1] for k in 1:length(sol.t)];
total_water_flux =
    [parent(sv.saveval[k].snow.total_water_flux)[1] for k in 1:length(sol.t)];
total_energy_flux =
    [parent(sv.saveval[k].snow.total_energy_flux)[1] for k in 1:length(sol.t)];
snow_water_flux =
    [parent(sv.saveval[k].snow.snow_water_flux)[1] for k in 1:length(sol.t)];
snow_energy_flux =
    [parent(sv.saveval[k].snow.snow_energy_flux)[1] for k in 1:length(sol.t)];



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
plot1b = plot(xticks = :none, ylabel = "Flux (inches/day)", xlim = (0, ndays))
plot!(
    plot1b,
    t ./ 3600 ./ 24,
    3600 * 24 * 39.3701 .* model.atmos.snow_precip.(t),
    label = "Snow",
    color = "red",
)
plot!(
    plot1b,
    t ./ 3600 ./ 24,
    3600 * 24 * 39.3701 .* model.atmos.liquid_precip.(t),
    label = "Rain",
    color = "green",
)
plot!(
    plot1b,
    t ./ 3600 ./ 24,
    3600 * 24 * 39.3701 .* evaporation,
    label = "Evaporation",
    color = "purple",
)
plot!(
    plot1b,
    t ./ 3600 ./ 24,
    3600 * 24 * 39.3701 .* water_runoff,
    label = "Runoff",
    color = "blue",
)
plot!(plot1b, ylim = [-4, 0.4])

plot1c = plot(xlim = (0, ndays))
plot!(
    plot1c,
    t ./ 3600 ./ 24,
    q_l,
    xlabel = "Time (days)",
    label = "",
    ylabel = "Liquid mass fraction",
)

plot(plot1a, plot1b, plot1c, layout = (3, 1), size = (800, 800))
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
    xticks = :none,
    ylabel = "Snow Energy per Area",
)

plot2c = plot(xlim = (0, ndays), ylim = (-500, 500))
plot!(
    plot2c,
    t ./ 3600 ./ 24,
    R_n,
    label = "Net Rad",
    xlabel = "Time (days)",
    ylabel = "Fluxes (W/m^2)",
)
plot!(plot2c, t ./ 3600 ./ 24, turbulent_energy_flux, label = "SHF+LHF")
plot!(plot2c, t ./ 3600 ./ 24, energy_runoff, label = "Energy Runoff")

plot(plot2a, plot2b, plot2c, layout = (3, 1), size = (800, 800))
savefig(joinpath(savedir, "snow_energy_content.png"))
