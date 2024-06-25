import SciMLBase
import ClimaTimeSteppers as CTS
using Plots
import ClimaParams as CP
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
using ClimaCore
using ClimaLand.Snow
using ClimaLand.Domains
using ClimaComms
import ClimaLand
import ClimaLand.Parameters as LP
using DataFrames
using NCDatasets
using Dates
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
using Thermodynamics
using Statistics
using Insolation
using DelimitedFiles

# Site-specific quantities
# Error if no site argument is provided
if length(ARGS) < 1
    @error("Please provide a site name as command line argument")
else
    SITE_NAME = ARGS[1]
end

climaland_dir = pkgdir(ClimaLand)

FT = Float32
param_set = LP.LandParameters(FT)
context = ClimaComms.context()

# This reads in the data and sets up the drivers, as well as computes the IC from the data
include(
    joinpath(climaland_dir, "experiments/standalone/Snow/process_snowmip.jl"),
)
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
    atmos = atmos,
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
saveat = FT.(collect(t0:3600:tf))
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

start_day = 1
days = start_day .+ floor.(t ./ 3600 ./ 24)
doys = days .% 365 # doesn't account for leap year

obs_swes = Vector{Union{Float64, Missing}}(missing, length(doys))
obs_swes[snow_data_avail] .= mass[snow_data_avail] ./ 1000

obs_tsnows = Vector{Union{Float64, Missing}}(missing, length(doys))
obs_tsnows[snow_data_avail] = T_snow .+ 273.15

obs_df = DataFrame(
    doy = doys,
    model_swe = S,
    obs_swe = obs_swes,
    model_tsnow = T,
    obs_tsnow = obs_tsnows,
)
function missingmean(x)
    return mean(skipmissing(x))
end

mean_obs_df = combine(
    groupby(obs_df, :doy),
    [:model_swe, :obs_swe, :model_tsnow, :obs_tsnow] .=> missingmean,
    renamecols = false,
)

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

plot1c = plot()
plot!(
    plot1c,
    mean_obs_df.doy,
    mean_obs_df.model_swe,
    label = "Model",
    xlabel = "Day of year",
    ylabel = "Average snow water",
)
scatter!(plot1c, mean_obs_df.doy, mean_obs_df.obs_swe, label = "Data")

plot(plot1a, plot1b, plot1c, layout = (3, 1), size = (800, 800))
savefig(joinpath(savedir, "snow_water_content_$(SITE_NAME).png"))

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
# Plot the observed atmospheric temperature
scatter!(
    plot2a,
    seconds[snow_data_avail] ./ 3600 ./ 24,
    Tair[snow_data_avail],
    label = "Atmospheric Temperature",
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

plot2c = plot()
plot!(
    plot2c,
    mean_obs_df.doy,
    mean_obs_df.model_tsnow,
    label = "Model",
    xlabel = "Day of year",
    ylabel = "Average snow temperature",
)
scatter!(plot2c, mean_obs_df.doy, mean_obs_df.obs_tsnow, label = "Data")

snow_energy_plot =
    plot(plot2a, plot2b, plot2c, layout = (3, 1), size = (800, 800))
savefig(joinpath(savedir, "snow_energy_content_$(SITE_NAME).png"))

plot3 = plot(
    xlabel = "Time (days)",
    ylabel = "Cumulative height (m)",
    xlim = (0, ndays),
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
savefig(joinpath(savedir, "water_fluxes_$(SITE_NAME).png"))
