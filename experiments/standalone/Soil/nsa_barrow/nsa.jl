import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation
using StatsBase

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
climalsm_dir = pkgdir(ClimaLSM)
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/standalone/Soil/nsa_barrow/met_drivers.jl",
    ),
)
include(
    joinpath(climalsm_dir, "experiments/standalone/Soil/nsa_barrow/domain.jl"),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/standalone/Soil/nsa_barrow/parameters.jl",
    ),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/standalone/Soil/nsa_barrow/simulation.jl",
    ),
)
# Now we set up the model. For the soil model, we pick
# a model type and model args:
params = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = Soil.κ_dry(ρp, soil_ν, κ_solid, κ_air),
    κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, soil_ν, κ_ice),
    κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, soil_ν, κ_liq),
    ρc_ds = ρc_ds,
    ν = soil_ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = vanGenuchten(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r = θ_r,
    earth_param_set = earth_param_set,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
);
Δz = minimum(ClimaCore.Fields.Δz_field(ClimaLSM.coordinates(domain).subsurface))
sources = (Soil.PhaseChange{FT}(Δz),)
top_bc = ClimaLSM.Soil.AtmosDrivenFluxBC(atmos, radiation)
zero_flux = FluxBC((p, t) -> 0.0)
boundary_conditions =
    (; top = top_bc, bottom = (water = Soil.FreeDrainage(), heat = zero_flux))
soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = domain,
    boundary_conditions = boundary_conditions,
    sources = sources,
)

Y, p, cds = initialize(soil)

#Initial conditions
Y.soil.ϑ_l = SWC_1[1 + Int(round(t0 / DATA_DT))] # Get soil water content at t0
# Both data and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = TS_3[1 + Int(round(t0 / DATA_DT))] # Get soil temperature at t0
ρc_s = volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, Ref(soil.parameters))
Y.soil.ρe_int =
    volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, Ref(soil.parameters))
set_initial_aux_state! = make_set_initial_aux_state(soil)
set_initial_aux_state!(p, Y, t0);

# Simulation
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)

exp_tendency! = make_exp_tendency(soil)
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction((T_exp!) = exp_tendency!),
    Y,
    (t0, tf),
    p,
);
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
)

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climalsm_dir, "experiments/standalone/Soil/nsa_barrow/")


# Soil water content

plt1 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
    label = "Model, 2.5cm",
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0.05, 0.55],
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    color = "blue",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "20cm",
)

Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 4] for k in 1:1:length(sol.t)],
    label = "56m",
)

Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC_1, label = "Data, ?cm")
Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC_2, label = "Data, ?cm")

plt2 = Plots.plot(
    arm_seconds ./ 3600 ./ 24,
    P .* (-1e3 * 24 * 3600),
    label = "Data",
    ylabel = "Precipitation [mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    margin = 10Plots.mm,
    ylim = [-200, 0],
    size = (1500, 400),
)
Plots.plot(plt2, plt1, layout = grid(2, 1, heights = [0.2, 0.8]))
Plots.savefig(joinpath(savedir, "soil_water_content.png"))

# Temp
plt3 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt3,
    daily,
    [parent(sv.saveval[k].soil.T)[end] for k in 1:1:length(sol.t)],
    label = "Model, 2.5cm",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "Days",
    ylabel = "T [m/m]",
    color = "blue",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt3,
    daily,
    [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:1:length(sol.t)],
    label = "20cm",
)

Plots.plot!(
    plt3,
    daily,
    [parent(sv.saveval[k].soil.T)[end - 6] for k in 1:1:length(sol.t)],
    label = "1m",
)

Plots.plot!(
    plt3,
    daily,
    [parent(sv.saveval[k].soil.T)[end - 11] for k in 1:1:length(sol.t)],
    label = "3.3m",
)

Plots.plot!(plt3, seconds ./ 3600 ./ 24, TS_1, label = "Data, 1")
Plots.plot!(plt3, seconds ./ 3600 ./ 24, TS_2, label = "Data, 2")
Plots.plot!(plt3, seconds ./ 3600 ./ 24, TS_3, label = "Data, 3")


# LHF
lhf =
    [parent(sv.saveval[k].soil.sfc_conditions.lhf)[1] for k in 1:length(sol.t)]
Plots.plot(daily, lhf, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, LE, label = "Data")

# SHF
shf =
    [parent(sv.saveval[k].soil.sfc_conditions.shf)[1] for k in 1:length(sol.t)]
Plots.plot(daily, shf, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, H, label = "Data")

soil_Rn = [parent(sv.saveval[k].soil.R_n)[1] for k in 1:length(sol.t)]
soil_Rn_data = @. (SW_IN - SW_OUT + LW_IN - LW_OUT)
Plots.plot(daily, -soil_Rn, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, soil_Rn_data, label = "Data")
