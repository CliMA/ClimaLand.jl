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
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
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
        "experiments/integrated/vaira_ranch/vaira_met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/vaira_ranch/vaira_domain.jl",
    ),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/vaira_ranch/vaira_parameters.jl",
    ),
)
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/vaira_ranch/vaira_simulation.jl",
    ),
)
# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_ps = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
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

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
    energy = Canopy.BigLeafEnergyModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args = (;
    parameters = AutotrophicRespirationParameters{FT}(;
        ne = ne,
        ηsl = ηsl,
        σl = σl,
        μr = μr,
        μs = μs,
        f1 = f1,
        f2 = f2,
    )
)
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters{FT}(;
        Ω = Ω,
        ld = ld,
        α_PAR_leaf = α_PAR_leaf,
        λ_γ_PAR = λ_γ_PAR,
        λ_γ_NIR = λ_γ_NIR,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
        n_layers = n_layers,
        ϵ_canopy = ϵ_canopy,
    )
)
# Set up conductance
conductance_args = (;
    parameters = MedlynConductanceParameters{FT}(;
        g1 = g1,
        Drel = Drel,
        g0 = g0,
    )
)
# Set up photosynthesis
photosynthesis_args = (;
    parameters = FarquharParameters{FT}(
        Canopy.C3();
        oi = oi,
        ϕ = ϕ,
        θj = θj,
        f = f,
        sc = sc,
        pc = pc,
        Vcmax25 = Vcmax25,
        Γstar25 = Γstar25,
        Kc25 = Kc25,
        Ko25 = Ko25,
        To = To,
        ΔHkc = ΔHkc,
        ΔHko = ΔHko,
        ΔHVcmax = ΔHVcmax,
        ΔHΓstar = ΔHΓstar,
        ΔHJmax = ΔHJmax,
        ΔHRd = ΔHRd,
    )
)
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)
# Canopy component args
canopy_component_args = (;
    autotrophic_respiration = autotrophic_respiration_args,
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
    energy = energy_args,
)
# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = canopy_domain)

# Integrated plant hydraulics and soil model
land_input = (atmos = atmos, radiation = radiation)
land = SoilCanopyModel{FT}(;
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)

#Initial conditions
Y.soil.ϑ_l = SWC_3[1 + Int(round(t0 / DATA_DT))] # Get soil water content at t0
# Both data and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = TS_3[1 + Int(round(t0 / DATA_DT))] # Get soil temperature at t0
ρc_s =
    volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, Ref(land.soil.parameters))
Y.soil.ρe_int =
    volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T_0,
        Ref(land.soil.parameters),
    )
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve(retention_model, ψ_leaf_0, plant_ν, plant_S_s)

Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction(plant_ν, S_l_ini)


Y.canopy.energy.T = TA[1 + Int(round(t0 / 1800))] # Get atmos temperature at t0

set_initial_aux_state! = make_set_initial_aux_state(land)
set_initial_aux_state!(p, Y, t0);

# Simulation
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)


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
savedir = joinpath(climalsm_dir, "experiments/integrated/vaira_ranch/")
# Number of datapoints per day
data_daily_points = Int64(86400 / DATA_DT)
# Number of model points per day
model_daily_points = Int64(86400 / n / dt)
# Scales data indices 0 to 24 in a day
data_daily_indices = range(0, step = DATA_DT / 3600, length = data_daily_points)
model_daily_indices =
    range(0, step = dt * n / 3600, length = model_daily_points)
# Number of data points per each model point (Ratio of data dt to model dt)
data_per_model = Int64(dt * n ÷ DATA_DT)

# This function will be used to compute averages over diurnal cycles. Input a
# data series to average over and the total number of days in the data series, 
# and it will return a vector of the average value at each timestamp in the 
# day averaged over every day in the series.
function diurnal_avg(series)
    num_days = Int64(N_days - N_days_spinup)
    daily_points = Int64(length(series) ÷ num_days) # Num datapoints per day
    daily_data = [
        series[i:1:(i + daily_points - 1)] for
        i in 1:daily_points:(daily_points * num_days)
    ]
    daily_avgs =
        [mean([daily_data[i][j] for i in 1:num_days]) for j in 1:daily_points]
    return daily_avgs
end
# Autotrophic Respiration
AR = [
    parent(sv.saveval[k].canopy.autotrophic_respiration.Ra)[1] for
    k in 1:length(sv.saveval)
]
AR_daily_avg_model = diurnal_avg(AR)

pltAR = Plots.plot(size = (1500, 400))

Plots.plot!(
    pltAR,
    model_daily_indices,
    AR_daily_avg_model .* 1e6,
    label = "Model",
    title = "AR [μmol/m^2/s]",
    xlabel = "Hour of day",
    ylabel = "Average over simulation",
    margin = 10Plots.mm,
)

Plots.savefig(joinpath(savedir, "AutoResp.png"))


# GPP
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in 1:length(sv.saveval)
]

GPP_daily_avg_model = diurnal_avg(model_GPP)
GPP_daily_avg_data =
    diurnal_avg(GPP[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])
RMSD =
    StatsBase.rmsd(
        GPP_daily_avg_model,
        GPP_daily_avg_data[1:data_per_model:end],
    ) * 1e6
R² =
    Statistics.cor(
        GPP_daily_avg_model,
        GPP_daily_avg_data[1:data_per_model:end],
    )^2

plt1 = Plots.plot(size = (800, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    GPP_daily_avg_data .* 1e6,
    label = "Data",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt1,
    model_daily_indices,
    GPP_daily_avg_model .* 1e6,
    xlabel = "Hour of day",
    ylabel = "GPP (mol/m²/s)",
    label = "Model",
    title = "GPP [μmol/m^2/s]: RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
)
Plots.savefig(joinpath(savedir, "GPP.png"))

# SW_OUT
SW_out_model = [parent(sv.saveval[k].SW_out)[1] for k in 1:length(sv.saveval)]
SW_out_model_avg = diurnal_avg(SW_out_model)
SW_out_data_avg =
    diurnal_avg(FT.(SW_OUT)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

RMSD = StatsBase.rmsd(SW_out_model_avg, SW_out_data_avg[1:data_per_model:end])
R² = Statistics.cor(SW_out_model_avg, SW_out_data_avg[1:data_per_model:end])^2

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    SW_out_data_avg,
    label = "Data",
    title = "Outgoing SW (W/m^2): RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
    xlabel = "Hour of day",
    margin = 10Plots.mm,
)
Plots.plot!(plt1, model_daily_indices, SW_out_model_avg, label = "Model")
Plots.savefig(joinpath(savedir, "SW_OUT.png"))

# LW_OUT
LW_out_model = [parent(sv.saveval[k].LW_out)[1] for k in 1:length(sv.saveval)]
LW_out_model_avg = diurnal_avg(LW_out_model)
LW_out_data_avg =
    diurnal_avg(FT.(LW_OUT)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])


RMSD = StatsBase.rmsd(LW_out_model_avg, LW_out_data_avg[1:data_per_model:end])
R² = Statistics.cor(LW_out_model_avg, LW_out_data_avg[1:data_per_model:end])^2

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    LW_out_data_avg,
    label = "Data",
    title = "Outgoing LW (W/m^2): RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
    xlabel = "Hour of day",
    margin = 10Plots.mm,
)
Plots.plot!(plt1, model_daily_indices, LW_out_model_avg, label = "Model")
Plots.savefig(joinpath(savedir, "LW_out.png"))


T =
    [
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
E =
    [parent(sv.saveval[k].soil_evap)[1] for k in 1:length(sol.t)] .* (1e3 * 24 * 3600)
measured_T = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)

ET_avg_model = diurnal_avg(T .+ E)
ET_avg_data =
    diurnal_avg(measured_T[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

RMSD = StatsBase.rmsd(ET_avg_model, ET_avg_data[1:data_per_model:end])
R² = Statistics.cor(ET_avg_model, ET_avg_data[1:data_per_model:end])^2

plt1 = Plots.plot(size = (800, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    ET_avg_data,
    label = "Data ET",
    xlabel = "Hour of day",
    ylabel = "ET",
    margins = 10Plots.mm,
    title = "Average Daily ET: RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
)
Plots.plot!(plt1, model_daily_indices, ET_avg_model, label = "Model ET")
Plots.savefig(joinpath(savedir, "ET.png"))

# Water stress factor
β = [parent(sv.saveval[k].canopy.hydraulics.β)[1] for k in 1:length(sol.t)]
plt1 =
    Plots.plot(size = (1500, 400), xlabel = "Day of year", margin = 10Plots.mm)
Plots.plot!(
    plt1,
    daily,
    β,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Moisture stress factor",
)
Plots.savefig(joinpath(savedir, "moisture_stress.png"))

# Stomatal conductance
g_stomata =
    [parent(sv.saveval[k].canopy.conductance.gs)[1] for k in 1:length(sol.t)]
plt1 =
    Plots.plot(size = (1500, 400), xlabel = "Day of year", margin = 10Plots.mm)
Plots.plot!(
    plt1,
    daily,
    g_stomata,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Stomatal conductance (mol/m^2/s)",
)
Plots.savefig(joinpath(savedir, "stomatal_conductance.png"))

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

Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC_1, label = "Data, 5cm")
Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC_2, label = "Data, 20cm")
Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC_3, label = "Data, 50cm")
plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
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

# Sensible Heat Flux

SHF_soil = [parent(sv.saveval[k].soil_shf)[1] for k in 1:length(sol.t)]
SHF_canopy =
    [parent(sv.saveval[k].canopy.energy.shf)[1] for k in 1:length(sol.t)]
SHF = SHF_soil + SHF_canopy
SHF_soil_avg_model = diurnal_avg(SHF_soil)
SHF_canopy_avg_model = diurnal_avg(SHF_canopy)
SHF_avg_model = diurnal_avg(SHF)
SHF_avg_data = diurnal_avg(H[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

RMSD = StatsBase.rmsd(SHF_avg_model, SHF_avg_data[1:data_per_model:end])
R² = Statistics.cor(SHF_avg_model, SHF_avg_data[1:data_per_model:end])^2


plt1 = Plots.plot(size = (800, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    SHF_avg_data,
    xlabel = "Hour of day",
    ylabel = "SHF (W/m²)",
    label = "Data",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    model_daily_indices,
    SHF_avg_model,
    label = "Model_total",
    title = "Sensible Heat Flux: RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
)
Plots.plot!(
    plt1,
    model_daily_indices,
    SHF_soil_avg_model,
    label = "SHF_soil",
    color = "blue",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt1,
    model_daily_indices,
    SHF_canopy_avg_model,
    label = "SHF_canopy",
    color = "cyan",
    margin = 10Plots.mm,
)
Plots.savefig(joinpath(savedir, "shf.png"))

# Latent Heat Flux

LHF_soil = [parent(sv.saveval[k].soil_lhf)[1] for k in 1:length(sol.t)]
LHF_canopy =
    [parent(sv.saveval[k].canopy.energy.lhf)[1] for k in 1:length(sol.t)]
LHF = LHF_soil + LHF_canopy
LHF_soil_avg_model = diurnal_avg(LHF_soil)
LHF_canopy_avg_model = diurnal_avg(LHF_canopy)
LHF_avg_model = diurnal_avg(LHF)
LHF_avg_data = diurnal_avg(LE[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

RMSD = StatsBase.rmsd(LHF_avg_model, LHF_avg_data[1:data_per_model:end])
R² = Statistics.cor(LHF_avg_model, LHF_avg_data[1:data_per_model:end])^2

plt1 = Plots.plot(size = (800, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    LHF_avg_data,
    xlabel = "Hour of day",
    ylabel = "LHF (W/m²)",
    label = "Data",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    model_daily_indices,
    LHF_avg_model,
    label = "Model_total",
    title = "Latent Heat Flux: RMSD = $(round(RMSD, digits = 2)), R² = $(round(R², digits = 2))",
)
Plots.plot!(
    plt1,
    model_daily_indices,
    LHF_soil_avg_model,
    label = "LHF_soil",
    color = "blue",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt1,
    model_daily_indices,
    LHF_canopy_avg_model,
    label = "LHF_canopy",
    color = "cyan",
    margin = 10Plots.mm,
)
Plots.savefig(joinpath(savedir, "lhf.png"))


# Cumulative ET
dt_model = sol.t[2] - sol.t[1]
dt_data = seconds[2] - seconds[1]
# Find which index in the data our simulation starts at:
idx = argmin(abs.(seconds .- sol.t[1]))
Plots.plot(
    seconds ./ 24 ./ 3600,
    cumsum(measured_T[:]) * dt_data,
    label = "Data ET",
)
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(P[:]) * dt_data * (1e3 * 24 * 3600),
    label = "Data P",
)
Plots.plot!(
    daily,
    cumsum(T .+ E) * dt_model .+ cumsum(measured_T[:])[idx] * dt_data,
    label = "Model ET",
)

Plots.plot!(ylabel = "∫ Water fluxes dt", xlabel = "Days", margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "cumul_p_et.png"))

# Leaf Water Potentials
lwp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[1] * 9800 for k in 1:length(sol.t)
]
ψ_soil = [
    sum(
        parent(sv.saveval[k].soil.ψ) .*
        root_distribution.(parent(cds.subsurface.z)),
    ) / sum(root_distribution.(parent(cds.subsurface.z))) * 9800 for
    k in 1:length(sol.t)
]


plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    daily,
    lwp,
    label = "Model, Leaf",
    title = "Water potentials",
    xlabel = "Day of year",
    margin = 10Plots.mm,
    xlim = [minimum(daily), maximum(daily)],
)
Plots.plot!(plt1, daily, ψ_soil, label = "Model, Mean soil")

Plots.savefig(joinpath(savedir, "leaf_water_potential.png"))

# Temperatures
soil_T_sfc = [parent(sv.saveval[k].soil.T)[end] for k in 1:length(sol.t)]
soil_T_sfc_avg = diurnal_avg(soil_T_sfc)

canopy_T = [
    parent(
        ClimaLSM.Canopy.canopy_temperature(
            land.canopy.energy,
            land.canopy,
            sol.u[k],
            sv.saveval[k],
            sol.t[k],
        ),
    )[1] for k in 1:length(sol.t)
]
canopy_T_avg = diurnal_avg(canopy_T)

TA_avg = diurnal_avg(FT.(TA)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])
TS_avg = diurnal_avg(FT.(TS_1)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    TS_avg,
    label = "Soil-D",
    title = "Temperature",
)
Plots.plot!(plt1, data_daily_indices, TA_avg, label = "Atmos-D")

Plots.plot!(plt1, model_daily_indices, soil_T_sfc_avg, label = "Soil-M-2.5cm")

Plots.plot!(plt1, model_daily_indices, canopy_T_avg, label = "Canopy-M")
Plots.plot!(plt1, xlabel = "Hour of day", ylabel = "Average over Simulation")
Plots.plot!(plt1, margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "temperature.png"))

# Soil Temperature
soil_T_1 = [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:length(sol.t)]# 20cm
soil_T_2 = [parent(sv.saveval[k].soil.T)[end - 3] for k in 1:length(sol.t)]#36cm
soil_T_3 = [parent(sv.saveval[k].soil.T)[end - 5] for k in 1:length(sol.t)] # 80cm

TA_avg = diurnal_avg(FT.(TA)[Int64(t_spinup ÷ DATA_DT + 1):Int64(tf ÷ DATA_DT)])

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    diurnal_avg(FT.(TS_1)[Int64(t_spinup ÷ DATA_DT + 1):Int64(tf ÷ DATA_DT)]),
    label = "Data L1",
    title = "Temperature",
)
Plots.plot!(
    plt1,
    data_daily_indices,
    diurnal_avg(FT.(TS_2)[Int64(t_spinup ÷ DATA_DT + 1):Int64(tf ÷ DATA_DT)]),
    label = "Data L2",
    title = "Temperature",
)
Plots.plot!(
    plt1,
    data_daily_indices,
    diurnal_avg(FT.(TS_3)[Int64(t_spinup ÷ DATA_DT + 1):Int64(tf ÷ DATA_DT)]),
    label = "Data L3",
    title = "Temperature",
)
#Plots.plot!(plt1, data_daily_indices, TA_avg, label = "Air")
Plots.plot!(
    plt1,
    model_daily_indices,
    diurnal_avg(soil_T_1),
    label = "Model L1",
)
Plots.plot!(
    plt1,
    model_daily_indices,
    diurnal_avg(soil_T_2),
    label = "Model L2",
)
Plots.plot!(
    plt1,
    model_daily_indices,
    diurnal_avg(soil_T_3),
    label = "Model L3",
)
Plots.plot!(plt1, xlabel = "Hour of day", ylabel = "Average over Simulation")
Plots.plot!(plt1, margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "soil_temperature.png"))



# Ground heat flux
Δz = parent(cds.subsurface.z)[end] - parent(cds.subsurface.z)[end - 2]
first_layer_flux = [
    -parent(sv.saveval[k].soil.κ)[1] * (
        parent(sv.saveval[k].soil.T)[end] -
        parent(sv.saveval[k].soil.T)[end - 2]
    ) / Δz for k in 1:length(sol.t)
]
G_model = [
    (
        parent(sv.saveval[k].soil_shf)[1] + parent(sv.saveval[k].soil_lhf)[1] -
        parent(sv.saveval[k].soil_LW_n)[1] -
        parent(sv.saveval[k].soil_SW_n)[1]
    ) for k in 1:length(sol.t)
]
canopy_G = [
    (
        parent(sv.saveval[k].canopy.energy.shf)[1] +
        parent(sv.saveval[k].canopy.energy.lhf)[1] -
        parent(sv.saveval[k].canopy.radiative_transfer.LW_n)[1] -
        parent(sv.saveval[k].canopy.radiative_transfer.SW_n)[1]
    ) for k in 1:length(sol.t)
]

G_model_avg = diurnal_avg(G_model)
canopy_G_avg = diurnal_avg(canopy_G)

G_data_avg = diurnal_avg(FT.(G)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])
Rn = (SW_IN .- SW_OUT) .+ (LW_IN .- LW_OUT)
G_alternate_data_avg = diurnal_avg(
    FT.(H .+ LE .- Rn)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)],
)
HplusL_avg =
    diurnal_avg(FT.(H .+ LE)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])
RminusG_avg =
    diurnal_avg(FT.(Rn .- G)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)])

plt1 = Plots.plot(data_daily_indices, -1 .* G_data_avg, label = "Data: -G")
Plots.plot!(
    plt1,
    data_daily_indices,
    G_alternate_data_avg,
    label = "Data: (H+L-Rn)_site",
    xlabel = "Hour of day",
)
Plots.plot!(plt1, ylabel = "Flux (W/m^2)", title = "Energy balance at the site")

plt2 = Plots.scatter(
    HplusL_avg,
    RminusG_avg,
    label = "Diurnally averaged data",
    xlabel = "H+L",
    ylabel = "R-G",
)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "energy_balance_data.png"))

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    -1 .* G_data_avg,
    label = "Data: -G",
    margins = 10Plots.mm,
    xlabel = "Hour of day",
)
Plots.plot!(
    plt1,
    data_daily_indices,
    canopy_G_avg,
    label = "Model: (H+L-Rn)_canopy",
    margins = 10Plots.mm,
    xlabel = "Hour of day",
)
Plots.plot!(
    plt1,
    model_daily_indices,
    G_model_avg,
    label = "Model: (H+L-Rn)_soil",
    title = "Ground Heat Flux [W/m^2]",
)


Plots.plot!(
    plt1,
    model_daily_indices,
    G_alternate_data_avg,
    label = "Data: (H+L-Rn)_site",
)
Plots.plot!(
    plt1,
    data_daily_indices,
    diurnal_avg(first_layer_flux),
    label = "Model: -κ∂T∂z|_10cm",
)
Plots.savefig(joinpath(savedir, "ground_heat_flux.png"))


# Aerodynamic resistance
r_ae = [parent(sv.saveval[k].canopy.energy.r_ae)[1] for k in 1:length(sol.t)]
r_ae = min.(r_ae, 1e4)
r_ae_avg = diurnal_avg(r_ae)
R = FT(LSMP.gas_constant(earth_param_set))

r_canopy =
    1.0 ./ [
        ClimaLSM.Canopy.upscale_leaf_conductance(
            g_stomata[k],
            parent(sv.saveval[k].canopy.hydraulics.area_index.leaf)[1],
            atmos.T(sv.t[k]),
            R,
            atmos.P(sv.t[k]),
        ) for k in 1:length(sol.t)
    ]

plt1 = Plots.plot(size = (1500, 400))
Plots.plot!(
    plt1,
    data_daily_indices,
    r_ae_avg,
    label = "Aerodynamic",
    margins = 10Plots.mm,
    xlabel = "Hour of day",
    title = "Aerodynamic Resistance [s/m]",
)

Plots.plot!(
    plt1,
    data_daily_indices,
    diurnal_avg(r_canopy),
    label = "Canopy",
    margins = 10Plots.mm,
    xlabel = "Hour of day",
    title = "Resistance [s/m]",
    yaxis = :log,
)
Plots.savefig(joinpath(savedir, "r_ae.png"))


# Run script with comand line argument "save" to save model output to CSV
if length(ARGS) ≥ 1 && ARGS[1] == "save"
    # Formats fields as semicolon seperated strings
    field_to_array = (field) -> join(parent(field), ';')
    # Recursively unpacks a nested NamedTuple of fields into an array of strings
    function unpack(tup, data)
        for entry in tup
            if entry isa NamedTuple
                unpack(entry, data)
            else
                push!(data, field_to_array(entry))
            end
        end
    end
    # Recursively extracts the names of all fields in a nested namedTuple
    function extract_names(nt, names)
        for entry in pairs(nt)
            if entry[2] isa NamedTuple
                extract_names(entry[2], names)
            else
                push!(names, entry[1])
            end
        end
    end
    # Collect unpacked data from each timestep into an array
    timestamps = [[]]
    push!(timestamps[1], "Timestep")
    extract_names(sv.saveval[1], timestamps[1])
    local cnt = 0
    for timestamp in sv.saveval
        cnt = cnt + 1
        save_data = Any[cnt]
        unpack(timestamp, save_data)
        push!(timestamps, save_data)
    end
    # Write all data to a csv file
    writedlm(joinpath(savedir, "model_output.csv"), timestamps, ',')
    @info "Saved model output to $(savedir)model_output.csv"
end

rm(joinpath(savedir, "Artifacts.toml"))
