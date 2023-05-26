using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation

using ClimaLSM
using ClimaLSM.Domains: LSMSingleColumnDomain
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
        "experiments/LSM/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_parameters.jl"))
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_domain.jl"))
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_simulation.jl"))

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain.subsurface
soil_ps = Soil.RichardsParameters{FT}(
    soil_ν,
    soil_vg_α,
    soil_vg_n,
    soil_vg_m,
    soil_K_sat,
    soil_S_s,
    θ_r,
)

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.RichardsModel{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    radiative_transfer = Canopy.BeerLambertModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
)
# Individual Component arguments
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = BeerLambertParameters{FT}(;
        Ω = Ω,
        ld = ld,
        ρ_leaf = ρ_leaf,
        λ_γ = λ_γ,
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
area_index = (root = RAI, stem = SAI, leaf = LAI)
K_sat_root = FT(K_sat_plant) # m/s
K_sat_stem = FT(K_sat_plant)
K_sat_leaf = FT(K_sat_plant)
K_sat = (root = K_sat_root, stem = K_sat_stem, leaf = K_sat_leaf)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters{FT}(
    area_index,
    K_sat,
    plant_vg_α,
    plant_vg_n,
    plant_vg_m,
    plant_ν,
    plant_S_s,
    root_distribution,
)

plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

# Canopy component args
canopy_component_args = (;
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
)
# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    LAI,
    h_stem + h_leaf,
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = land_domain.surface)

# Integrated plant hydraulics and soil model
# The default is diagnostics transpiration. In this case, we are testing with prescribed
# but will change that soon.
land_input =
    (transpiration = transpiration, atmos = atmos, radiation = radiation)
land = SoilPlantHydrologyModel{FT}(;
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
Y.soil.ϑ_l = FT(0.4)
ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        [ψ_stem_0, ψ_leaf_0],
        plant_ν,
        plant_S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l[i] .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

update_aux! = make_update_aux(land)
update_aux!(p, Y, 0.0)

# Simulation
sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
prob = ODEProblem(exp_tendency!, Y, (t0, tf), p);
cb = SavingCallback(
    (u, t, integrator) -> copy(integrator.p),
    sv;
    saveat = hourly,
)
sol = solve(
    prob,
    timestepper;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = hourly,
)

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climalsm_dir, "experiments/LSM/ozark")
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for k in 1:length(sol.t)
]

plt1 = Plots.plot(size = (500, 700))
Plots.plot!(
    plt1,
    daily,
    model_GPP .* 1e6,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "GPP [μmol/mol]",
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    GPP .* 1e6,
    label = "Data",
    margins = 10Plots.mm,
    lalpha = 0.3,
)
plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    FT.(SW_IN),
    label = "Data",
    xlim = [minimum(daily), maximum(daily)],
    ylabel = " SW radiation (W/m^2)",
    size = (500, 700),
)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "GPP.png"))


T =
    [
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
measured_T = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)
plt1 = Plots.plot(size = (500, 700))
Plots.plot!(
    plt1,
    daily,
    T,
    label = "Modeled T",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "Vapor Flux [mm/day]",
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    measured_T,
    label = "Data (vapor flux)",
    margins = 10Plots.mm,
    lalpha = 0.3,
)

plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    VPD,
    label = "Measured",
    ylabel = "VPD (Pa)",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 700),
    margins = 10Plots.mm,
)
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "ET.png"))
plt1 = Plots.plot(size = (500, 700))

β = [parent(sv.saveval[k].canopy.hydraulics.β)[1] for k in 1:length(sol.t)]
plt1 = Plots.plot(
    daily,
    β,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "Moisture Stress",
    size = (500, 700),
)
g_stomata =
    [parent(sv.saveval[k].canopy.conductance.gs)[1] for k in 1:length(sol.t)]
plt2 = Plots.plot(
    daily,
    g_stomata,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "Stomatal conductance",
    ylim = [0, 0.3],
    size = (500, 700),
)
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "stomatal_conductance.png"))


plt1 = Plots.plot(size = (500, 700))
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
    label = "10cm",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0.05, 0.5],
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    margins = 10Plots.mm,
    color = "blue",
)

plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "20cm",
    color = "red",
)

plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "30cm",
    color = "purple",
)

Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC, label = "Data")
plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    P .* (-1e3 * 24 * 3600),
    label = "Data",
    ylabel = "Precipitation [mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 700),
)
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "soil_water_content.png"))

root_stem_flux = [
    sum(sv.saveval[k].root_extraction) .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]

stem_leaf_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[1] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
leaf_air_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[2] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]


plt1 = Plots.plot(
    daily,
    root_stem_flux,
    label = "Soil-root-stem flux",
    ylabel = "Within plant fluxes[mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 700),
    margins = 10Plots.mm,
)
Plots.plot!(plt1, daily, stem_leaf_flux, label = "Stem-leaf flux")
Plots.plot!(plt1, daily, leaf_air_flux, label = "Leaf-air flux")

lwp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[2] * 9800 for k in 1:length(sol.t)
]
swp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[1] * 9800 for k in 1:length(sol.t)
]

plt2 = Plots.plot(
    daily,
    lwp,
    label = "Leaf",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "days",
    ylabel = "Water potential (Pa)",
    size = (500, 700),
    left_margin = 10Plots.mm,
)

Plots.plot!(plt2, daily, swp, label = "Stem")

Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "plant_hydraulics.png"))


dt_model = sol.t[2] - sol.t[1]
dt_data = seconds[2] - seconds[1]
Plots.plot(daily, cumsum(T) * dt_model, label = "Model T")
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(measured_T[:]) * dt_data,
    label = "Data ET",
)
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(P[:]) * dt_data * (1e3 * 24 * 3600),
    label = "Data P",
)
Plots.plot!(ylabel = "∫ Water fluxes dt", xlabel = "Days", margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "cumul_p_et.png"))


rm(joinpath(savedir, "Artifacts.toml"))
