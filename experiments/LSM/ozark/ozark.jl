using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation

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
        "experiments/LSM/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_parameters.jl"))
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_domain.jl"))
include(joinpath(climalsm_dir, "experiments/LSM/ozark/ozark_simulation.jl"))

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelements)
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
        ψc = ψc,
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

canopy_model_args =
    (; parameters = shared_params, domain = ClimaLSM.Point(; z_sfc = zmax))

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
Y.soil.ϑ_l = FT(0.35)
p_stem_0 = FT(-1e6 / 9800)
p_leaf_0 = FT(-2e6 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        [p_stem_0, p_leaf_0],
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
    saveat = halfhourly,
)
sol = solve(
    prob,
    timestepper;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = halfhourly,
)

# Plotting
hours = sol.t ./ 3600
savedir = joinpath(climalsm_dir, "experiments/LSM/ozark")
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for k in 1:length(sol.t)
]
Plots.plot(
    hours,
    model_GPP .* 1e6,
    label = "Model",
    xlim = [minimum(hours), maximum(hours)],
    xlabel = "hours",
    ylabel = "GPP [μmol/mol]",
)
Plots.plot!(seconds ./ 3600, GPP .* 1e6, label = "Data")
Plots.savefig(joinpath(savedir, "GPP.png"))
T =
    [
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
measured_T = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)
Plots.plot(
    hours,
    T,
    label = "Model",
    xlim = [minimum(hours), maximum(hours)],
    xlabel = "hours",
    ylabel = "T [mm/day]",
)
Plots.plot!(seconds ./ 3600, measured_T, label = "Data")
Plots.savefig(joinpath(savedir, "T.png"))

plt1 = Plots.plot()
Plots.plot!(
    plt1,
    hours,
    [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
    label = "10cm",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [minimum(hours), maximum(hours)],
    ylim = [0.2, 0.5],
    xlabel = "Hours",
    ylabel = "SWC [m/m]",
)

plot!(
    plt1,
    hours,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "20cm",
)

plot!(
    plt1,
    hours,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "30cm",
)

plot!(
    plt1,
    hours,
    [parent(sol.u[k].soil.ϑ_l)[end - 3] for k in 1:1:length(sol.t)],
    label = "40cm",
)
Plots.plot!(plt1, seconds ./ 3600, SWC, label = "Data")
plt2 = Plots.plot(
    seconds ./ 3600,
    P .* (1e3 * 24 * 3600),
    label = "Data",
    ylabel = "Precipitation [mm/day]",
    xlabel = "Hours",
    xlim = [minimum(hours), maximum(hours)],
    ylim = [0, 200],
)
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "SWC.png"))
