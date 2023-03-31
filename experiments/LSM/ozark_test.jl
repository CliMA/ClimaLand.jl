using DiffEqCallbacks
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

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
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include("./ozark_drivers_from_data.jl")

# Now we set up the model. For the soil model, we pick
# a model type and model args:
nelements = 10
zmin = FT(-2)
zmax = FT(0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelements)
soil_ν = FT(0.55) # m3/m3, guess 
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(1.5) # unitless
soil_vg_α = FT(0.10682) # inverse meters. From Natan (10.9/MPa)
soil_vg_m = FT(1) - FT(1) / soil_vg_n # unitless
θ_r = FT(0) # m3/m3, guess 
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
radiative_transfer_args = (; parameters = BeerLambertParameters{FT}())
# Set up conductance
conductance_args = (; parameters = MedlynConductanceParameters{FT}())
# Set up photosynthesis
photosynthesis_args = (; parameters = FarquharParameters{FT}(Canopy.C3()))
# Set up plant hydraulics
RAI = FT(2) # m2/m2
SAI = FT(1) # m2/m2
LAI = FT(2) # m2/m2
h_stem = FT(10) # m, average full grown white oaks are ~20 m in height
h_leaf = FT(10) # m
area_index = (root = RAI, stem = SAI, leaf = LAI)
# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
n_stem = Int64(1)
n_leaf = Int64(1)
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [FT(0.0), h_stem, h_stem + h_leaf]
K_sat_plant = 1.8e-8 # m/s, guess, close to Natan
K_sat_root = FT(K_sat_plant) # m/s
K_sat_stem = FT(K_sat_plant)
K_sat_leaf = FT(K_sat_plant)
K_sat = (root = K_sat_root, stem = K_sat_stem, leaf = K_sat_leaf)
plant_vg_α = FT(0.002) # 1/m, matches Natan (fitted vG to Weibull)
plant_vg_n = FT(4.2) # unitless, matches Natan (fitted vG to Weibull)
plant_vg_m = FT(1) - FT(1) / plant_vg_n
plant_ν = FT(0.7) # guess, m3/m3
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(0.5)

# 0.5 is from Natan
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
z0_m = h_stem + h_leaf + FT(10)
z0_b = h_stem + h_leaf + FT(10)
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    LAI,
    h_stem + h_leaf,
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args =
    (; parameters = shared_params, domain = ClimaLSM.Point(; z_sfc = FT(0)))

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
ode! = make_ode_function(land)

#Initial conditions
Y.soil.ϑ_l = FT(0.3)
p_stem_0 = FT(-0.001e6 / 9800)
p_leaf_0 = FT(-0.002e6 / 9800)

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
t0 = FT(0);
N_days = 365
tf = FT(3600 * 24 * N_days)
dt = FT(100);

sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
daily = Array(2:(3600 * 24):(N_days * 3600 * 24))
hourly = Array(2:(3600):(N_days * 3600 * 24))
prob = ODEProblem(ode!, Y, (t0, tf), p);
cb =
    SavingCallback((u, t, integrator) -> copy(integrator.p), sv; saveat = daily)
sol =
    solve(prob, RK4(); dt = dt, saveat = daily, callback = cb, adaptive = false)

# Diagnostics
root_depths = parent(cds.subsurface.z)[:]
ϑ_stem =
    [parent(sol.u[k].canopy.hydraulics.ϑ_l[1])[1] for k in 1:1:length(sol.t)] # m3 m-3
ϑ_leaf =
    [parent(sol.u[k].canopy.hydraulics.ϑ_l[2])[1] for k in 1:1:length(sol.t)] # m3 m-3
S_l_stem = PlantHydraulics.effective_saturation.(plant_ν, ϑ_stem) # m3 m-3
S_l_leaf = PlantHydraulics.effective_saturation.(plant_ν, ϑ_leaf) # m3 m-3
p_stem =
    PlantHydraulics.water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        S_l_stem,
        plant_ν,
        plant_S_s,
    ) .* 0.0098 # m to MPa
p_leaf =
    PlantHydraulics.water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        S_l_leaf,
        plant_ν,
        plant_S_s,
    ) .* 0.0098 # m to MPa

ψ_soil_10_cm = [parent(sv.saveval[k].soil.ψ)[end] for k in 1:1:length(sol.t)] # m3 m-3
ψ_soil_20_cm =
    [parent(sv.saveval[k].soil.ψ)[end - 1] for k in 1:1:length(sol.t)] # m3 m-3
ψ_soil_30_cm =
    [parent(sv.saveval[k].soil.ψ)[end - 2] for k in 1:1:length(sol.t)] # m3 m-3

# Leaf and stem absolute pressure as function of time, Natan modeled leaf pressure, 
# and leaf pressure data at Ozark for all trees, and for white oak, for 1 year (2005)
plot1 = plot(
    sol.t ./ 3600 ./ 24,
    p_stem,
    label = "CliMA stem",
    color = :red,
    dpi = 400,
)

plot!(
    seconds ./ 3600 ./ 24,
    ylim = [-4, 2],
    our_year_lwp,
    label = "Natan leaf",
    color = :blue,
    xlabel = "t (days)",
    ylabel = "ψ plant (MPa)",
    legend = :topright,
    dpi = 400,
)

plot!(
    sol.t ./ 3600 ./ 24,
    p_leaf,
    label = "CliMA leaf",
    color = :green,
    dpi = 400,
)

plot!(
    unique_tree_observation_date_white_oak,
    mean_white_oak,
    color = :orange,
    label = "obs. mean white oak",
    dpi = 400,
)

plot!(
    day_of_year_of_observed_predawn_lwp_2005,
    observed_predawn_lwp_2005,
    seriestype = :scatter,
    markersize = 1,
    color = :black,
    label = "obs. Ozark, all species",
    legend = :topright,
    legendfontsize = 3,
    dpi = 400,
)

plot!(
    tree_observation_date_white_oak,
    white_oak,
    seriestype = :scatter,
    markerstrokecolor = :orange,
    color = :orange,
    markersize = 1,
    label = "obs. individual white oaks",
    legend = :topright,
    dpi = 400,
)

# Leaf absolute pressure as function of time for 6 Ozark tree species
plot2 = plot(
    tree_observation_date_white_oak,
    white_oak,
    xlabel = "t (days)",
    ylabel = "ψ leaf (MPa)",
    ylim = [-4, 2],
    label = "obs. white oak",
    legend = :topright,
    legendfontsize = 3,
    dpi = 400,
)

plot!(
    tree_observation_date_black_oak,
    black_oak,
    label = "obs. black oak",
    dpi = 400,
)

plot!(
    tree_observation_date_eastern_redcedar,
    eastern_redcedar,
    label = "obs. eastern redcedar",
    dpi = 400,
)

plot!(
    tree_observation_date_shagbark_hickory,
    shagbark_hickory,
    markersize = 1,
    label = "obs. shagbark hickory",
    dpi = 400,
)

plot!(
    tree_observation_date_sugar_maple,
    sugar_maple,
    label = "obs. sugar maple",
    dpi = 400,
)

plot!(
    tree_observation_date_white_ash,
    white_ash,
    label = "obs. white ash",
    dpi = 400,
)

# Augmented liquid fraction in stem and leaf as function of time
plot3 = plot(
    sol.t ./ 3600 ./ 24,
    ϑ_stem,
    xlabel = "t (days)",
    ylabel = "ϑ (m3 m-3)",
    color = :red,
    label = "ϑ stem",
    dpi = 400,
    legend = :bottomright,
    legendfontsize = 3,
)
plot!(sol.t ./ 3600 ./ 24, ϑ_leaf, color = :green, label = "ϑ leaf")
plot!(
    sol.t ./ 3600 ./ 24,
    plant_ν .* ones(length(sol.t), 1),
    color = :black,
    label = "plant porosity",
)

# SWC observation
plot4 = plot(
    seconds ./ 3600 ./ 24,
    our_year_swc_obs_surface,
    color = :black,
    xlabel = "t (days)",
    ylabel = "ϑ soil",
    label = "obs, sfc",
    legend = :bottomleft,
    legendfontsize = 3,
    dpi = 400,
)

plot!(
    seconds ./ 3600 ./ 24,
    our_year_swc_surface,
    color = :blue,
    label = "Natan, sfc",
)

plot!(
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
    label = "10cm",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
)

plot!(
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "20cm",
)

plot!(
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "30cm",
)

plot!(
    sol.t ./ 3600 ./ 24,
    [parent(sol.u[k].soil.ϑ_l)[end - 3] for k in 1:1:length(sol.t)],
    label = "40cm",
)

# Precipitation
plot5 = plot(
    sol.t / 3600 ./ 24,
    precipitation_function.(sol.t) * 3600 * 24 * 1e3,
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    label = "Fluxnet precipitation",
    dpi = 400,
    legendfontsize = 4,
    xlabel = "t (days)",
    ylabel = "Precip (mm/day)",
    legend = :bottomright,
)

plot!(xlabel = "t (days)", ylabel = "Precipitation (mm/day)")
plot(plot1, plot2, plot3, plot4, plot5, layout = 6)
savefig("./experiments/LSM/integrated_plant_soil_diagnostics_1.png")

# Transpiration
plot6 = plot(
    sol.t / 3600 ./ 24,
    transpiration_function.(sol.t) * 3600 * 24 * 1e3, # in (mm/day)
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    label = "Fluxnet ET",
    dpi = 400,
    xlabel = "t (days)",
    ylabel = "ET (mm/day)",
    legendfontsize = 4,
    legend = :topleft,
)

# Roots to stem flux
sum_flux_out_roots = [
    sum(
        flux.(
            root_depths,
            h_stem / 2,
            parent(sv.saveval[k].soil.ψ),
            p_stem[k],
            plant_vg_α,
            plant_vg_n,
            plant_vg_m,
            plant_ν,
            plant_S_s,
            K_sat[:root],
            K_sat[:stem],
        ) .* root_distribution.(root_depths) .* (
            vcat(root_depths, [0.0])[2:end] -
            vcat(root_depths, [0.0])[1:(end - 1)]
        ),
    ) for k in 1:1:length(sol.t)
] # in (m/s)

plot7 = plot(
    sol.t / 3600 ./ 24,
    sum_flux_out_roots * 3600 * 24 * 1e3, # in (mm/day)
    xlabel = "t (days)",
    ylabel = "Flux (mm/day)",
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    label = "roots to stem",
    dpi = 400,
    legendfontsize = 4,
    legend = :bottomright,
)

flux_out_stem = [
    flux(
        h_stem / 2,
        h_stem + h_leaf / 2,
        p_stem[k],
        p_leaf[k],
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        plant_ν,
        plant_S_s,
        K_sat[:stem],
        K_sat[:leaf],
    ) for k in 1:1:length(sol.t)
]

# Flux out stem
plot8 = plot(
    sol.t / 3600 ./ 24,
    flux_out_stem * 3600 * 24 * 1e3, # in (mm/day)
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    xlabel = "t (days)",
    ylabel = "Flux (mm/day)",
    label = "stem to leaf",
    legend = :topright,
    dpi = 400,
    legendfontsize = 4,
)

# Transpiration * LAI
plot9 = plot(
    sol.t / 3600 ./ 24,
    transpiration_function.(sol.t) * 3600 * 24 * 1e3 * LAI, # in (mm/day)
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    label = "LAI * Fluxnet ET",
    dpi = 400,
    legendfontsize = 4,
    legend = :topleft,
)
plot!(xlabel = "t (days)", ylabel = "LAI * ET (mm/day)")

# Flux out roots * RAI
plot10 = plot(
    sol.t / 3600 ./ 24,
    sum_flux_out_roots .* RAI * 3600 * 24 * 1e3, # in (mm/day)
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    xlabel = "t (days)",
    ylabel = "RAI * Flux (mm/day)",
    label = "roots to stem",
    legend = :topright,
    dpi = 400,
    legendfontsize = 4,
)

# Flux out stem * SAI
plot11 = plot(
    sol.t / 3600 ./ 24,
    flux_out_stem .* SAI * 3600 * 24 * 1e3, # in (mm/day)
    xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    xlabel = "t (days)",
    ylabel = "SAI * Flux (mm/day)",
    label = "stem to leaf",
    legend = :topright,
    dpi = 400,
    legendfontsize = 4,
)

plot(plot6, plot7, plot8, plot9, plot10, plot11, layout = 6)
savefig("./experiments/LSM/integrated_plant_soil_diagnostics_2.png")

# K(p)
function weibull(K_max::FT, b::FT, c::FT, p::FT) where {FT}
    K = K_max * exp(-((-p / b)^c))
    return K
end

function logistic(K_max::FT, a::FT, b::FT, p::FT) where {FT}
    K = K_max * (a + 1) / a * (1 - 1 / (1 + a * exp(b * p))) # Conductivity
    return K
end

a_stem = FT(0.1)
b_stem = FT(0.17 * 0.0098) # MPa-1 to m-1
b_natan = FT(-4 / 0.0098) # MPa to m. Natan's P63 for xylem conductance, see Holtzman clima parameter table march 2022, "The maximum xylem conductance is scaled by a PLC factor equal to exp(-(psi / xylem P63)^(Weibull exponent))"
c_natan = FT(4) # Natan's Weibull exponent for steepness of PLC and beta curves
mol_H2O_to_m3 = FT(18e-6)
K_sat_Natan = FT(10 * mol_H2O_to_m3 * 0.0098) # Conductance Natan (10 mol/s/m^2/MPa) to m^3/m^2/s/m
p_test = collect(range(-15, 2.0, 100)) ./ 0.0098
S_l_test =
    inverse_water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        p_test,
        plant_ν,
        plant_S_s,
    )
K_p_Natan = weibull.(K_sat_Natan * h_stem, b_natan, c_natan, p_test) # Conductivity
K_p_original_ozark_test =
    logistic.(K_sat_Natan * h_stem, a_stem, b_stem, p_test) # Conductivity
K_p_updated_ozark_test =
    hydraulic_conductivity.(K_sat_plant, plant_vg_m, S_l_test) # Conductivity

plot12 = plot(
    p_test * 0.0098, # MPa
    K_p_Natan * 3600 * 24 * 1e3,
    xlabel = "ψ (MPa)",
    ylabel = "K (mm/day)",
    dpi = 400,
    legendfontsize = 3,
    label = "K(ψ) Natan's Ozark test (Weibull)",
)

plot!(
    p_test * 0.0098, # MPa
    K_p_original_ozark_test * 3600 * 24 * 1e3,
    xlabel = "ψ (MPa)",
    ylabel = "K (mm/day)",
    dpi = 400,
    legendfontsize = 4,
    label = "K(ψ) original Ozark test (logistic)",
)

plot13 = plot(
    p_test * 0.0098, # MPa
    K_p_updated_ozark_test * 3600 * 24 * 1e3,
    xlabel = "ψ (MPa)",
    ylabel = "K (mm/day)",
    label = "K(ψ) updated Ozark test (vG)",
    dpi = 400,
    legendfontsize = 3,
)

p_test = collect(range(-15, 2, 100)) ./ 0.0098 # MPa to m
S_l_test =
    inverse_water_retention_curve.(
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        p_test,
        plant_ν,
        plant_S_s,
    )

plot14 = plot(
    S_l_test,
    p_test * 0.0098, # MPa
    ylabel = "ψ (MPa)",
    xlabel = "S_l (m3/m3)",
    dpi = 400,
    legendfontsize = 3,
    label = "ψ(S_l) with vG params used for vG K(P)",
)

plot(plot12, plot13, plot14, layout = 6)
savefig("./experiments/LSM/integrated_plant_soil_diagnostics_3.png")
