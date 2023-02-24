using DiffEqCallbacks
using ArtifactWrappers
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using ClimaCore
import CLIMAParameters as CP
using DelimitedFiles
using Dierckx
using Plots
using Statistics
using Dates


if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.PlantHydraulics
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

const FT = Float64

# Fluxnet driver data (precipitation, ET) and test data for comparison (soil moisture, leaf water potential)
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/3d03opwwczvxia0fks38uxr62qne57f0.csv",
    filename = "Ozark_test_param_set.csv",
)
dataset = ArtifactWrapper(@__DIR__, "driver", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "Ozark_test_param_set.csv")
driver_data = readdlm(data, ',', skipstart = 1)

# Natan's model output data for comparison
af2 = ArtifactFile(
    url = "https://caltech.box.com/shared/static/qkm7lrqaphlehe6p5hvccldhie6gxo0v.csv",
    filename = "holtzman_clima_output_april1.csv",
)
dataset2 = ArtifactWrapper(@__DIR__, "Natan", ArtifactFile[af2]);
dataset_path2 = get_data_folder(dataset2);
data2 = joinpath(dataset_path2, "holtzman_clima_output_april1.csv")
Natan_data = readdlm(data2, ',', skipstart = 1)

# 2005 data
precip_flux = driver_data[17569:2:35088, 5] # m/s
et_flux = driver_data[17569:2:35088, 7] # m/s
tree_species = driver_data[125:440, 15] # names of species
observed_predawn_lwp_2005 = driver_data[125:440, 16] # MPa
white_oak = observed_predawn_lwp_2005[tree_species .== "white oak"] # MPa
black_oak = observed_predawn_lwp_2005[tree_species .== "black oak"] # MPa
eastern_redcedar =
    observed_predawn_lwp_2005[tree_species .== "eastern redcedar"] # MPa
shagbark_hickory =
    observed_predawn_lwp_2005[tree_species .== "shagbark hickory"] # MPa
sugar_maple = observed_predawn_lwp_2005[tree_species .== "sugar maple"] # MPa
white_ash = observed_predawn_lwp_2005[tree_species .== "white ash"] # MPa
day_of_year_of_observed_predawn_lwp_2005 = driver_data[125:440, 14] # unitless, number of day in year
tree_observation_date_white_oak =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "white oak"] # MPa
tree_observation_date_black_oak =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "black oak"] # MPa
tree_observation_date_eastern_redcedar =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "eastern redcedar"] # MPa
tree_observation_date_shagbark_hickory =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "shagbark hickory"] # MPa
tree_observation_date_sugar_maple =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "sugar maple"] # MPa
tree_observation_date_white_ash =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "white ash"] # MPa
unique_tree_observation_date_white_oak = unique(tree_observation_date_white_oak) #  unitless, number of day in year
mean_white_oak = [
    mean(
        white_oak[tree_observation_date_white_oak .== unique_tree_observation_date_white_oak[k]],
    ) for k in 1:1:length(unique_tree_observation_date_white_oak)
]

t = Array(0:(60 * 60):((length(precip_flux) - 1) * 60 * 60)) # s
p_spline = Spline1D(t, -precip_flux) # m/s
et_spline = Spline1D(t, et_flux) # m/s

# Natan's data
natan_et = Natan_data[2:end, 15] .* 18 / 1e3 / 1e3 # m^3/m^2/s
swc_column = Natan_data[2:end, 20] # m^3/m^3
swc_surface = Natan_data[2:end, 19] # m^3/m^3
swc_obs_surface = Natan_data[2:end, 12] # m^3/m^3
lwp = Natan_data[2:end, 18] # MPa
dates = Natan_data[2:end, 1]
Natan_data = nothing
dates_julia = tryparse.(DateTime, dates)
our_year = dates_julia[Dates.year.(dates_julia) .== 2005]
seconds = Dates.value.(our_year .- our_year[1]) ./ 1000
our_year_swc_column = FT.(swc_column[Dates.year.(dates_julia) .== 2005])
our_year_swc_surface = FT.(swc_surface[Dates.year.(dates_julia) .== 2005])
our_year_swc_obs_surface =
    FT.(swc_obs_surface[Dates.year.(dates_julia) .== 2005])
our_year_lwp = FT.(lwp[Dates.year.(dates_julia) .== 2005])
our_year_et = FT.(natan_et[Dates.year.(dates_julia) .== 2005])
natan_et_spline = Spline1D(seconds, our_year_et)

# Pressure head takes in S_l; so assuming our_year_swc_obs_surface approx = ϑ_soil
soil_ν = FT(0.55) # guess, m3/m3
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3 * 0.0098) # 1/MPa to 1/m, Range is 0.001/MPa -> 1/MPa. Range from: Christoffersen 2016, Fig. S2.3, Scholz 2007; TRY, Xylem Functional Traits (XFT) Database, datasetID 241, traitID 1098. Calculated from from stem stem water capacitance/water density.
soil_vg_n = FT(1.5) # unitless
soil_vg_α = FT(0.10682) # 1/m, from Natan (10.9/MPa)
soil_vg_m = FT(1) - FT(1) / soil_vg_n
θ_r = FT(0) # m3/m3
ψ_observation =
    pressure_head.(
        soil_vg_α,
        soil_vg_n,
        soil_vg_m,
        θ_r,
        our_year_swc_obs_surface,
        soil_ν,
        soil_S_s,
    ) # m
sp_spline = Spline1D(seconds, ψ_observation) # m

precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
transpiration_function(t::FT) where {FT} = et_spline(t) # m/s
prescribed_soil_pressure_function(t::FT) where {FT} = sp_spline(t) # m/s

saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
earth_param_set = create_lsm_parameters(FT)

# Plant hydraulics params
# Two compartments of equal volume and height to begin
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

# 0.5 is from Natan
function root_distribution(z::T) where {T}
    return T(1.0 / 0.5) * exp(z / T(0.5)) # 1/m
end
earth_param_set = create_lsm_parameters(FT)
zmin = FT(-2.0)
zmax = FT(0.0)

plant_hydraulics_domain = ClimaLSM.Domains.Point(; z_sfc = zmax)

plant_hydraulics_ps =
    PlantHydraulics.PlantHydraulicsParameters{FT, typeof(earth_param_set)}(
        area_index,
        K_sat,
        plant_vg_α,
        plant_vg_n,
        plant_vg_m,
        plant_ν,
        plant_S_s,
        root_distribution,
        earth_param_set,
    )

transpiration =
    PrescribedTranspiration{FT}((t::FT) -> transpiration_function(t))
plant_hydraulics_args = (
    domain = plant_hydraulics_domain,
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

nelements = 10
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

# Integrated plant hydraulics and soil model
land_args = (
    precipitation = precipitation_function,
    transpiration = transpiration_function,
)

land = SoilPlantHydrologyModel{FT}(;
    land_args = land_args,
    soil_model_type = Soil.RichardsModel{FT},
    soil_args = soil_args,
    vegetation_model_type = PlantHydraulics.PlantHydraulicsModel{FT},
    vegetation_args = plant_hydraulics_args,
)
Y, p, cds = initialize(land)
ode! = make_ode_function(land)

# IC from equilibrium run
ic = [
    0.355686,
    0.354263,
    0.352855,
    0.351462,
    0.350083,
    0.348718,
    0.347368,
    0.346031,
    0.344708,
    0.343398,
]
root_depths =
    sort(unique(parent(ClimaLSM.Domains.coordinates(soil_args.domain).z)[:]))
ic_spline = Spline1D(root_depths, ic)
Y.soil.ϑ_l = ic_spline.(cds.subsurface.z)

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
    Y.vegetation.ϑ_l[i] .= augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

update_aux! = make_update_aux(land)
update_aux!(p, Y, 0.0)

# Sim
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
ϑ_stem = [parent(sol.u[k].vegetation.ϑ_l[1])[1] for k in 1:1:length(sol.t)] # m3 m-3
ϑ_leaf = [parent(sol.u[k].vegetation.ϑ_l[2])[1] for k in 1:1:length(sol.t)] # m3 m-3
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
savefig("./experiments/LSM/coupled_plant_soil_diagnostics_1.png")

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
savefig("./experiments/LSM/coupled_plant_soil_diagnostics_2.png")

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
savefig("./experiments/LSM/coupled_plant_soil_diagnostics_3.png")
