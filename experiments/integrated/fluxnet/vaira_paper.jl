import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates
using ClimaDiagnostics
using ClimaUtilities

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
site_ID = "US-Var"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    g1,
    Drel,
    g0,
    Vcmax25,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    G_Function,
    n_stem,
    n_leaf,
    h_leaf,
    h_stem,
    h_canopy,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# Construct the ClimaLand domain to run the simulation on
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# Set up the timestepping information for the simulation
dt = Float64(450) # 7.5 minutes

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
(start_date, _) = FluxnetSimulations.get_data_dates(site_ID, time_offset)
stop_date = start_date + Year(1)
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
)


# Now we set up the model. For the soil model, we pick
# a model type and package up parameters.
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)

forcing = (; atmos, radiation)
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

soil = Soil.EnergyHydrology{FT}(
    soil_domain,
    forcing,
    toml_dict;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo = soil_albedo,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
    retention_parameters,
    composition_parameters,
    S_s = soil_S_s,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
)

# Soil microbes model
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers, toml_dict)

# Now we set up the canopy model, one component at a time.
# Set up radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = G_Function,
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    surface_domain,
    toml_dict;
    radiation_parameters,
    ϵ_canopy,
)

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(surface_domain, toml_dict; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(surface_domain, toml_dict; photosynthesis_parameters)

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = land_domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# Set up optimal LAI model (loads spatially varying GSL and A0_annual)
lai_model = Canopy.OptimalLAIModel{FT}(surface_domain, toml_dict)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
RAI = maxLAI * f_root_to_shoot
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    surface_domain,
    toml_dict;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
)

height = h_stem + h_leaf
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

# Set up the energy model
energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

# Combine the components into a CanopyModel
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    lai_model,
    hydraulics,
    energy,
    biomass,
)

# Snow model
snow = Snow.SnowModel(
    FT,
    surface_domain,
    forcing,
    toml_dict,
    dt;
    prognostic_land_components,
)

# Integrated plant hydraulics, soil, and snow model
land = LandModel{FT}(canopy, snow, soil, soilco2);
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
dt_save = Second(1800)
updateat = Second(dt)
saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, dt_save)
sv = saving_cb.affect!.saved_values
user_callbacks = (saving_cb,)
solver_kwargs = (; saveat = dt_save)
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    set_ic!,
    updateat,
    solver_kwargs,
    user_callbacks,
    diagnostics = (),
)
@time sol = solve!(simulation);

function compute_monthly_avg(x, dates)
    moy = month.(dates)
    monthly_avg = zeros(12)
    for i in 1:12
        moy_mask = moy .== i
        monthly_avg[i] = mean(x[moy_mask])
    end
    return monthly_avg
end

comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir =
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/$(site_ID)/out")
mkpath(savedir)

spinup = Day(20)
data_id_post_spinup =
    (comparison_data.UTC_datetime .> start_date .+ spinup) .&&
    (comparison_data.UTC_datetime .<= stop_date)
data_dates = comparison_data.UTC_datetime[data_id_post_spinup]
model_dates = Second.(float.(sol.t)) .+ start_date
model_id_post_spinup =
    (model_dates .> start_date .+ spinup) .&& (model_dates .<= stop_date)
model_dates = model_dates[model_id_post_spinup]
model_idx1 = findfirst(model_id_post_spinup)
model_idxend = findlast(model_id_post_spinup)
# Extract model output from the saved output
θ1 = [
    parent(sol.u[k].soil.ϑ_l .+ sol.u[k].soil.θ_i)[end] for
    k in model_idx1:1:model_idxend
]
θ3 = [
    parent(sol.u[k].soil.ϑ_l .+ sol.u[k].soil.θ_i)[end - 10] for
    k in model_idx1:1:model_idxend
]
T_soil1 = [parent(sv.saveval[k].soil.T)[end] for k in model_idx1:1:model_idxend]; # 2cm
T_soil5 =
    [parent(sv.saveval[k].soil.T)[end - 16] for k in model_idx1:1:model_idxend]
SW_u = [parent(sv.saveval[k].SW_u)[1] for k in model_idx1:1:model_idxend]
LHF = [
    parent(
        sv.saveval[k].snow.turbulent_fluxes.lhf .+
        sv.saveval[k].soil.turbulent_fluxes.lhf .+
        sv.saveval[k].canopy.turbulent_fluxes.lhf,
    )[1] for k in model_idx1:1:model_idxend
]
SHF = [
    parent(
        sv.saveval[k].canopy.turbulent_fluxes.shf .+
        sv.saveval[k].snow.turbulent_fluxes.shf .+
        sv.saveval[k].soil.turbulent_fluxes.shf,
    )[1] for k in model_idx1:1:model_idxend
]

GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in model_idx1:1:model_idxend
]


# Monthly averages
GPP = GPP .* 1e6 # mol to μmol
GPP_data = comparison_data.gpp[data_id_post_spinup] .* 1e6
GPP_model_monthly = compute_monthly_avg(GPP, model_dates)
GPP_data_monthly = compute_monthly_avg(GPP_data, data_dates)

# 
SW_u_model_monthly = compute_monthly_avg(SW_u, model_dates)
SW_u_data = comparison_data.swu[data_id_post_spinup]
SW_u_data_monthly = compute_monthly_avg(SW_u_data, model_dates)

LHF_data = comparison_data.lhf[data_id_post_spinup]
LHF_model_monthly = compute_monthly_avg(LHF, model_dates)
LHF_data_monthly = compute_monthly_avg(LHF_data, data_dates)

SHF_data = comparison_data.shf[data_id_post_spinup]
SHF_model_monthly = compute_monthly_avg(SHF, model_dates)
SHF_data_monthly = compute_monthly_avg(SHF_data, data_dates)
fig = Figure(size = (800, 800), fontsize = 24)
ax4 = Axis(
    fig[1, 1],
    xlabel = "",
    ylabel = "SWᵤ (W/m²)",
    xticks = 1:12,
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax4,
    collect(1:12),
    SW_u_model_monthly,
    label = "Model",
    color = "blue",
    linewidth = 3,
)
lines!(
    ax4,
    collect(1:12),
    SW_u_data_monthly,
    label = "Data",
    color = "orange",
    linewidth = 3,
)
@show mean(abs.(SW_u_model_monthly .- SW_u_data_monthly))
axislegend(ax4, position = :rt, framevisible = false)
ax3 = Axis(
    fig[3, 1],
    xlabel = "",
    ylabel = "H (W/m²)",
    xticks = 1:12,
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax3,
    collect(1:12),
    SHF_model_monthly,
    label = "",
    color = "blue",
    linewidth = 3,
)
lines!(
    ax3,
    collect(1:12),
    SHF_data_monthly,
    label = "",
    color = "orange",
    linewidth = 3,
)
@show mean(abs.(SHF_model_monthly .- SHF_data_monthly))
ax2 = Axis(
    fig[2, 1],
    xlabel = "",
    ylabel = "LE (W/m²)",
    xticks = 1:12,
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax2,
    collect(1:12),
    LHF_model_monthly,
    label = "",
    color = "blue",
    linewidth = 3,
)
lines!(
    ax2,
    collect(1:12),
    LHF_data_monthly,
    label = "",
    color = "orange",
    linewidth = 3,
)
@show mean(abs.(LHF_model_monthly .- LHF_data_monthly))
ax1 = Axis(
    fig[4, 1],
    xlabel = "Month of year",
    ylabel = "GPP (μmol/m²/s)",
    xticks = (1:12, Dates.monthabbr.(1:12)),
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax1,
    collect(1:12),
    GPP_model_monthly,
    label = "",
    color = "blue",
    linewidth = 3,
)
lines!(
    ax1,
    collect(1:12),
    GPP_data_monthly,
    label = "",
    color = "orange",
    linewidth = 3,
)
@show mean(abs.(GPP_model_monthly .- GPP_data_monthly))
CairoMakie.save(joinpath(savedir, "vaira_monthly.png"), fig)

fig = Figure(size = (800, 800), fontsize = 24)
ax0 = Axis(
    fig[4, 1],
    xlabel = "Date",
    xgridvisible = false,
    ygridvisible = false,
    yticks = [270, 290, 310],
)
ylims!(ax0, 270, 320)
lines!(ax0, model_dates, T_soil1, label = "", color = "blue", linewidth = 3)

lines!(
    ax0,
    data_dates,
    comparison_data.tsoil[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 3,
)
text!(ax0, model_dates[1], 310, text = "2 cm")

ax1 = Axis(
    fig[3, 1],
    xlabel = "",
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yticks = [270, 290, 310],
)
ylims!(ax1, 270, 320)
lines!(ax1, model_dates, T_soil5, label = "", color = "blue", linewidth = 3)

lines!(
    ax1,
    data_dates,
    comparison_data.tsoil5[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 3,
)
text!(ax1, model_dates[1], 310, text = "32 cm")


ax2 = Axis(
    fig[2, 1],
    xlabel = "",
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yticks = [0, 0.2, 0.4],
)
ylims!(ax2, 0.0, 0.5)
lines!(ax2, model_dates, θ1, color = "blue", label = "", linewidth = 3)
lines!(
    ax2,
    data_dates,
    comparison_data.swc[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 3,
)
text!(ax2, model_dates[1], 0, text = "0-2 cm")

ax3 = Axis(
    fig[1, 1],
    xlabel = "",
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yticks = [0, 0.2, 0.4],
)
ylims!(ax3, 0.0, 0.5)
lines!(ax3, model_dates, θ3, color = "blue", label = "Model", linewidth = 3)
lines!(
    ax3,
    data_dates,
    comparison_data.swc3[data_id_post_spinup],
    label = "Data",
    color = "orange",
    linewidth = 3,
)
text!(ax3, model_dates[1], 0, text = "20 cm")

ylabel01 = Label(fig[1:2, 0], "Volumetric Water (m/m)", rotation = pi / 2)
ylabel23 = Label(fig[3:4, 0], "Temperature (K)", rotation = pi / 2)

CairoMakie.save(joinpath(savedir, "vaira_soil.png"), fig)
