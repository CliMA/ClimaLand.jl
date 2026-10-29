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
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Diagnostics: at_depth
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil,)
site_ID = "US-Var"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
zmin = FT(-5.0)
zmax = FT(0.0)
nelements = 15
dz_tuple = FT.((0.5, 0.04))
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
    SAI,
    f_root_to_shoot,
    rooting_depth,
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
(start_date, stop_date) = FluxnetSimulations.get_data_dates(site_ID, time_offset)
#stop_date = start_date + Year(4)
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
 #   retention_parameters,
 #   composition_parameters,
)
# Read in LAI from MODIS data
surface_space = land_domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
RAI = maxLAI * f_root_to_shoot
height = h_canopy
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)


ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict; soil_params = (;ν = soil.parameters.ν, θ_r = soil.parameters.θ_r, hydrology_cm = soil.parameters.hydrology_cm))
# Combine the components into a CanopyModel
radiation_parameters = (; Ω = 1.0,  α_NIR_leaf = 0.11,  α_PAR_leaf = 0.35, G_Function = CLMGFunction(-0.3))
rad = Canopy.BeerLambertModel{FT}(ClimaLand.Domains.obtain_surface_domain(land_domain), toml_dict; radiation_parameters)
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    #biomass,
    radiative_transfer = rad,
    soil_moisture_stress
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
lake = nothing
soilco2 = nothing
land = LandModel{FT}(canopy, snow, soil, soilco2, lake);
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
function reset_soil_ic!(Y,p,t0,land)
    evaluate!(p.drivers.T, atmos.T, t0)
    @. Y.soil.ϑ_l = land.soil.parameters.ν*0.9
    @. Y.soil.θ_i = 0

    (; θ_r, ν, ρc_ds, earth_param_set) = soil.parameters
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        land.soil.parameters.ρc_ds,
        land.soil.parameters.earth_param_set,
    )
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        p.drivers.T,
        land.soil.parameters.earth_param_set,
    )
end

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
#reset_soil_ic!(simulation._integrator.u, simulation._integrator.p, ClimaUtilities.TimeManager.ITime(0, Second(1), simulation.start_date), land)

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

spinup = Month(3)
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
# Model values at the sensor depths, interpolated between layer centers
sensor_depths = FluxnetSimulations.get_sensor_depths(FT, Val(site_ID_val))
tsoil_depths = sensor_depths.tsoil
swc_depths = sensor_depths.swc
θ1 = [
    at_depth(sol.u[k].soil.ϑ_l .+ sol.u[k].soil.θ_i, swc_depths[1]) for
    k in model_idx1:1:model_idxend
]
θ3 = [
    at_depth(sol.u[k].soil.ϑ_l .+ sol.u[k].soil.θ_i, swc_depths[3]) for
    k in model_idx1:1:model_idxend
]
T_soil1 = [
    at_depth(sv.saveval[k].soil.T, tsoil_depths[1]) for
    k in model_idx1:1:model_idxend
]
T_soil5 = [
    at_depth(sv.saveval[k].soil.T, tsoil_depths[5]) for
    k in model_idx1:1:model_idxend
]
SW_u = [parent(sv.saveval[k].SW_u)[1] for k in model_idx1:1:model_idxend]
LHF = [
    parent(
        sv.saveval[k].snow.turbulent_fluxes.lhf .*  sv.saveval[k].snow.snow_cover_fraction.+
        sv.saveval[k].soil.turbulent_fluxes.lhf .* (1 .-sv.saveval[k].snow.snow_cover_fraction) .+
        sv.saveval[k].canopy.turbulent_fluxes.lhf,
    )[1] for k in model_idx1:1:model_idxend
]
SHF = [
    parent(
        sv.saveval[k].canopy.turbulent_fluxes.shf .+
        sv.saveval[k].snow.turbulent_fluxes.shf .* sv.saveval[k].snow.snow_cover_fraction .+
        sv.saveval[k].soil.turbulent_fluxes.shf .* (1 .-sv.saveval[k].snow.snow_cover_fraction),
    )[1] for k in model_idx1:1:model_idxend
]

# Ground
LHF_g = [
    parent(
        sv.saveval[k].snow.turbulent_fluxes.lhf .* sv.saveval[k].snow.snow_cover_fraction.+
        sv.saveval[k].soil.turbulent_fluxes.lhf .* (1 .-sv.saveval[k].snow.snow_cover_fraction)
    )[1] for k in model_idx1:1:model_idxend
]
SHF_g = [
    parent(
        sv.saveval[k].snow.turbulent_fluxes.shf .* sv.saveval[k].snow.snow_cover_fraction.+
        sv.saveval[k].soil.turbulent_fluxes.shf.* (1 .-sv.saveval[k].snow.snow_cover_fraction)
    )[1] for k in model_idx1:1:model_idxend
]

# Canopy
LHF_c = [
    parent(
        sv.saveval[k].canopy.turbulent_fluxes.lhf,
    )[1] for k in model_idx1:1:model_idxend
]
SHF_c = [
    parent(
        sv.saveval[k].canopy.turbulent_fluxes.shf
    )[1] for k in model_idx1:1:model_idxend
]

GPP = [
    parent(get_GPP(sv.saveval[k], land.canopy.photosynthesis))[1] for
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
LHFc_model_monthly = compute_monthly_avg(LHF_c, model_dates)
LHFg_model_monthly = compute_monthly_avg(LHF_g, model_dates)
LHF_data_monthly = compute_monthly_avg(LHF_data, data_dates)

SHF_data = comparison_data.shf[data_id_post_spinup]
SHF_model_monthly = compute_monthly_avg(SHF, model_dates)
SHFc_model_monthly = compute_monthly_avg(SHF_c, model_dates)
SHFg_model_monthly = compute_monthly_avg(SHF_g, model_dates)
SHF_data_monthly = compute_monthly_avg(SHF_data, data_dates)

LAI_monthly = compute_monthly_avg([
    parent(
        sv.saveval[k].canopy.biomass.area_index.leaf,
    )[1] for k in model_idx1:1:model_idxend
], model_dates)
msf_monthly = compute_monthly_avg([
    parent(
        sv.saveval[k].canopy.soil_moisture_stress.βm,
    )[1] for k in model_idx1:1:model_idxend
], model_dates)
fig = Figure(size = (800, 1000), fontsize = 24)
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
    SHFc_model_monthly,
    label = "",
    color = "green",
    linewidth = 3,
)
lines!(
    ax3,
    collect(1:12),
    SHFg_model_monthly,
    label = "",
    color = "brown",
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
    LHFc_model_monthly,
    label = "",
    color = "green",
    linewidth = 3,
)
lines!(
    ax2,
    collect(1:12),
    LHFg_model_monthly,
    label = "",
    color = "brown",
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

ax0 = Axis(
    fig[5, 1],
    xlabel = "Month of year",
    ylabel = "LAI (m²/m²)",
    xticks = (1:12, Dates.monthabbr.(1:12)),
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax0,
    collect(1:12),
    LAI_monthly,
    label = "",
    color = "blue",
    linewidth = 3,
)

ax0b = Axis(
    fig[6, 1],
    xlabel = "Month of year",
    ylabel = "βm",
    xticks = (1:12, Dates.monthabbr.(1:12)),
    xgridvisible = false,
    ygridvisible = false,
)
lines!(
    ax0b,
    collect(1:12),
    msf_monthly,
    label = "",
    color = "blue",
    linewidth = 3,
)
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
lines!(ax0, model_dates, T_soil1, label = "", color = "blue", linewidth = 1)
#lines!(ax0, model_dates, [parent(sv.saveval[k].drivers.T)[1] for k in model_idx1:1:model_idxend], label = "", color = "green", linewidth = 1)

lines!(
    ax0,
    data_dates,
    comparison_data.tsoil[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 1,
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
lines!(ax1, model_dates, T_soil5, label = "", color = "blue", linewidth = 1)

lines!(
    ax1,
    data_dates,
    comparison_data.tsoil5[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 1,
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
lines!(ax2, model_dates, θ1, color = "blue", label = "", linewidth = 1)
lines!(
    ax2,
    data_dates,
    comparison_data.swc[data_id_post_spinup],
    label = "",
    color = "orange",
    linewidth = 1,
)
text!(ax2, model_dates[1], 0, text = "Surface")

ax3 = Axis(
    fig[1, 1],
    xlabel = "",
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false,
    yticks = [0, 0.2, 0.4],
)
ylims!(ax3, 0.0, 0.5)
lines!(ax3, model_dates, θ3, color = "blue", label = "Model", linewidth = 1)
lines!(
    ax3,
    data_dates,
    comparison_data.swc3[data_id_post_spinup],
    label = "Data",
    color = "orange",
    linewidth = 1,
)
text!(ax3, model_dates[1], 0, text = "20 cm")

ylabel01 = Label(fig[1:2, 0], "Volumetric Water (m/m)", rotation = pi / 2)
ylabel23 = Label(fig[3:4, 0], "Temperature (K)", rotation = pi / 2)

CairoMakie.save(joinpath(savedir, "vaira_soil.png"), fig)

SW_net =  [parent(sv.saveval[k].drivers.SW_d .- sv.saveval[k].SW_u)[1] for k in model_idx1:1:model_idxend]
LW_net =  [parent(sv.saveval[k].drivers.LW_d .- sv.saveval[k].LW_u)[1] for k in model_idx1:1:model_idxend]
SW_n_canopy = [parent(sv.saveval[k].canopy.radiative_transfer.SW_n)[1] for k in model_idx1:1:model_idxend]
LW_n_canopy = [parent(sv.saveval[k].canopy.radiative_transfer.LW_n)[1] for k in model_idx1:1:model_idxend]
R_n_soil = -1 .* [parent(sv.saveval[k].soil.R_n)[1] for k in model_idx1:1:model_idxend]
id = model_dates .< stop_date .&& model_dates .> start_date 
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, model_dates[id], SW_net[id])
lines!(ax, model_dates[id], SW_n_canopy[id])
ax2 = Axis(fig[2,1])
lines!(ax2, model_dates[id], LW_net[id])
lines!(ax2, model_dates[id], LW_n_canopy[id])
ax3 = Axis(fig[3,1])
lines!(ax3, model_dates[id], (LW_net .+ SW_net)[id])
lines!(ax3, model_dates[id], (LW_n_canopy .+ SW_n_canopy)[id])
lines!(ax3, model_dates[id], R_n_soil[id])
ax4 = Axis(fig[4,1])
lines!(ax4, model_dates[id], [parent(sv.saveval[k].canopy.biomass.area_index.leaf)[1] for k in model_idx1:1:model_idxend][id])
ax5 = Axis(fig[5,1])
lines!(ax5, model_dates[id], [parent(sv.saveval[k].canopy.soil_moisture_stress.βm)[1] for k in model_idx1:1:model_idxend][id])
save(joinpath(savedir,"some_other_plots.png"), fig)


fig2 = Figure()
ax = Axis(fig2[1,1])
lines!(ax, model_dates[id], SHF_c[id])
lines!(ax, model_dates[id], SHF_g[id])
ax2 = Axis(fig2[2,1])
lines!(ax2, model_dates[id], LHF_c[id])
lines!(ax2, model_dates[id], LHF_g[id])
save(joinpath(savedir,"some_other_plots2.png"), fig2)


error = T_soil1 .- comparison_data.tsoil[data_id_post_spinup];
days = start_date .+Month.(0:1:12*4);
monthly_error = zeros(length(days))
monthly_soil_rn = zeros(length(days))
monthly_soil_bowen = zeros(length(days))
monthly_rn = zeros(length(days))
monthly_bowen = zeros(length(days))
monthly_q = zeros(length(days))
monthly_qc = zeros(length(days))
monthly_T = zeros(length(days))
monthly_Tc = zeros(length(days))
monthly_Tleaf = zeros(length(days))
monthly_lai = zeros(length(days))
T_air = [parent(sv.saveval[k].drivers.T)[1] for k in model_idx1:1:model_idxend]
q_air = [parent(sv.saveval[k].drivers.q)[1] for k in model_idx1:1:model_idxend]
T_c = [parent(sv.saveval[k].canopy.turbulent_fluxes.T_sfc)[1] for k in model_idx1:1:model_idxend]
T_leaf = [parent(sol.u[k].canopy.energy.T)[1] for k in model_idx1:1:model_idxend]
q_c = [parent(sv.saveval[k].canopy.turbulent_fluxes.q_sfc)[1] for k in model_idx1:1:model_idxend]

lai_ts  = [parent(sv.saveval[k].canopy.biomass.area_index.leaf)[1] for k in model_idx1:1:model_idxend]
monthly_swc = zeros(length(days))
for (i,month) in enumerate(days)
    q = (model_dates .< (month + Month(1))) .&& (model_dates .>= month)
    monthly_error[i] = mean(error[q])
    monthly_soil_bowen[i] = mean(LHF_g[q])/mean(SHF_g[q] .+ LHF_g[q])
    monthly_soil_rn[i] = mean(R_n_soil[q])
    monthly_bowen[i] = mean(LHF[q])/mean(SHF[q] .+ LHF[q])
    monthly_rn[i] = mean(SW_net[q] .+ LW_net[q])
    monthly_q[i] = mean(q_air[q])
    monthly_T[i] = mean(T_air[q])
    monthly_qc[i] = mean(q_c[q])
    monthly_Tc[i] = mean(T_c[q])
    monthly_lai[i] = mean(lai_ts[q])
    monthly_swc[i] = mean(θ1[q])
    monthly_Tleaf[i] = mean(T_leaf[q])
end

fig = Figure()
ax = Axis(fig[1,1], title = "Monthly mean error")
lines!(ax, days, monthly_error)
ax2 = Axis(fig[2,1], title = "LHF/(LHF+SHF)")
lines!(ax2, days, monthly_soil_bowen)
lines!(ax2, days, monthly_bowen)
ylims!(ax2, (0,1))
ax3 = Axis(fig[3,1], title = "R_n")
lines!(ax3, days, monthly_soil_rn)
lines!(ax3, days, monthly_rn)
save(joinpath(savedir,"errors.png"), fig)

fig2 = Figure()
ax = Axis(fig2[1,1], title = "q")
lines!(ax, days, monthly_q .- monthly_qc, label = "Atmos - canopy air")
ax2 = Axis(fig2[2,1], title = "T")
lines!(ax2, days, monthly_T .- monthly_Tc, label = "Atmos - Canopy Air")
ax3 = Axis(fig2[3,1], title = "LAI")
lines!(ax3, days, monthly_lai)
ax4 = Axis(fig2[4,1], title = "SWC")
lines!(ax4, days, monthly_swc)

 save(joinpath(savedir,"errors2.png"), fig2)

