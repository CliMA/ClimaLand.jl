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
using ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

site_ID = "US-MOz"
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
    n_stem,
    n_leaf,
    h_leaf,
    h_stem,
    h_canopy,
    z_0m,
    z_0b,
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
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
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
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers, toml_dict)

# Now we set up the canopy model, one component at a time.
# Set up radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
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
conductance = PModelConductance{FT}(toml_dict)

# Set up photosynthesis
is_c3 = FT(1)
photosynthesis = PModel{FT}(surface_domain, toml_dict; is_c3)

# Set up soil moisture stress using soil retention parameters
soil_moisture_stress = PiecewiseMoistureStressModel{FT}(
    land_domain,
    toml_dict;
    soil_params = retention_parameters,
)

# Set up plant hydraulics
# Read in LAI from MODIS data
surface_space = land_domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
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
    z_0m,
    z_0b,
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    soil_moisture_stress,
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
# Callbacks
output_vars = ["gpp", "shf", "lhf", "swu", "lwu", "swc", "swe", "tsoil", "sco2", "so2", "soc"]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    reduction_period = :halfhourly,
);

simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    user_callbacks = (),
    set_ic!,
    updateat = Second(dt), # How often we want to update the drivers
    diagnostics = diags,
)

@time solve!(simulation)


using Logging

io = open("logfile.txt", "w")
logger = ConsoleLogger(io)

with_logger(logger) do
    solve!(simulation)
end

close(io)


#=
short_names = [d.variable.short_name for d in simulation.diagnostics] # short_name_X_average e.g.
diag_names = [d.output_short_name for d in simulation.diagnostics] # short_name_X_average e.g.
diag_units = [d.variable.units for d in simulation.diagnostics]


model_time, co2 = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    diag_names[9],
)
model_time, o2_a = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    diag_names[10],
)
model_time, soc = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    diag_names[11],
)

dn = diag_names
unit = diag_units
sn = short_names


save_Δt = model_time[2] - model_time[1] # in seconds since the start_date. if model_time is an Itime, the epoch should be start_date
import ClimaUtilities.TimeManager: ITime, date
function time_to_date(t::ITime, start_date)
    start_date != t.epoch &&
        @warn("$(start_date) is different from the simulation time epoch.")
    return isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
end
model_dates = time_to_date.(model_time, start_date)


fig = Figure()
ax = Axis(fig[1, 1])
p = lines!(ax, model_dates, co2, label = "co2")
ax.ylabel = "co2"
axislegend()
save("co2.png", fig)

fig = Figure()
ax = Axis(fig[1, 1])
p = lines!(ax, model_dates, o2_a, label = "o2_a")
ax.ylabel = "o2_a"
axislegend()
save("o2_a.png", fig)

fig = Figure()
ax = Axis(fig[1, 1])
p = lines!(ax, model_dates, soc, label = "soc")
ax.ylabel = "soc"
axislegend()
save("soc.png", fig)

####################



  # Get diagnostic names
  diag_names = [d.output_short_name for d in simulation.diagnostics]

  # Set indices manually based on your earlier code
  co2_idx = 9
  o2a_idx = 10
  soc_idx = 11

  # Get number of layers from the field
  field = first(values(simulation.diagnostics[1].output_writer[diag_names[co2_idx]]))
  n_layers = length(parent(field)[:, 1])

  # Time conversion function
  import ClimaUtilities.TimeManager: ITime, date
  function time_to_date(t::ITime, start_date)
      start_date != t.epoch &&
          @warn("$(start_date) is different from the simulation time epoch.")
      return isnothing(t.epoch) ? start_date + t.counter * t.period : date(t)
  end

  # Extract data for all layers
  function get_all_layers(writer, diag_name, n_layers)
      times = collect(keys(writer[diag_name]))
      sort_indices = sortperm(times)
      sorted_times = times[sort_indices]

      # Get data for each layer
      layer_data = []
      for layer in 1:n_layers
          time, vals = ClimaLand.Diagnostics.diagnostic_as_vectors(
              writer, diag_name; layer = layer
          )
          push!(layer_data, vals)
      end

      return sorted_times, layer_data
  end

  # Get all data
  model_time_co2, co2_layers = get_all_layers(
      simulation.diagnostics[1].output_writer, diag_names[co2_idx], n_layers
  )
  model_time_o2a, o2a_layers = get_all_layers(
      simulation.diagnostics[1].output_writer, diag_names[o2a_idx], n_layers
  )
  model_time_soc, soc_layers = get_all_layers(
      simulation.diagnostics[1].output_writer, diag_names[soc_idx], n_layers
  )

  # Convert times to dates
  model_dates_co2 = time_to_date.(model_time_co2, start_date)
  model_dates_o2a = time_to_date.(model_time_o2a, start_date)
  model_dates_soc = time_to_date.(model_time_soc, start_date)

  # Better color scheme
  colors = reverse(cgrad(:RdYlBu, n_layers, categorical = true))

  # ============================================================================
  # TIMESERIES PLOTS
  # ============================================================================

  layer_step = max(1, n_layers ÷ 10)
  layers_to_plot = 1:layer_step:n_layers
  if !(n_layers in layers_to_plot)
      layers_to_plot = [layers_to_plot..., n_layers]
  end

  # CO2 timeseries at selected depths
  fig_co2_ts = Figure(size = (1200, 600))
  ax_co2 = Axis(fig_co2_ts[1, 1],
      xlabel = "Date",
      ylabel = "CO2 (kg C/m³)",
      title = "CO2 Time Series at Selected Depths"
  )
  for layer in layers_to_plot
      lines!(ax_co2, model_dates_co2, co2_layers[layer],
             color = colors[layer],
             linewidth = 2,
             label = layer == n_layers ? "Surface" : layer == 1 ? "Bottom" : "Layer $layer")
  end
  Legend(fig_co2_ts[1, 2], ax_co2, "Depth", framevisible = true)
  save("co2_timeseries_all_depths.png", fig_co2_ts)

  # O2_a timeseries at selected depths
  fig_o2a_ts = Figure(size = (1200, 600))
  ax_o2a = Axis(fig_o2a_ts[1, 1],
      xlabel = "Date",
      ylabel = "O2_a (fraction)",
      title = "O2_a Time Series at Selected Depths",
      yscale = log10
  )
  for layer in layers_to_plot
      valid_vals = o2a_layers[layer] .> 1e-40
      if any(valid_vals)
          lines!(ax_o2a, model_dates_o2a[valid_vals], o2a_layers[layer][valid_vals],
                 color = colors[layer],
                 linewidth = 2,
                 label = layer == n_layers ? "Surface" : layer == 1 ? "Bottom" : "Layer $layer")
      end
  end
  Legend(fig_o2a_ts[1, 2], ax_o2a, "Depth", framevisible = true)
  save("o2a_timeseries_all_depths.png", fig_o2a_ts)

  # SOC timeseries at selected depths
  fig_soc_ts = Figure(size = (1200, 600))
  ax_soc = Axis(fig_soc_ts[1, 1],
      xlabel = "Date",
      ylabel = "SOC (kg C/m³)",
      title = "SOC Time Series at Selected Depths"
  )
  for layer in layers_to_plot
      lines!(ax_soc, model_dates_soc, soc_layers[layer],
             color = colors[layer],
             linewidth = 2,
             label = layer == n_layers ? "Surface" : layer == 1 ? "Bottom" : "Layer $layer")
  end
  Legend(fig_soc_ts[1, 2], ax_soc, "Depth", framevisible = true)
  save("soc_timeseries_all_depths.png", fig_soc_ts)

  # ============================================================================
  # VERTICAL PROFILE PLOTS - WITH TEMPORAL EVOLUTION (10 snapshots)
  # ============================================================================

  # Get depth coordinates
  z_profile = parent(land.soil.domain.fields.z)[:, 1]

  # Select 10 time snapshots evenly distributed
  n_times = length(model_time_co2)
  time_indices = unique(round.(Int, range(1, n_times, length=10)))
  n_snapshots = length(time_indices)

  # Color gradient for time evolution (early = blue, late = red)
  time_colors = cgrad(:viridis, n_snapshots, categorical = true)

  # CO2 vertical profiles over time
  fig_co2_prof = Figure(size = (900, 800))
  ax_co2_prof = Axis(fig_co2_prof[1, 1],
      xlabel = "CO2 (kg C/m³)",
      ylabel = "Depth (m)",
      title = "CO2 Vertical Profile Evolution"
  )
  for (i, tidx) in enumerate(time_indices)
      # Get profile at this time
      profile = [co2_layers[layer][tidx] for layer in 1:n_layers]
      date_label = Dates.format(model_dates_co2[tidx], "yyyy-mm-dd")
      lines!(ax_co2_prof, profile, z_profile,
             linewidth = 2,
             color = time_colors[i],
             label = date_label)
  end
  Legend(fig_co2_prof[1, 2], ax_co2_prof, "Date", framevisible = true)
  save("co2_profile_evolution.png", fig_co2_prof)

  # O2_a vertical profiles over time
  fig_o2a_prof = Figure(size = (900, 800))
  ax_o2a_prof = Axis(fig_o2a_prof[1, 1],
      xlabel = "O2_a (fraction)",
      ylabel = "Depth (m)",
      title = "O2_a Vertical Profile Evolution",
      xscale = log10
  )
  for (i, tidx) in enumerate(time_indices)
      # Get profile at this time
      profile = [o2a_layers[layer][tidx] for layer in 1:n_layers]
      # Filter for positive values
      valid_idx = profile .> 1e-40
      if any(valid_idx)
          date_label = Dates.format(model_dates_o2a[tidx], "yyyy-mm-dd")
          lines!(ax_o2a_prof, profile[valid_idx], z_profile[valid_idx],
                 linewidth = 2,
                 color = time_colors[i],
                 label = date_label)
      end
  end
  Legend(fig_o2a_prof[1, 2], ax_o2a_prof, "Date", framevisible = true)
  save("o2a_profile_evolution.png", fig_o2a_prof)

  # SOC vertical profiles over time
  fig_soc_prof = Figure(size = (900, 800))
  ax_soc_prof = Axis(fig_soc_prof[1, 1],
      xlabel = "SOC (kg C/m³)",
      ylabel = "Depth (m)",
      title = "SOC Vertical Profile Evolution"
  )
  for (i, tidx) in enumerate(time_indices)
      # Get profile at this time
      profile = [soc_layers[layer][tidx] for layer in 1:n_layers]
      date_label = Dates.format(model_dates_soc[tidx], "yyyy-mm-dd")
      lines!(ax_soc_prof, profile, z_profile,
             linewidth = 2,
             color = time_colors[i],
             label = date_label)
  end
  Legend(fig_soc_prof[1, 2], ax_soc_prof, "Date", framevisible = true)
  save("soc_profile_evolution.png", fig_soc_prof)

  # ============================================================================
  # FINAL TIMESTEP PROFILES (as before, for comparison)
  # ============================================================================

  # Get final profiles
  co2_profile = parent(simulation._integrator.u.soilco2.C)[:, 1]
  o2a_profile = parent(simulation._integrator.u.soilco2.O2_a)[:, 1]
  soc_profile = parent(simulation._integrator.u.soilco2.SOC)[:, 1]

  # CO2 vertical profile (final)
  fig_co2_prof_final = Figure(size = (700, 800))
  ax_co2_prof_final = Axis(fig_co2_prof_final[1, 1],
      xlabel = "CO2 (kg C/m³)",
      ylabel = "Depth (m)",
      title = "CO2 Vertical Profile (Final Time)"
  )
  lines!(ax_co2_prof_final, co2_profile, z_profile, linewidth = 3, color = :dodgerblue)
  scatter!(ax_co2_prof_final, co2_profile, z_profile, markersize = 12, color = :dodgerblue, strokewidth = 1, strokecolor = :black)
  save("co2_profile_final.png", fig_co2_prof_final)

  # O2_a vertical profile (final)
  fig_o2a_prof_final = Figure(size = (700, 800))
  ax_o2a_prof_final = Axis(fig_o2a_prof_final[1, 1],
      xlabel = "O2_a (fraction)",
      ylabel = "Depth (m)",
      title = "O2_a Vertical Profile (Final Time)",
      xscale = log10
  )
  valid_o2 = o2a_profile .> 1e-40
  lines!(ax_o2a_prof_final, o2a_profile[valid_o2], z_profile[valid_o2], linewidth = 3, color = :crimson)
  scatter!(ax_o2a_prof_final, o2a_profile[valid_o2], z_profile[valid_o2], markersize = 12, color = :crimson, strokewidth = 1, strokecolor =
  :black)
  save("o2a_profile_final.png", fig_o2a_prof_final)

  # SOC vertical profile (final)
  fig_soc_prof_final = Figure(size = (700, 800))
  ax_soc_prof_final = Axis(fig_soc_prof_final[1, 1],
      xlabel = "SOC (kg C/m³)",
      ylabel = "Depth (m)",
      title = "SOC Vertical Profile (Final Time)"
  )
  lines!(ax_soc_prof_final, soc_profile, z_profile, linewidth = 3, color = :forestgreen)
  scatter!(ax_soc_prof_final, soc_profile, z_profile, markersize = 12, color = :forestgreen, strokewidth = 1, strokecolor = :black)
  save("soc_profile_final.png", fig_soc_prof_final)

=#








comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/fluxnet/US-MOz/pmodel/out",
)
mkpath(savedir)
LandSimVis.make_diurnal_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["gpp", "shf", "lhf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
    comparison_data,
)
LandSimVis.make_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["swc", "tsoil", "swe"],
    spinup_date = start_date + Day(20),
    comparison_data,
)
