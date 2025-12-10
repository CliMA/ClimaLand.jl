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
drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(soil_domain, drivers, toml_dict)

# Now we set up the canopy model, one component at a time.
# Set up radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf,
    #τ_PAR_leaf,
    α_NIR_leaf,
    #τ_NIR_leaf,
)
radiative_transfer = Canopy.BeerLambertModel{FT}(
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

# Set up optimal LAI model
lai_model =
    Canopy.OptimalLAIModel{FT}(Canopy.OptimalLAIParameters{FT}(toml_dict))

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
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    lai_model,
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
output_vars = [
    "gpp",
    "shf",
    "lhf",
    "swu",
    "lwu",
    "swc",
    "swe",
    "tsoil",
    "sco2",
    "so2",
    "soc",
    "scms",
    "lai",
    "lai_pred",
]
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

#=
using Logging

io = open("logfile.txt", "w")
logger = ConsoleLogger(io)

with_logger(logger) do
    solve!(simulation)
end

close(io)
=#

#=
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
=#

# Need to solve!
# can also look at simulation._integrator.p

short_names = [d.variable.short_name for d in simulation.diagnostics] # short_name_X_average e.g.
diag_names = [d.output_short_name for d in simulation.diagnostics] # short_name_X_average e.g.
diag_units = [d.variable.units for d in simulation.diagnostics]
i = 14
dn = diag_names[i]
unit = diag_units[i]
sn = short_names[i]
model_time, LAI_opt = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    dn,
)
LAI_opt


i = 13
dn = diag_names[i]
unit = diag_units[i]
sn = short_names[i]
model_time, LAI_obs = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    dn,
)
LAI_obs

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
p = lines!(ax, model_dates, LAI_obs, label = "LAI obs")
p2 = lines!(ax, model_dates, LAI_opt, label = "LAI opt")
ax.ylabel = "LAI m2 m-2"
axislegend()
save("test.png", fig)


# Thoughts:
# it makes no sense to have LAI < 0
# maybe params (α, z, m) are not adequate with our units of A (GPP)
# Solutions:
# we could convert A internally
# and force L to be 0 min
