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

# Parse site_ID from command line
if length(ARGS) < 1
    error("Usage: julia --project run_any_fluxnet.jl <SITE_ID> (e.g. US-MOz)")
end
site_ID = ARGS[1]

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# Get the default domain info for a generic 10m soil column, 20 layers
(; dz_tuple, nelements, zmin, zmax) = FluxnetSimulations.get_domain_info(FT)

# Get location info from FLUXNET metadata
(; time_offset, lat, long, atmos_h) = FluxnetSimulations.get_location(site_ID)

# Construct the ClimaLand domain to run the simulation on
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# Get PFT from CLM dominant PFT map
pft = ClimaLand.Canopy.clm_dominant_pft(land_domain.space.surface)

# Get parameters from global maps + PFT + metadata
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_hydrology_cm,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_albedo,
    Ω,
    χl,
    G_Function,
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
) = FluxnetSimulations.get_parameters(FT, site_ID, land_domain, pft)

# Set up the timestepping information for the simulation
dt = Float64(450) # 7.5 minutes

# Get data dates and forcing from the flux tower site
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
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

# Set up the soil model with CLM albedo from global maps
forcing = (; atmos, radiation)
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = soil_hydrology_cm,
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)

soil = Soil.EnergyHydrology{FT}(
    land_domain,
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
soilco2 =
    Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, toml_dict)

# Canopy model components
# Radiative transfer
radiation_parameters = (;
    Ω,
    G_Function,
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

# PModel conductance
conductance = PModelConductance{FT}(toml_dict)

# PModel photosynthesis - auto-detect C3/C4 from CLM maps
surface_space = land_domain.space.surface
clm_photo_params = ClimaLand.Canopy.clm_photosynthesis_parameters(surface_space)
is_c3 = ClimaCore.Fields.field2array(clm_photo_params.is_c3)[1]
photosynthesis = PModel{FT}(surface_domain, toml_dict; is_c3)

# Soil moisture stress using soil retention parameters
soil_moisture_stress = PiecewiseMoistureStressModel{FT}(
    land_domain,
    toml_dict;
    soil_params = retention_parameters,
)

# Plant hydraulics
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
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

# Energy model
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

# Integrated land model
land = LandModel{FT}(canopy, snow, soil, soilco2)
set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)

# Diagnostics
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
]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.DictWriter(),
    output_vars,
    reduction_period = :halfhourly,
)

simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    user_callbacks = (),
    set_ic!,
    updateat = Second(dt),
    diagnostics = diags,
)

@time solve!(simulation)

# Save comparison plots
comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/fluxnet/$(site_ID)/pmodel_generic/out",
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
