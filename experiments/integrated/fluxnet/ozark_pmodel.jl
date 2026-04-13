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
using ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :soil, :soilco2)#

site_ID = "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

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

forcing = (; atmos, radiation)
soil = EnergyHydrology{FT}(land_domain, forcing, toml_dict; prognostic_land_components,            additional_sources = (ClimaLand.RootExtraction{FT}(),));
# Set up conductance
conductance = PModelConductance{FT}(toml_dict)

# Set up photosynthesis
fractional_c3 = FT(1)
photosynthesis = PModel{FT}(surface_domain, toml_dict; fractional_c3)

# Set up soil moisture stress using soil retention parameters
soil_moisture_stress = ExperimentalMSModel{FT}(
    land_domain,
    toml_dict;
    soil_params = soil.parameters
)

# Set up plant hydraulics
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
SAI = FT(0)
RAI = maxLAI * FT(3.5)
surface_space = land_domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
biomass =
    Canopy.PrescribedBiomassModel{FT}(land_domain, LAI, toml_dict; SAI, RAI)
ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    canopy_forcing,
    toml_dict;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
    biomass,
);

# Integrated plant hydraulics, soil, and snow model

land = SoilCanopyModel{FT}(forcing, LAI, toml_dict, land_domain; canopy, soil);
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
    "msf",
    "swc",
    "tsoil",
    "lai",
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

@time solve!(simulation);

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
    short_names = ["swc", "tsoil", "lai", "msf", "gpp", "lhf", "shf", "swu", "lwu"],
    spinup_date = start_date + Day(20),
    comparison_data,
)

