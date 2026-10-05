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
import ClimaLand.Simulations: LandSimulation, solve!

import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, Statistics
import ClimaLand.LandSimVis as LandSimVis
using Flux, StaticArrays, JLD2, Adapt, InteractiveUtils

ClimaComms.@import_required_backends
NeuralSnow =
    Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt).NeuralSnow;
const FT = Float64
toml_dict = LP.create_toml_dict(FT)
climaland_dir = pkgdir(ClimaLand)
prognostic_land_components = (:canopy, :snow, :soil)

site_ID = "US-Var"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

(;rooting_depth, h_canopy,) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

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


# Now we set up the canopy model, one component at a time.
# Set up radiative transfer
radiation_parameters = (;
    Ω= FT(1),
    G_Function = CLMGFunction(FT(0)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.18),
)
radiative_transfer = Canopy.BeerLambertModel{FT}(
    surface_domain,
    toml_dict;
    radiation_parameters,
)

surface_space = land_domain.space.surface;
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
RAI = maxLAI
SAI = FT(0)
height = h_canopy
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
(; ν, θ_r) = Soil.rosetta_soil_vangenuchten_parameters(land_domain.space.subsurface, FT)
soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(surface_domain, toml_dict; soil_params = (;ν, θ_r))
# Combine the components into a CanopyModel
canopy = Canopy.CanopyModel{FT}(
    surface_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    radiative_transfer,
    biomass,
    soil_moisture_stress
)
forcing= (;atmos, radiation);
density = NeuralSnow.NeuralDepthModel(toml_dict, Δt = dt)
α_snow = NeuralSnow.NeuralAlbedoModel(toml_dict, land_domain.space.surface, Δt = dt)
snow = ClimaLand.Snow.SnowModel(
    FT,
    surface_domain,
    forcing,
    toml_dict,
    dt;
    prognostic_land_components,
    density,
    α_snow,
)

# Integrated plant hydraulics, soil, and snow model
    land = LandModel{FT}(forcing, LAI, toml_dict, land_domain, dt;canopy,prognostic_land_components, snow)

set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)

# Callbacks
output_vars = [
    "sif",
    "gs",
    "gpp",
    "ct",
    "swu",
    "lwu",
    "et",
    "msf",
    "shf",
    "lhf",
    "rn",
    "swe",
    "swc",
    "tsoil",
    "si",
    "snowc",
    "snd",
    "lai",
    "ghf",
    "soilshf",
    "soillhf",
    "soilrn","snowtsfc", "snowtb", "tair"
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
    set_ic!,
    updateat = dt, # How often we want to update the drivers.
    diagnostics = diags,
)
@time solve!(simulation)

comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir =
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/$(site_ID)/out")
mkpath(savedir)
LandSimVis.make_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = output_vars,
    spinup_date = start_date + Day(20),
    comparison_data,
)
# The observed soil temperature is compared at the depth of its sensor when
# that is documented for the site; otherwise the top model layer is shown.
tsoil_depths = FluxnetSimulations.get_sensor_depths(FT, Val(site_ID_val)).tsoil
LandSimVis.make_timeseries(
    land_domain,
    diags,
    start_date;
    savedir,
    short_names = ["tsoil"],
    spinup_date = start_date + Day(20),
    comparison_data,
    depth = isnothing(tsoil_depths) ? nothing : tsoil_depths[1],
)
