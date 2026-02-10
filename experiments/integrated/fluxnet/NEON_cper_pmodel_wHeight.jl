
# ── 1. Environment & Imports ──────────────────────────────────────────────────
#import Pkg
#Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".buildkite"))

import ClimaLand
import SciMLBase
using ClimaCore
import ClimaComms
import ClimaParams as CP
using Dates
using Insolation

using ClimaLand
using ClimaLand: LandModel
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeManager: date

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CairoMakie
CairoMakie.activate!()
using Statistics
using Logging
import Random
using CSV
using DataFrames
using TOML

# ── 2. Configuration ─────────────────────────────────────────────────────────
const FT = Float64
site_ID = "NEON-cper"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)
toml_dict = LP.create_toml_dict(FT)

# Simulation parameters
spinup_days = 20
SOC_init = FT(2.0)  # Initial soil organic carbon [kg/m²]
dt = Float64(450)   # Time step: 7.5 minutes

# Domain information and site location from FLUXNET metadata
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
spinup_date = start_date + Day(spinup_days)

(; rooting_depth, n_stem, n_leaf, h_leaf, h_stem, f_root_to_shoot, SAI) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# Spatial domain with lon/lat coordinates (enables spatial parameter lookup)
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val)) #CHECK THAT

# ── 3. Atmospheric & Radiation Forcing ───────────────────────────────────────
# Time series from NEON measurements: temperature, precipitation, wind, radiation
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

# ── 4. Load & Process NEON Observations ──────────────────────────────────────
#csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
#obs_df = CSV.read(csv_path, DataFrame)

# ── 5. Model Component Setup ─────────────────────────────────────────────────
# Prognostic components (full land model with snow and soilco2)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# Forcing NamedTuple
forcing = (; atmos, radiation)

# Full LandModel with global-style spatial parameter lookup
# The LandModel constructor automatically sets up:
# - Soil: Van Genuchten, composition, albedo, runoff from spatial data
# - Canopy: pmodel photosynthesis and conductance, soil moisture stress
# - Snow: snow model with default parameters but α_snow
# - SoilCO2: soil CO2 model (included via prognostic_land_components)

# LAI from MODIS (global-style spatial data)
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
# Get the maximum LAI at this site over the first year of the simulation
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long);
RAI = maxLAI * f_root_to_shoot

# Construct the P model manually since it is not a default
photosynthesis = PModel{FT}(land_domain, toml_dict)
conductance = PModelConductance{FT}(toml_dict)
# Use the soil moisture stress function based on soil moisture only
soil_moisture_stress =
    ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict)

ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

height = h_stem + h_leaf
biomass =
    Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)

hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    toml_dict;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
)
# Canopy model with custom components
canopy = ClimaLand.Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
    hydraulics,
    biomass,
)

# ── 7. Snow Model ────────────────────────────────────────────────────────────
# Snow model with zenith angle dependent albedo
α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
snow = Snow.SnowModel(
    FT,
    canopy_domain,
    forcing,
    toml_dict,
    dt;
    prognostic_land_components,
    α_snow,
    #scf,
)

# ── 8. Integrated Land Model ────────────────────────────────────────────────
# LandModel combines all components (soil, canopy, snow, SoilCO2)
# with spatial parameters from global datasets
land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    dt;
    prognostic_land_components,
    snow,
    canopy,
)

# Custom set_ic! that overrides SOC to SOC_init (2.0 instead of default 5.0)
base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site_ID,
    start_date,
    time_offset,
    land,
)
function custom_set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    Y.soilco2.SOC .= SOC_init  # Set SOC to 2.0 kg/m²
end

# Diagnostics: sco2_ppm at halfhourly resolution
output_writer = ClimaDiagnostics.Writers.DictWriter()

short_names_1D = [
    "sif",
    "ra",
    "gs",
    "gpp",
    "ct",
    "swu",
    "lwu",
    "er",
    "et",
    "msf",
    "shf",
    "lhf",
    "rn",
]
short_names_2D = ["swc", "tsoil", "si", "sco2", "soc", "so2", "sco2_ppm"]
output_vars = [short_names_1D..., short_names_2D...]

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    reduction_period = :halfhourly,
)

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = ClimaDiagnostics.Writers.NetCDFWriter(land_domain.space.subsurface, "/Users/evametz/Documents/PostDoc/Projekte/CliMA/Siteruns/NEON-CPER/20260210_NEON-cper_pmodel_wHeight/output/"),
    output_vars,
    reduction_period = :halfhourly,
);

# Build and run simulation
simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    land;
    set_ic! = custom_set_ic!,
    updateat = Second(dt),
    diagnostics = diags,
)
solve!(simulation)
#=

using Logging

io = open("logfile5.txt", "w")
logger = ConsoleLogger(io)

with_logger(logger) do
    solve!(simulation)
end

close(io)

comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
savedir =
    #joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/US-MOz/pft/out")
    "/Users/evametz/Documents/PostDoc/Projekte/CliMA/Siteruns/FirstTries/NEON_cper_pft/NEON-cper-withERA/pft/"
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
    short_names = ["swc", "tsoil", "gpp","swu","sco2"],
    spinup_date = start_date + Day(20),
    comparison_data,
)=#