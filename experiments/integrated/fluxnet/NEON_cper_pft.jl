
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
toml_dict = LP.create_toml_dict(FT) #CHECK THAT


spinup_days = 20
SOC_init = FT(2.0)
dt = Float64(450)  # 7.5 minutes

# Dates - get location from FLUXNET metadata
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
spinup_date = start_date + Day(spinup_days)


# Build domain with longlat - this enables spatial parameter lookup
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val)) #CHECK THAT

# NEON site forcing
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

# LAI from MODIS (global-style)
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# ── 4. Load & Process NEON Observations ──────────────────────────────────────
#csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
#obs_df = CSV.read(csv_path, DataFrame)

# Prognostic components (full land model with snow and soilco2)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
# Forcing NamedTuple
forcing = (; atmos, radiation)

# Full LandModel with global-style spatial parameter lookup
# The LandModel constructor automatically sets up:
# - Soil: Van Genuchten, composition, albedo, runoff from spatial data
# - Canopy: default component models with spatial parameters
# - Snow: snow model with default parameters
# - SoilCO2: soil CO2 model (included via prognostic_land_components)
land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain,
    dt;
    prognostic_land_components,
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
    Y.soilco2.SOC .= SOC_init
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
output_vars = ["sco2_ppm"] #short_names_1D..., short_names_2D...]

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
    output_writer = ClimaDiagnostics.Writers.NetCDFWriter(land_domain.space.subsurface, "/Users/evametz/Documents/PostDoc/Projekte/CliMA/Siteruns/NEON-CPER/20260205_NEON-cper_new/output/"),
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