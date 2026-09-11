# # Global full land (snow+soil+canopy) run

# The code sets up the ClimaLand land model  on a spherical domain,
# forcing with ERA5 data, but does not actually run the simulation.
# To run the simulation, we strongly recommend using a GPU.

# First we import a lot of packages:
import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities
import Interpolations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaParams as CP
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis;

# Set the simulation float type, determine the
# context (MPI or on a single node), and device type (CPU or GPU).
# Create a default output directory for diagnostics.
const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "land_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir);
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict);

# Set timestep, start_date, stop_date:
Δt = 450.0
start_date = DateTime(2008)
stop_date = DateTime(2009);

# Create the domain - this is intentionally low resolution,
# about 4.5 x 4.5 degrees horizontally, to avoid allocating a lot
# of memory when building the documentation.
# By default and for testing runs we use `nelements = (101, 15)`,
# which is about 0.9 x 0.9 degrees horizontally with 15 layers vertically.
nelements = (20, 7)
domain = ClimaLand.Domains.global_domain(FT; context, nelements);

# Low-resolution forcing data from ERA5 is used here,
# but high-resolution should be used for production runs.
era5_ncdata_path = ClimaLand.Artifacts.era5_land_forcing_data2008_path(;
    context,
    lowres = true,
)
forcing = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    domain.space.surface,
    start_date,
    earth_param_set,
    FT;
    max_wind_speed = 25.0,
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
);

# MODIS LAI is prescribed for the canopy model:
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);
# Make the model:
model = ClimaLand.LandModel{FT}(forcing, LAI, toml_dict, domain, Δt);
simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    outdir,
);
# Run the simulation and make plots:
#
# ``` julia
# ClimaLand.Simulations.solve!(simulation)
# LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
# LandSimVis.make_heatmaps(simulation;date = stop_date, savedir = root_path)
# LandSimVis.make_leaderboard_plots(simulation, savedir = root_path)
# ```
