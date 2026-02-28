# # Global bucket run

# The code sets up and runs the bucket model  on a spherical domain,
# using ERA5 data.

# First we import a lot of packages:
import ClimaComms
using ClimaCore
using ClimaUtilities
import Interpolations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
import ClimaParams as CP
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis;

# Set the simulation float type, determine the
# context (MPI or on a single node), and device type. Create
# a default output directory for diagnostics.
const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "bucket_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir);

# Set timestep, start_date, stop_date:
Δt = 900.0
start_date = DateTime(2008)
stop_date = DateTime(2009);

# Create the domain - this is intentionally low resolution,
# about 4.5 degrees x 4.5 degrees, to run quickly
# when making the documentation on CPU.
nelements = (20, 7)
depth = FT(3.5)
dz_tuple = FT.((1.0, 0.05))
domain =
    ClimaLand.Domains.global_domain(FT; context, nelements, depth, dz_tuple);

# Parameters:
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
α_snow = FT(0.8)
albedo = PrescribedBaregroundAlbedo{FT}(α_snow, domain.space.surface)
bucket_parameters = BucketModelParameters(
    toml_dict;
    albedo,
    σS_c = FT(0.2),
    W_f = FT(0.2),
    z_0m = FT(1e-3),
    z_0b = FT(1e-3),
    κ_soil = FT(1.5),
    ρc_soil = FT(2e6),
    τc = FT(float(Δt)),
);

# Low-resolution forcing data from ERA5 is used here,
# but high-resolution should be used for production runs.
era5_ncdata_path = ClimaLand.Artifacts.era5_land_forcing_data2008_path(;
    context,
    lowres = true,
)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    domain.space.surface,
    start_date,
    earth_param_set,
    FT;
    max_wind_speed = 25.0,
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
);

# Make the model:
bucket = BucketModel(
    parameters = bucket_parameters,
    domain = domain,
    atmosphere = atmos,
    radiation = radiation,
);

# Create a function which sets the initial conditions.
# This should have the argument structure (Y,p,t, model)
# in order to be used by the `LandSimulation` struct, below:
function set_ic!(Y, p, t, bucket)
    coords = ClimaCore.Fields.coordinate_field(Y.bucket.T)
    T_sfc_0 = 271.0
    @. Y.bucket.T = T_sfc_0 + 40 * cosd(coords.lat)^4
    Y.bucket.W .= 0.15
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0
end

# Define timestepper and ODE algorithm
timestepper = CTS.RK4()
timestepper = CTS.ExplicitAlgorithm(timestepper);

# Create the simulation and solve it:
simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    bucket;
    set_ic!,
    timestepper,
    outdir,
);

solve!(simulation);

# Make some plots:
short_names = ["lhf", "shf", "wsoil"]

LandSimVis.make_annual_timeseries(
    simulation;
    savedir = ".",
    short_names,
    plot_stem_name = "bucket_annual_timeseries",
)
# ![](shf_bucket_annual_timeseries.png)
# ![](lhf_bucket_annual_timeseries.png)
# ![](wsoil_bucket_annual_timeseries.png)

LandSimVis.make_heatmaps(
    simulation;
    savedir = ".",
    short_names,
    date = stop_date,
    plot_stem_name = "bucket_heatmap",
)
# ![](shf_bucket_heatmap.png)
# ![](lhf_bucket_heatmap.png)
# ![](wsoil_bucket_heatmap.png)
