# # Global bucket run

# The code sets up and runs the bucket model  on a spherical domain,
# using ERA5 data.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 5 in vertical
# Soil depth: 3.5 m
# Simulation duration: 365 d
# Timestep: 3600 s
# Timestepper: RK4
# Atmos forcing update: every 3 hours
import ClimaComms
using ClimaCore
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime
using ClimaUtilities
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
import ClimaLand
import ClimaLand.Parameters as LP

import ClimaLand.Simulations: LandSimulation, solve!
using Statistics
using Dates
import NCDatasets

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

const FT = Float64;
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "bucket_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_model(FT, start_date, domain, earth_param_set, Δt, toml_dict)
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method,
        regridder_type,
    )

    albedo = PrescribedBaregroundAlbedo(toml_dict, surface_space)
    bucket_parameters =
        BucketModelParameters(toml_dict, albedo = albedo, τc = FT(float(Δt)))
    bucket = BucketModel(
        parameters = bucket_parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = radiation,
    )
    return bucket
end

Δt = 900.0
start_date = DateTime(2008)
stop_date = DateTime(2010)
nelements = (101, 7)
@info "Run: Global Bucket Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
# Domain
depth = FT(3.5)
dz_tuple = FT.((1.0, 0.05))
domain =
    ClimaLand.Domains.global_domain(FT; context, nelements, depth, dz_tuple)

# Parameters
toml_dict = LP.create_toml_dict(FT)
params = LP.LandParameters(toml_dict)

# Model
model = setup_model(FT, start_date, domain, params, Δt, toml_dict)

# IC function
function set_ic!(Y, p, t, bucket)
    temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4
    # Set temperature IC including anomaly, based on atmospheric setup
    T_sfc_0 = 271.0
    cds = ClimaCore.Fields.coordinate_field(Y.bucket.T)
    @. Y.bucket.T = T_sfc_0 + temp_anomaly_amip(cds)
    Y.bucket.W .= 0.15
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0
end

# Define timestepper and ODE algorithm
timestepper = CTS.RK4()
timestepper = CTS.ExplicitAlgorithm(timestepper)

# Create the simulation
simulation = LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    set_ic!,
    timestepper,
    outdir,
)
@info "Run: Global Bucket Model"
@info "Resolution: $nelements"
@info "Timestep: $Δt s"
@info "Start Date: $start_date"
@info "Stop Date: $stop_date"
solve!(simulation);

short_names = ["tsfc", "lhf", "shf", "wsoil"]
LandSimVis.make_annual_timeseries(simulation; savedir = root_path, short_names)
LandSimVis.make_heatmaps(
    simulation;
    savedir = root_path,
    short_names,
    date = stop_date,
)
