# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 1 day (gpu) or 2 hours (cpu)
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours
# No user callbacks or diagnostics

# See "benchmark_sim.jl" for details on what is benchmarked

# When run with buildkite on clima, without Nisght, this code compares with the previous
# best time saved at the top of this file

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")


import ClimaComms
ClimaComms.@import_required_backends
import CUDA
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime
import ClimaParams as CP
using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates

const FT = Float64;

######################################################################
# This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME = 1.32
######################################################################

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
include("benchmark_sim.jl")
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "soil_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]

function setup_simulation()
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    time_interpolation_method = LinearInterpolation()
    start_date = DateTime(2008)
    duration =
        device isa ClimaComms.CPUSingleThreaded ? Dates.Hour(2) : Dates.Day(1)
    stop_date = start_date + duration
    Δt = 450.0
    nelements = (101, 15)
    domain = ClimaLand.Domains.global_domain(FT; context, nelements)
    params = LP.LandParameters(FT)
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    forcing = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        domain.space.surface,
        start_date,
        params,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method,
    )
    model = ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, toml_dict;)
    simulation = LandSimulation(
        start_date,
        stop_date,
        Δt,
        model;
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

run_benchmarks(device, setup_simulation, profiler, PREVIOUS_GPU_TIME, outdir)
