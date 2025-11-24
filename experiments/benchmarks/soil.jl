# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 6 hours
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours
# No user callbacks or diagnostics

# See "benchmark_utils.jl" for details on what is benchmarked

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
const PREVIOUS_GPU_TIME_S = 0.326
## This result is from a benchmark ran with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 8.76
######################################################################

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
include("benchmark_utils.jl")
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "soil_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_soil()
    toml_dict = LP.create_toml_dict(FT)
    time_interpolation_method = LinearInterpolation()
    start_date = DateTime(2008)
    duration = Hour(6)
    stop_date = start_date + duration
    Δt = 450.0
    nelements = (101, 15)
    domain = ClimaLand.Domains.global_domain(FT; context, nelements)
    params = LP.LandParameters(toml_dict)
    forcing = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        domain.space.surface,
        toml_dict,
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

reference_time =
    device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S :
    PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_soil, device, reference_time, outdir)
