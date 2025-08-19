# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 2-20 years, based on LONGER_RUN setting
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours

# You can choose between the profiling type
# by requestion --profiler nsight or --profiler integrated. If not provided,
# the integrated profiler is used, and if benchmarking on cpu, flamegraphs are created.

# When run with buildkite on clima, without Nisght, this code also compares with the previous best time
# saved at the top of this file

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")


import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaTimeSteppers as CTS
using ClimaUtilities.ClimaArtifacts

using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates
using Test
using ArgParse

const FT = Float64;

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--profiler"
        help = "Profiler option: nsight or flamegraph"
        arg_type = String
        default = "integrated"
    end
    return parse_args(s)
end


const PREVIOUS_GPU_TIME = 1.32

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "soil_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]

include("benchmark_sim.jl")

function setup_simulation()
    time_interpolation_method = LinearInterpolation()
    start_date = DateTime(2008)
    stop_date = start_date + Dates.Day(1)
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
    model = ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, params)
    simulation =
        LandSimulation(start_date, stop_date, Δt, model; diagnostics = [])
    return simulation
end

run_benchmarks(device, setup_simulation, profiler, PREVIOUS_GPU_TIME, outdir)
