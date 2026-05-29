# Snowy Land Benchmarks

# The code sets up and profiles the snowy land model, which
# includes soil, canopy, and snow.

# Simulation Setup
# Number of spatial elements: 180, 360, 15
# Soil depth: 50 m
# Simulation duration: 6 hours
# Timestep: 900 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every Newton iteration
# forcing update: every 3 hours
# No user callbacks or diagnostics
# See "benchmark_utils.jl" for details on what is benchmarked

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using Dates
using CSV
const FT = Float64

######################################################################
## This result is from a benchmark run on an A100 on the clima cluster
const PREVIOUS_GPU_TIME_S = 0.515
## This result is from a benchmark run with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 15.8
######################################################################

include("benchmark_utils.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.AbstractCPUDevice ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_snowyland()
    toml_dict = LP.create_toml_dict(FT)
    start_date = DateTime(2008)
    duration = Dates.Hour(6)
    stop_date = start_date + duration
    Δt = 900.0
    time_interpolation_method = LinearInterpolation(PeriodicCalendar())
    nelements = (180, 360, 15)
    earth_param_set = LP.LandParameters(toml_dict)
    domain = ClimaLand.Domains.global_box_domain(FT; nelements)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    # Forcing data
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        time_interpolation_method = time_interpolation_method,
        context,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Construct land model with all default components
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
    )


    # Set initial conditions

    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        land;
        set_ic! = ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
            land,
        ),
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

reference_time =
    device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S :
    PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_snowyland, device, reference_time, outdir)
