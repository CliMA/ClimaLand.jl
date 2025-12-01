# Snowy Land Benchmarks

# The code sets up and profiles the snowy land model, which
# includes soil, canopy, and snow.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 6 hours
# Timestep: 450 s
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
## This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME_S = 0.7
## This result is from a benchmark ran with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 14.9
######################################################################

include("benchmark_utils.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_snowyland()
    toml_dict = LP.create_toml_dict(FT)
    start_date = DateTime(2008)
    duration = Dates.Hour(6)
    stop_date = start_date + duration
    Δt = 450.0
    time_interpolation_method = LinearInterpolation(PeriodicCalendar())
    nelements = (101, 15)
    earth_param_set = LP.LandParameters(toml_dict)
    domain = ClimaLand.Domains.global_domain(FT; nelements = nelements)
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
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt)


    # Set initial conditions
    function set_ic!(Y, p, t0, land)
        (; θ_r, ν, ρc_ds, earth_param_set) = land.soil.parameters
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
        Y.soil.θ_i .= 0
        T = FT(276.85)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                earth_param_set,
            )
        Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air

        plant_ν = land.canopy.hydraulics.parameters.ν
        Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
        evaluate!(Y.canopy.energy.T, atmos.T, t0)

        Y.snow.S .= 0
        Y.snow.S_l .= 0
        Y.snow.U .= 0
    end
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        land;
        set_ic!,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

reference_time =
    device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S :
    PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_snowyland, device, reference_time, outdir)
