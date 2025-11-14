# # Global run of RichardsModel

# The code sets up and runs RichardsModel for 7 days on a spherical domain,
# using prescribed precipitation from ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# See "benchmark_utils.jl" for details on what is benchmarked

# When run with buildkite on clima, without Nisght, this code compares with the previous
# best time saved at the top of this file

# When run with buildkite on clima, without Nisght, this code also compares with the previous best time
# saved at the top of this file
# Simulation Setup
# This simulation setup is taken from the richards_runoff.jl experiment.
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 7 days
# Timestep: 1800 s (30 min)
# Timestepper: ARS111
# Fixed number of iterations: 2
# Jacobian update: Every Newton iteration
# Precipitation data update: every timestep
# No user callbacks or diagnostics
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
using Dates
using Test
import NCDatasets

import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaTimeSteppers as CTS
using ClimaParams

import ClimaLand
import ClimaLand.Parameters
import Interpolations

import Profile, ProfileCanvas

const FT = Float64;

######################################################################
## This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME_S = 0.57
## This result is from a benchmark ran with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 1.14
######################################################################

include("benchmark_utils.jl")


context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "richards_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)


function setup_richards()
    FT = Float64
    start_date = DateTime(2008)
    duration = Dates.Week(1)
    Δt = FT(1800.0)
    stop_date = start_date + duration
    updateat = Array(start_date:Dates.Second(2Δt):stop_date)
    nelements = (101, 15)
    dz_tuple = (10.0, 0.1)
    domain = ClimaLand.Domains.global_domain(FT; nelements, dz_tuple)
    surface_space = domain.space.surface
    regridder_type = :InterpolationsRegridder
    era5_ncdata_path =
        ClimaLand.Artifacts.find_era5_year_paths(start_date, stop_date; context)
    interpolation_method = Interpolations.Linear()

    # Below, the preprocess_func argument is used to
    # 1. Convert precipitation to be negative (as it is downwards)
    # 2. Convert mass flux to equivalent liquid water flux
    # Precipitation:
    precip = TimeVaryingInput(
        era5_ncdata_path,
        "mtpr",
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 1000),
    )
    forcing = (;
        atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip),
    )
    model = ClimaLand.Soil.RichardsModel{FT}(domain, forcing)
    function hydrostatic_profile(
        lat::FT,
        z::FT,
        ν::FT,
        θ_r::FT,
        α::FT,
        n::FT,
        S_s::FT,
    ) where {FT}
        m = 1 - 1 / n
        zmin = FT(-50.0)
        zmax = FT(0.0)

        z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5)
        if z > z_∇
            S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
            ϑ_l = S * (ν - θ_r) + θ_r
        else
            ϑ_l = -S_s * (z - z_∇) + ν
        end
        return FT(ϑ_l)
    end
    function set_ic!(Y, p, t0, model)
        z = ClimaLand.get_domain(model).fields.z
        lat =
            ClimaLand.Domains.get_lat(ClimaLand.get_domain(model).space.surface)
        hydrology_cm = model.parameters.hydrology_cm
        ν = model.parameters.ν
        θ_r = model.parameters.S_s
        S_s = model.parameters.θ_r
        # Set initial state values
        vg_α = hydrology_cm.α
        vg_n = hydrology_cm.n
        Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s)
    end

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 2,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        model;
        set_ic!,
        timestepper = ode_algo,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

reference_time =
    device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S :
    PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_richards, device, reference_time, outdir)
