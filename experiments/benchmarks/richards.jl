# # Global run of RichardsModel

# The code sets up and runs RichardsModel for 7 days on a spherical domain,
# using prescribed precipitation from ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# You can choose between the profiling type
# by requestion --profiler nsight or --profiler integrated. If not provided,
# the integrated profiler is used, and if benchmarking on cpu, flamegraphs are created.

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
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
using Dates
using Test
import NCDatasets

import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaParams

import ClimaLand
import ClimaLand.Parameters
using ClimaLand:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_imp_tendency,
    make_jacobian,
    make_set_initial_cache

import Profile, ProfileCanvas
using ArgParse

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
const FT = Float64;

const PREVIOUS_GPU_TIME = 0.6

include("benchmark_sim.jl")

regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "richards_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]


function setup_simulation()
    FT = Float64
    start_date = DateTime(2008)
    duration = Dates.Week(1)
    Δt = FT(1800.0)
    end_date = start_date + duration
    updateat = Array(start_date:Dates.Second(2Δt):end_date)
    nelements = (101, 15)
    dz_tuple = (10.0, 0.1)
    domain = ClimaLand.Domains.global_domain(FT; nelements, dz_tuple)
    surface_space = domain.space.surface
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)

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
        end_date,
        Δt,
        model;
        updateat = updateat,
        set_ic!,
        timestepper = ode_algo,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

run_benchmarks(
    device,
    setup_simulation,
    "integrated",
    PREVIOUS_GPU_TIME,
    outdir,
)
