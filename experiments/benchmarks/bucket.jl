# # Global bucket run using temporal map albedo

# The code sets up and runs the bucket for 7 days using albedo read in from a file
# containing temporally-varying data over the globe, and analytic atmospheric and radiative
# forcings.

# You can choose between the profiling type
# by requestion --profiler nsight or --profiler integrated. If not provided,
# the integrated profiler is used, and if benchmarking on cpu, flamegraphs are created.

# When run with buildkite on clima, without Nisght, this code also compares with the previous best time
# saved at the bottom of this file
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
using Dates
using Test
import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

import ClimaTimeSteppers as CTS
import NCDatasets
import ClimaCore
@show pkgversion(ClimaCore)
import ClimaLand
using ClimaParams
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedSurfaceAlbedo
using ClimaLand:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_set_initial_cache,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes

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


const PREVIOUS_GPU_TIME = 0.116

include("benchmark_sim.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

earth_param_set = ClimaLand.Parameters.LandParameters(FT);

parsed_args = parse_commandline()
profiler = parsed_args["profiler"]
outdir = "bucket_benchmark_$(device_suffix)"
@info "device: $device"
!ispath(outdir) && mkpath(outdir)

function setup_simulation()
    # We set up the problem in a function so that we can make multiple copies (for profiling)

    # Set up simulation domain
    dz_tuple = FT.((1.0, 0.05))
    depth = FT(3.5)
    nelements = (200, 7)
    bucket_domain =
        ClimaLand.Domains.global_domain(FT; nelements, dz_tuple, depth)
    start_date = DateTime(2005)
    end_date = start_date + Week(1)
    Δt = 3600.0

    # Initialize parameters
    σS_c = FT(0.2)
    W_f = FT(0.15)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(0.7)
    ρc_soil = FT(2e6)
    τc = FT(3600)

    surface_space = bucket_domain.space.surface
    # Construct albedo parameter object using temporal map
    albedo = PrescribedSurfaceAlbedo{FT}(start_date, surface_space)

    bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc)

    # Precipitation:
    precip = (t) -> 0
    snow_precip = (t) -> -5e-7 * (float(t) < 1 * 86400)
    # Diurnal temperature variations:
    T_atmos = (t) -> 275.0 + 5.0 * sin(2.0 * π * float(t) / 86400 - π / 2)
    # Constant otherwise:
    u_atmos = (t) -> 3.0
    q_atmos = (t) -> 0.001
    h_atmos = FT(2)
    P_atmos = (t) -> 101325

    bucket_atmos = PrescribedAtmosphere(
        TimeVaryingInput(precip),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        start_date,
        h_atmos,
        earth_param_set,
    )

    # Prescribed radiation -- a prescribed downwelling SW diurnal cycle, with a
    # peak at local noon, and a prescribed downwelling LW radiative
    # flux, assuming the air temperature is on average 275 degrees
    # K with a diurnal amplitude of 5 degrees K:
    SW_d = (t) -> max(1361 * sin(2π * float(t) / 86400 - π / 2), 0.0)
    LW_d =
        (t) ->
            5.67e-8 * (275.0 + 5.0 * sin(2.0 * π * float(t) / 86400 - π / 2))^4
    bucket_rad = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(SW_d),
        TimeVaryingInput(LW_d),
        start_date,
    )

    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )


    function set_ic!(Y, _, _, _)
        Y.bucket.T .= FT(270)
        Y.bucket.W .= FT(0.05)
        Y.bucket.Ws .= FT(0.0)
        Y.bucket.σS .= FT(0.08)
    end


    updateat = collect(start_date:Second(3Δt):end_date)
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        end_date,
        Δt,
        model;
        updateat = updateat,
        set_ic!,
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
