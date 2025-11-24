# # Global bucket run using temporal map albedo
# No user callbacks or diagnostics

# The code sets up and runs the bucket for 7 days using albedo read in from a file
# containing temporally-varying data over the globe, and analytic atmospheric and radiative
# forcings.
# See "benchmark_utils.jl" for details on what is benchmarked

delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
using Dates
import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaParams as CP

import ClimaCore
import ClimaLand
using ClimaParams
import ClimaLand.Parameters as LP
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedSurfaceAlbedo
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes

const FT = Float64;

######################################################################
## This result is from a benchmark ran on an A100 on the clima cluster
const PREVIOUS_GPU_TIME_S = 0.117
## This result is from a benchmark ran with a single process on the clima cluster
const PREVIOUS_CPU_TIME_S = 0.697
######################################################################
include("benchmark_utils.jl")

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

outdir = "bucket_benchmark_$(device_suffix)"
@info "device: $device"
!ispath(outdir) && mkpath(outdir)

function setup_bucket()
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

    toml_dict = LP.create_toml_dict(FT)
    bucket_parameters = BucketModelParameters(toml_dict; albedo, z_0m, z_0b, τc)

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
        toml_dict,
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
        set_ic!,
        user_callbacks = (),
        diagnostics = [],
    )
    return simulation
end

reference_time =
    device isa ClimaComms.AbstractCPUDevice ? PREVIOUS_CPU_TIME_S :
    PREVIOUS_GPU_TIME_S
profile_and_benchmark(setup_bucket, device, reference_time, outdir)
