# # Global bucket run using temporal map albedo

# The code sets up and runs the bucket for 7 days using albedo read in from a file
# containing temporally-varying data over the globe, and analytic atmospheric and radiative
# forcings.

# This code runs the bucket multiple times and collects statistics for execution time and
# allocations
#
# When run with buildkite on clima, this code also compares with the previous best time
# saved at the bottom of this file

import SciMLBase
using Dates
using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

import ClimaTimeSteppers as CTS
import NCDatasets
using ClimaCore
import ClimaComms
import ClimaLand
using ClimaParams
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedSurfaceAlbedo
using ClimaLand.Domains: coordinates, Column
using ClimaLand:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_set_initial_cache,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes

import Profile, ProfileCanvas

const FT = Float64;

context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

earth_param_set = ClimaLand.Parameters.LandParameters(FT);
outdir = "bucket_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_prob(t0, tf, Δt)
    # We set up the problem in a function so that we can make multiple copies (for profiling)

    # Set up simulation domain
    soil_depth = FT(3.5)
    bucket_domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6.3781e6),
        depth = soil_depth,
        nelements = (100, 10),
        npolynomial = 1,
        dz_tuple = FT.((1.0, 0.05)),
    )
    ref_time = DateTime(2005)

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
    albedo = PrescribedSurfaceAlbedo{FT}(ref_time, t0, surface_space)

    bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc)

    # Precipitation:
    precip = (t) -> 0
    snow_precip = (t) -> -5e-7 * (t < 1 * 86400)
    # Diurnal temperature variations:
    T_atmos = (t) -> 275.0 + 5.0 * sin(2.0 * π * t / 86400 - π / 2)
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
        ref_time,
        h_atmos,
        earth_param_set,
    )

    # Prescribed radiation -- a prescribed downwelling SW diurnal cycle, with a
    # peak at local noon, and a prescribed downwelling LW radiative
    # flux, assuming the air temperature is on average 275 degrees
    # K with a diurnal amplitude of 5 degrees K:
    SW_d = (t) -> max(1361 * sin(2π * t / 86400 - π / 2), 0.0)
    LW_d = (t) -> 5.67e-8 * (275.0 + 5.0 * sin(2.0 * π * t / 86400 - π / 2))^4
    bucket_rad = PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(SW_d),
        TimeVaryingInput(LW_d),
        ref_time,
    )


    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    Y, p, _coords = initialize(model)

    Y.bucket.T .= FT(270)
    Y.bucket.W .= FT(0.05)
    Y.bucket.Ws .= FT(0.0)
    Y.bucket.σS .= FT(0.08)

    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, t0)
    exp_tendency! = make_exp_tendency(model)
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction((T_exp!) = exp_tendency!, (dss!) = ClimaLand.dss!),
        Y,
        (t0, tf),
        p,
    )
    updateat = collect(t0:(3Δt):tf)
    updatefunc = ClimaLand.make_update_drivers(bucket_atmos, bucket_rad)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb)

    return prob, cb
end

function setup_and_solve_problem()
    # We profile the setup phase as well here. This is not intended, but it is the easiest
    # to set up for both CPU/GPU at the same time
    t0 = 0.0
    tf = 7 * 86400
    Δt = 3600.0
    prob, cb = setup_prob(t0, tf, Δt)
    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)
    return nothing
end

# Warm up
setup_and_solve_problem()

@info "Starting profiling"
# Stop when we profile for MAX_PROFILING_TIME_SECONDS or MAX_PROFILING_SAMPLES
MAX_PROFILING_TIME_SECONDS = 500
MAX_PROFILING_SAMPLES = 100
time_now = time()
timings_s = Float64[]
while (time() - time_now) < MAX_PROFILING_TIME_SECONDS &&
    length(timings_s) < MAX_PROFILING_SAMPLES
    push!(timings_s, ClimaComms.@elapsed device setup_and_solve_problem())
end
num_samples = length(timings_s)
average_timing_s = round(sum(timings_s) / num_samples, sigdigits = 3)
max_timing_s = round(maximum(timings_s), sigdigits = 3)
min_timing_s = round(minimum(timings_s), sigdigits = 3)
std_timing_s = round(
    sum(((timings_s .- average_timing_s) .^ 2) / num_samples),
    sigdigits = 3,
)
@info "Num samples: $num_samples"
@info "Average time: $(average_timing_s) s"
@info "Max time: $(max_timing_s) s"
@info "Min time: $(min_timing_s) s"
@info "Standard deviation time: $(std_timing_s) s"
@info "Done profiling"

Profile.@profile setup_and_solve_problem()
results = Profile.fetch()
flame_file = joinpath(outdir, "flame_$device_suffix.html")
ProfileCanvas.html_file(flame_file, results)
@info "Saved compute flame to $flame_file"

Profile.Allocs.@profile sample_rate = 0.1 setup_and_solve_problem()
results = Profile.Allocs.fetch()
profile = ProfileCanvas.view_allocs(results)
alloc_flame_file = joinpath(outdir, "alloc_flame_$device_suffix.html")
ProfileCanvas.html_file(alloc_flame_file, profile)
@info "Saved allocation flame to $alloc_flame_file"

if ClimaComms.device() isa ClimaComms.CUDADevice
    import CUDA
    CUDA.@profile external = true setup_and_solve_problem()
end

if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark"
    PREVIOUS_BEST_TIME = 5.4
    if average_timing_s > 1.1PREVIOUS_BEST_TIME
        @info "Possible performance regression, previous average time was $(PREVIOUS_BEST_TIME)"
    elseif average_timing_s < 0.8PREVIOUS_BEST_TIME
        @info "Possible significant performance improvement, please update PREVIOUS_BEST_TIME in $(@__DIR__)"
    end
    @testset "Performance" begin
        @test 0.8PREVIOUS_BEST_TIME <= average_timing_s ≤ 1.1PREVIOUS_BEST_TIME
    end
end
