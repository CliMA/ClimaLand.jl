# # Global run of RichardsModel

# The code sets up and runs RichardsModel for 7 days on a spherical domain,
# using prescribed precipitation from ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# This code also assesses performance, either via Nsight or by running the
# model multiple times and collecting statistics for execution time and allocations
# to make flame graphs. You can choose between the two
# by requestion --profiler nsight or --profiler flamegraph. If not provided,
# flamegraphs are created.

# When run with buildkite on clima, without Nisght, this code also compares with the previous best time
# saved at the bottom of this file

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
        default = "flamegraph"
    end
    return parse_args(s)
end
const FT = Float64;

regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "richards_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_prob(t0, tf, Δt; nelements = (101, 15))
    dz_tuple = (10.0, 0.1)
    domain = ClimaLand.Domains.global_domain(FT; nelements, dz_tuple)
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2008)
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

    Y, p, cds = initialize(model)
    z = ClimaCore.Fields.coordinate_field(cds.subsurface).z
    lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
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

    hydrology_cm = model.parameters.hydrology_cm
    ν = model.parameters.ν
    θ_r = model.parameters.S_s
    S_s = model.parameters.θ_r
    # Set initial state values
    vg_α = hydrology_cm.α
    vg_n = hydrology_cm.n
    Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s)
    # Create model update functions
    set_initial_cache! = make_set_initial_cache(model)
    exp_tendency! = make_exp_tendency(model)
    imp_tendency! = make_imp_tendency(model)
    jacobian! = make_jacobian(model)

    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    updateat = Array(t0:(2Δt):tf)
    drivers = ClimaLand.get_drivers(model)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb)

    return prob, cb
end

function setup_simulation(; greet = false)
    t0 = 0.0
    tf = 3600.0 * 24 * 7
    Δt = 1800.0
    nelements = (101, 15)
    if greet
        @info "Run: Global RichardsModel"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    prob, cb = setup_prob(t0, tf, Δt; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 2,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    return prob, ode_algo, Δt, cb
end

parsed_args = parse_commandline()
profiler = parsed_args["profiler"]
prob, ode_algo, Δt, cb = setup_simulation(; greet = true)
@info "Starting profiling with $profiler"
if profiler == "flamegraph"
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)

    # Stop when we profile for MAX_PROFILING_TIME_SECONDS or MAX_PROFILING_SAMPLES
    MAX_PROFILING_TIME_SECONDS = 500
    MAX_PROFILING_SAMPLES = 100
    time_now = time()
    timings_s = Float64[]
    while (time() - time_now) < MAX_PROFILING_TIME_SECONDS &&
        length(timings_s) < MAX_PROFILING_SAMPLES
        lprob, lode_algo, lΔt, lcb = setup_simulation()
        push!(
            timings_s,
            ClimaComms.@elapsed device SciMLBase.solve(
                lprob,
                lode_algo;
                dt = lΔt,
                callback = lcb,
            )
        )
    end
    num_samples = length(timings_s)
    average_timing_s = round(sum(timings_s) / num_samples, sigdigits = 3)
    max_timing_s = round(maximum(timings_s), sigdigits = 3)
    min_timing_s = round(minimum(timings_s), sigdigits = 3)
    std_timing_s = round(
        sqrt(sum(((timings_s .- average_timing_s) .^ 2) / num_samples)),
        sigdigits = 3,
    )
    @info "Num samples: $num_samples"
    @info "Average time: $(average_timing_s) s"
    @info "Max time: $(max_timing_s) s"
    @info "Min time: $(min_timing_s) s"
    @info "Standard deviation time: $(std_timing_s) s"
    @info "Done profiling"

    if ClimaComms.device() isa ClimaComms.CUDADevice
        import CUDA
        lprob, lode_algo, lΔt, lcb = setup_simulation()
        p = CUDA.@profile SciMLBase.solve(
            lprob,
            lode_algo;
            dt = lΔt,
            callback = lcb,
        )
        # use "COLUMNS" to set how many horizontal characters to crop:
        # See https://github.com/ronisbr/PrettyTables.jl/issues/11#issuecomment-2145550354
        envs = ("COLUMNS" => 120,)
        withenv(envs...) do
            io = IOContext(
                stdout,
                :crop => :horizontal,
                :limit => true,
                :displaysize => displaysize(),
            )
            show(io, p)
        end
        println()
    else # Flame graphs can be misleading on GPU
        prob, ode_algo, Δt, cb = setup_simulation()
        Profile.@profile SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)
        results = Profile.fetch()
        flame_file = joinpath(outdir, "flame_$device_suffix.html")
        ProfileCanvas.html_file(flame_file, results)
        @info "Saved compute flame to $flame_file"

        prob, ode_algo, Δt, cb = setup_simulation()
        Profile.Allocs.@profile sample_rate = 0.005 SciMLBase.solve(
            prob,
            ode_algo;
            dt = Δt,
            callback = cb,
        )
        results = Profile.Allocs.fetch()
        profile = ProfileCanvas.view_allocs(results)
        alloc_flame_file = joinpath(outdir, "alloc_flame_$device_suffix.html")
        ProfileCanvas.html_file(alloc_flame_file, profile)
        @info "Saved allocation flame to $alloc_flame_file"
    end

    if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark"
        PREVIOUS_BEST_TIME = 0.6
        if average_timing_s > PREVIOUS_BEST_TIME + std_timing_s
            @info "Possible performance regression, previous average time was $(PREVIOUS_BEST_TIME)"
        elseif average_timing_s < PREVIOUS_BEST_TIME - std_timing_s
            @info "Possible significant performance improvement, please update PREVIOUS_BEST_TIME in $(@__DIR__)"
        end
        @testset "Performance" begin
            @test PREVIOUS_BEST_TIME - 2std_timing_s <=
                  average_timing_s <=
                  PREVIOUS_BEST_TIME + 2std_timing_s
        end
    end
elseif profiler == "nsight"
    integrator = SciMLBase.init(prob, ode_algo; dt = Δt, callback = cb)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
else
    @error("Profiler choice not supported.")
end
