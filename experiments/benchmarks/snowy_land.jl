# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, for 6 hours on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
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
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 6 hours
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every Newton iteration
# Atmos forcing update: every 3 hours
delete!(ENV, "JULIA_CUDA_MEMORY_POOL")

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import CUDA
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
import GeoMakie
using Dates
import Profile, ProfileCanvas
using Test
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

const FT = Float64
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_prob(t0, tf, Δt; outdir = outdir, nelements = (101, 15))
    earth_param_set = LP.LandParameters(FT)
    domain = ClimaLand.Domains.global_domain(FT; nelements = nelements)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    start_date = DateTime(2008)

    # Forcing data
    era5_ncdata_path = ClimaLand.Artifacts.era5_land_forcing_data2008_path(;
        context,
        lowres = true,
    )
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
        context = nothing,
        start_date = start_date + Second(t0),
        end_date = start_date + Second(t0) + Second(tf),
    )
    LAI = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = time_interpolation_method,
    )

    # Overwrite some defaults for the canopy model
    # Energy model
    ac_canopy = FT(2.5e3)
    energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

    # Roughness lengths
    hydraulics = Canopy.PlantHydraulicsModel{FT}(surface_domain, LAI)
    h_canopy = hydraulics.compartment_surfaces[end]
    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        earth_param_set;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        energy,
        hydraulics,
        z_0m = z0_m,
        z_0b = z0_b,
    )

    # Construct land model with all default components
    land = LandModel{FT}(forcing, LAI, earth_param_set, domain, Δt; canopy)

    Y, p, cds = initialize(land)

    # Set initial conditions
    (; θ_r, ν, ρc_ds) = land.soil.parameters
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
        Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T, earth_param_set)
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air

    plant_ν = land.canopy.hydraulics.parameters.ν
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    evaluate!(Y.canopy.energy.T, atmos.T, t0)

    Y.snow.S .= 0
    Y.snow.S_l .= 0
    Y.snow.U .= 0

    set_initial_cache! = make_set_initial_cache(land)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_imp_tendency(land)
    jacobian! = ClimaLand.make_jacobian(land)
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

    updateat = Array(t0:(3600 * 3):tf)
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    return prob, cb
end

function setup_simulation(; greet = false)
    t0 = 0.0
    tf = 60 * 60.0 * 6
    Δt = 450.0
    nelements = (101, 15)
    if greet
        @info "Run: Global Soil-Canopy-Snow Model"
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
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    return prob, ode_algo, Δt, cb
end
parsed_args = parse_commandline()
profiler = parsed_args["profiler"]
@info "Starting profiling with $profiler"
if profiler == "flamegraph"
    prob, ode_algo, Δt, cb = setup_simulation(; greet = true)
    @info "Starting profiling"
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
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
    else # Flame graphs can be misleading on gpus, so we only create them on CPU
        prob, ode_algo, Δt, cb = setup_simulation()
        Profile.@profile SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)
        results = Profile.fetch()
        flame_file = joinpath(outdir, "flame_$device_suffix.html")
        ProfileCanvas.html_file(flame_file, results)
        @info "Saved compute flame to $flame_file"

        prob, ode_algo, Δt, cb = setup_simulation()
        Profile.Allocs.@profile sample_rate = 0.0025 SciMLBase.solve(
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
    if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark" &&
       ClimaComms.device() isa ClimaComms.CUDADevice
        PREVIOUS_BEST_TIME = 0.67
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
    prob, ode_algo, Δt, cb = setup_simulation(; greet = true)
    integrator = SciMLBase.init(prob, ode_algo; dt = Δt, callback = cb)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
    SciMLBase.step!(integrator)
else
    @error("Profiler choice not supported.")
end
