# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, for 6 hours on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

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
import ClimaTimeSteppers as CTS
using ClimaCore
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

const FT = Float64
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "snowy_land_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_prob(t0, tf, Δt; outdir = outdir, nelements = (101, 15))
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    other_soil_params = (; f_over, R_sb)
    
    α_snow = FT(0.67)
    scalar_snow_params = (; α_snow,Δt)
    
    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    h_leaf = FT(1.0)
    scalar_canopy_params = (; ac_canopy,K_sat_plant, a, ψ63, Weibull_param, plant_ν, plant_S_s, h_leaf)

    earth_param_set = LP.LandParameters(FT),
    
    domain = ClimaLand.global_domain(FT; nelements = nelements, context = context)
    surface_space = domain.space.surface
    start_date = DateTime(2008)
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
            context,
            lowres = true,
        )
    forcing = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )
    LAI = ClimaLand.prescribed_lai_modis(ClimaLand.Artifacts.modis_lai_forcing_data2008_path(; context),
                                         domain.space.surface,
                                         start_date)
    land = global_land_model(FT,
                             scalar_soil_params,
                             scalar_canopy_params,
                             scalar_snow_params,
                             earth_param_set;
                             context = context,
                             domain = domain,
                             forcing = forcing,
                             LAI = LAI
                             )
    
    Y, p, cds = initialize(land)

    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
    Y.soil.θ_i .= 0
    T = FT(276.85)
    ρc_s =
        Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            soil_params.ρc_ds,
            soil_params.earth_param_set,
        )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            soil_params.earth_param_set,
        )
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
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

# Warm up and greet
prob, ode_algo, Δt, cb = setup_simulation(; greet = true);
SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)

@info "Starting profiling"
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
end

if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark"
    PREVIOUS_BEST_TIME = 3.98
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
