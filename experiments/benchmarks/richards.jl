# # Global run of RichardsModel

# The code sets up and runs RichardsModel for 7 days on a spherical domain,
# using prescribed precipitation from ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# This code runs the model multiple times and collects statistics for execution time and
# allocations
#
# When run with buildkite on clima, this code also compares with the previous best time
# saved at the bottom of this file

# Simulation Setup
# This simulation setup is taken from the richards_runoff.jl experiment.
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 7 days
# Timestep: 1800 s (30 min)
# Timestepper: ARS111
# Maximum iterations: 2
# Convergence criterion: 1e-8
# Jacobian update: Every Newton iteration
# Precipitation data update: every timestep

import SciMLBase
using Dates
using Test
import NCDatasets
import Interpolations

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaParams

import ClimaLand
import ClimaLand.Parameters
using ClimaLand:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_imp_tendency,
    make_tendency_jacobian,
    make_set_initial_cache

import Profile, ProfileCanvas

const FT = Float64;

regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

outdir = "richards_benchmark_$(device_suffix)"
!ispath(outdir) && mkpath(outdir)

function setup_prob(t0, tf, Δt; nelements = (101, 15))
    # We set up the problem in a function so that we can make multiple copies (for profiling)

    # Set up simulation domain
    soil_depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6.3781e6),
        depth = soil_depth,
        nelements = nelements,
        npolynomial = 1,
        dz_tuple = FT.((10.0, 0.1)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    ref_time = DateTime(2021)

    # Read in f_max data and land sea mask
    infile_path = ClimaLand.Artifacts.topmodel_data_path()
    f_max =
        SpaceVaryingInput(infile_path, "fmax", surface_space; regridder_type)
    mask = SpaceVaryingInput(
        infile_path,
        "landsea_mask",
        surface_space;
        regridder_type,
    )

    oceans_to_zero(field, mask) = mask > 0.5 ? field : eltype(field)(0)
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )
    soil_params_artifact_path =
        ClimaLand.Artifacts.soil_params_artifact_folder_path(; context)

    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    # We use this mask to set values of these parameters over the ocean, in order
    # to keep them in the physical range
    function mask_vg(var, value)
        if var < 1e-8
            return value
        else
            return var
        end
    end
    vg_α = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_2.5x2.5x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    vg_α .= mask_vg.(vg_α, 1e-3)
    # We use this mask to set values of this parameter over the ocean, in order
    # to keep it in the physical range
    function mask_vg_n(var, value)
        if var < 1
            return value
        else
            return var
        end
    end
    vg_n = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGn_map_gupta_etal2020_2.5x2.5x4.nc",
        ),
        "n",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    vg_n .= mask_vg_n.(vg_n, 1.001)
    vg_fields_to_hcm_field(α::FT, n::FT) where {FT} =
        ClimaLand.Soil.vanGenuchten{FT}(; @NamedTuple{α::FT, n::FT}((α, n))...)
    hydrology_cm = vg_fields_to_hcm_field.(vg_α, vg_n)

    θ_r = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "residual_map_gupta_etal2020_2.5x2.5x4.nc",
        ),
        "θ_r",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    ν = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "porosity_map_gupta_etal2020_2.5x2.5x4.nc",
        ),
        "ν",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    ν .= mask_vg.(ν, 1.0)
    K_sat = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "ksat_map_gupta_etal2020_2.5x2.5x4.nc",
        ),
        "Ksat",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)

    soil_params = ClimaLand.Soil.RichardsParameters(;
        hydrology_cm = hydrology_cm,
        ν = ν,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )

    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)

    # Below, the preprocess_func argument is used to
    # 1. Convert precipitation to be negative (as it is downwards)
    # 2. Convert accumulations over an hour to a rate per second
    ref_time = DateTime(2021)
    t_start = 0.0
    # Precipitation:
    precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "tp",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
    )
    atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)
    bottom_bc = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
    bc = (;
        top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model),
        bottom = bottom_bc,
    )
    model = ClimaLand.Soil.RichardsModel{FT}(;
        parameters = soil_params,
        domain = domain,
        boundary_conditions = bc,
        sources = (),
        lateral_flow = false,
    )

    Y, p, t = initialize(model)
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
    function hydrostatic_profile(
        lat::FT,
        z::FT,
        ν::FT,
        θ_r::FT,
        α::FT,
        n::FT,
        S_s::FT,
        fmax,
    ) where {FT}
        m = 1 - 1 / n
        zmin = FT(-50.0)
        zmax = FT(0.0)

        z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5 * (fmax - 0.35) / 0.7)
        if z > z_∇
            S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
            ϑ_l = S * (ν - θ_r) + θ_r
        else
            ϑ_l = -S_s * (z - z_∇) + ν
        end
        return FT(ϑ_l)
    end

    # Set initial state values
    Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s, f_max)
    @. Y.soil.ϑ_l = oceans_to_zero(Y.soil.ϑ_l, mask)

    # Create model update functions
    set_initial_cache! = make_set_initial_cache(model)
    exp_tendency! = make_exp_tendency(model)
    imp_tendency! = make_imp_tendency(model)
    tendency_jacobian! = make_tendency_jacobian(model)

    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs = (;
        jac_prototype = ClimaLand.Soil.ImplicitEquationJacobian(Y),
        Wfact = tendency_jacobian!,
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
    updatefunc = ClimaLand.make_update_drivers(atmos, nothing)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    cb = SciMLBase.CallbackSet(driver_cb)

    return prob, cb
end

function setup_and_solve_problem(; greet = false)
    # We profile the setup phase as well here. This is not intended, but it is the easiest
    # to set up for both CPU/GPU at the same time
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
    norm_condition = CTS.MaximumError(FT(1e-8))
    conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 2,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
            convergence_checker = conv_checker,
        ),
    )

    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)
    return nothing
end

# Warm up and greet
setup_and_solve_problem(; greet = true)

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

Profile.Allocs.@profile sample_rate = 0.01 setup_and_solve_problem()
results = Profile.Allocs.fetch()
profile = ProfileCanvas.view_allocs(results)
alloc_flame_file = joinpath(outdir, "alloc_flame_$device_suffix.html")
ProfileCanvas.html_file(alloc_flame_file, profile)
@info "Saved allocation flame to $alloc_flame_file"

if ClimaComms.device() isa ClimaComms.CUDADevice
    import CUDA
    CUDA.@profile setup_and_solve_problem()
end
#=
if get(ENV, "BUILDKITE_PIPELINE_SLUG", nothing) == "climaland-benchmark"
    PREVIOUS_BEST_TIME = 9.3
    if average_timing_s > PREVIOUS_BEST_TIME + std_timing_s
        @info "Possible performance regression, previous average time was $(PREVIOUS_BEST_TIME)"
    elseif average_timing_s < PREVIOUS_BEST_TIME - std_timing_s
        @info "Possible significant performance improvement, please update PREVIOUS_BEST_TIME in $(@__DIR__)"
    end
    @testset "Performance" begin
        @test PREVIOUS_BEST_TIME - std_timing_s <=
              average_timing_s ≤
              PREVIOUS_BEST_TIME + std_timing_s
    end
end
=#
