# # Global bucket run

# The code sets up and runs the bucket model  on a spherical domain,
# using ERA5 data.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 5 in vertical
# Soil depth: 3.5 m
# Simulation duration: 365 d
# Timestep: 3600 s
# Timestepper: RK4
# Atmos forcing update: every 3 hours
import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime
using ClimaDiagnostics
using ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using ClimaUtilities
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
import GeoMakie
using CairoMakie
using Dates
import NCDatasets

using Poppler_jll: pdfunite

const FT = Float64;
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "bucket_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_prob(
    t0,
    tf,
    Δt,
    start_date;
    outdir = outdir,
    nelements = (101, 7),
)
    earth_param_set = LP.LandParameters(FT)
    radius = FT(6378.1e3)
    depth = FT(3.5)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = nelements,
        dz_tuple = FT.((1.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # Forcing data
    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
    era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
        regridder_type = regridder_type,
    )
    # Set up parameters
    σS_c = FT(0.2)
    W_f = FT(0.2)
    z_0m = FT(1e-3)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    τc = FT(float(Δt))
    α_snow = FT(0.8)
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space)
    bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc)
    bucket = BucketModel(
        parameters = bucket_parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = radiation,
    )

    temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4
    Y, p, cds = initialize(bucket)
    # Set temperature IC including anomaly, based on atmospheric setup
    T_sfc_0 = FT(271.0)
    @. Y.bucket.T = T_sfc_0 + temp_anomaly_amip(cds.subsurface)
    Y.bucket.W .= FT(0.15)
    Y.bucket.Ws .= FT(0.0)
    Y.bucket.σS .= FT(0.0)

    set_initial_cache! = make_set_initial_cache(bucket)
    set_initial_cache!(p, Y, t0)
    exp_tendency! = make_exp_tendency(bucket)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLand.dss!),
        Y,
        (t0, tf),
        p,
    )

    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(bucket)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics
    # num_points is the resolution of the output diagnostics
    # These are currently chosen to get a 1:1 ration with the number of
    # simulation points, ~101x101x4x4
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        subsurface_space,
        outdir;
        start_date,
        num_points = (570, 285, 50),
    )

    diags = ClimaLand.default_diagnostics(
        bucket,
        start_date;
        output_writer = nc_writer,
        average_period = :monthly,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(nancheck_freq, start_date, t0, Δt)

    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb, nancheck_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes # hours in seconds
    days = 24hours # days in seconds
    years = 366days # years in seconds - 366 to make sure we capture at least full years
    tf = 1years # 2 years in seconds
    Δt = 900.0
    start_date = DateTime(2008)
    nelements = (101, 7)
    if greet
        @info "Run: Global Bucket Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)
    prob, cb = setup_prob(t0, tf, Δt, start_date; nelements)

    # Define timestepper and ODE algorithm
    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    return nothing
end

setup_and_solve_problem(; greet = true);
# read in diagnostics and make some plots!
#### ClimaAnalysis ####
short_names = ["tsfc", "lhf", "shf", "wsoil"]
include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "long_runs",
        "figures_function.jl",
    ),
)
make_figures(root_path, outdir, short_names)
