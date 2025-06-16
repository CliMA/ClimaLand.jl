import SciMLBase
import ClimaComms
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
import Interpolations
using Insolation
import TOML

using ClimaDiagnostics
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime, date
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
import ClimaLand.Parameters as LP

import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL

using Statistics
using Dates
import NCDatasets

FT = Float64
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "bucket_longrun"

function setup_prob(
    t0,
    tf,
    Δt,
    params,
    model_start,
    outdir;
    nelements = (101, 7),
)
    time_interpolation_method = LinearInterpolation()
    regridder_type = :InterpolationsRegridder
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
    era5_ncdata_path =
        ClimaLand.Artifacts.find_era5_year_paths(date(t0), date(tf); context)

    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        model_start,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )

    # Set up parameters
    p_names = collect(keys(params))
    p_values = [params[name]["value"] for name in p_names]
    params = (; zip(Symbol.(p_names), p_values)...)

    (; κ_soil, ρc_soil, f_bucket, W_f, p, z_0m) = params
    z_0b = z_0m
    τc = FT(float(Δt))
    α_snow = FT(0.8)
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space)
    bucket_parameters = BucketModelParameters(
        FT;
        albedo,
        z_0m,
        z_0b,
        τc,
        f_bucket,
        p,
        W_f,
        κ_soil,
        ρc_soil,
    )
    bucket = BucketModel(
        parameters = bucket_parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = radiation,
    )

    temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4
    Y, p, cds = ClimaLand.initialize(bucket)
    # Set temperature IC including anomaly, based on atmospheric setup
    T_sfc_0 = FT(271.0)
    @. Y.bucket.T = T_sfc_0 + temp_anomaly_amip(cds.subsurface)
    Y.bucket.W .= FT(f_bucket * W_f)
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

    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        subsurface_space,
        outdir;
        start_date = model_start,
    )

    diags = ClimaLand.default_diagnostics(
        bucket,
        model_start;
        output_writer = nc_writer,
        average_period = :monthly,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb), nc_writer
end

function CAL.forward_model(iteration, member)
    ensemble_member_path = path_to_ensemble_member(caldir, iteration, member)
    params_path = parameter_path(caldir, iteration, member)
    params = TOML.parsefile(params_path) # should load a Dict, that needs to be converted to namedtuple

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)

    model_start = start_date + Year(iteration)
    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes # hours in seconds
    days = 24hours # days in seconds
    years = 367days # years in seconds - 366 to make sure we capture at least full years
    t0 = 0.0
    tf = t0 + 2years
    Δt = 450.0

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)

    t0 = ITime(t0, epoch = model_start)
    tf = ITime(tf, epoch = model_start)
    Δt = ITime(Δt, epoch = model_start)
    t0, tf, Δt = promote(t0, tf, Δt)

    prob, cb, nc_writer =
        setup_prob(t0, tf, Δt, params, model_start, diagdir; nelements)

    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    close(nc_writer)
    return nothing
end
