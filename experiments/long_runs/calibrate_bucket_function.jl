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
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
using Insolation

using ClimaDiagnostics
using ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
using Dates
import NCDatasets

const FT = Float64;
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "bucket_longrun_$(device_suffix)"

function setup_prob(t0, tf, Δt,params, outdir; nelements = (101, 7))
    time_interpolation_method = LinearInterpolation(PeriodicCalendar())
    regridder_type = :InterpolationsRegridder
    earth_param_set = LP.LandParameters(FT)
    radius = FT(6378.1e3)
    depth = FT(3.5)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = nelements,
        npolynomial = 1,
        dz_tuple = FT.((1.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2021)
    # Forcing data
    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)    # Precipitation:
    precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "rf",
        surface_space;
        reference_date = start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
        method = time_interpolation_method,
    )

    snow_precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "sf",
        surface_space;
        reference_date = start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
        method = time_interpolation_method,
    )

    u_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "ws",
        surface_space;
        reference_date = start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    q_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "q",
        surface_space;
        reference_date = start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    P_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "sp",
        surface_space;
        reference_date = start_date,
        regridder_type,
        method = time_interpolation_method,
    )

    T_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "t2m",
        surface_space;
        reference_date = start_date,
        regridder_type,
        method = time_interpolation_method,
    )
    h_atmos = FT(10)

    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        earth_param_set,
    )

    # Prescribed radiation
    SW_d = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "ssrd",
        surface_space;
        reference_date = start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
        method = time_interpolation_method,
    )
    LW_d = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "strd",
        surface_space;
        reference_date = start_date,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
        method = time_interpolation_method,
    )

    function zenith_angle(
        t,
        start_date;
        latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
        longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
        insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
    ) where {FT}
        # This should be time in UTC
        current_datetime = start_date + Dates.Second(round(t))

        # Orbital Data uses Float64, so we need to convert to our sim FT
        d, δ, η_UTC =
            FT.(
                Insolation.helper_instantaneous_zenith_angle(
                    current_datetime,
                    start_date,
                    insol_params,
                )
            )

        Insolation.instantaneous_zenith_angle.(
            d,
            δ,
            η_UTC,
            longitude,
            latitude,
        ).:1
    end
    radiation =
        PrescribedRadiativeFluxes(FT, SW_d, LW_d, start_date; θs = zenith_angle)

    # Set up parameters
    z_0b = z_0m
    τc = FT(Δt)
    α_snow = FT(0.8)
    (; κ_soil, ρc_soil, f_bucket, W_f, p, z_0m) = params
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space)
    bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc, f_bucket, σS_c, p, W_f, κ_soil, ρc_soil)
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
    Y.bucket.W .= FT(f_bucket*W_f)
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

    updateat = Array(t0:(3600 * 3):tf)
    drivers = ClimaLand.get_drivers(bucket)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics

    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(subsurface_space, outdir)

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
    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb)
end

function bucket_turbulent_fluxes(params)
    t0 = 0.0
    tf = 60 * 60.0 * 24 * 365
    Δt = 900.0
    nelements = (101, 7)
    output_id = ?#
    diagnostics_outdir = joinpath(root_path, "global_diagnostics", output_id)
    outdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)
    prob, cb = setup_prob(t0, tf, Δt, params, outdir; nelements)

    timestepper = CTS.RK4()
    ode_algo = CTS.ExplicitAlgorithm(timestepper)
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)

    # read in output
    simdir = ClimaAnalysis.SimDir(outdir)
    short_names =
        ["lhf", "shf"]
    timeseries_list = [[],[]]
    # First read in a monthly average of one variable to obtain size of data
    var = get(simdir; "lhf")
    
    for short_name, timeseries in zip(short_names,timeseries_list)
        var = get(simdir; short_name)
        kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
        times = ClimaAnalysis.times(var)
        for t in times
            data = ClimaAnalysis.slice(var, time = t; kwargs...)
            push!(timeseries, data)
        end
    end
    return timeseries_list
end

# not sure about mask, output_id
