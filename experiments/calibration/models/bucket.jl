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

using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
    ::Type{ClimaLand.Bucket.BucketModel},
) where {FT}
    surface_space = domain.space.surface

    # Forcing data
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
        use_lowres_forcing = true,
    )

    albedo = PrescribedBaregroundAlbedo(toml_dict, surface_space)
    bucket_parameters =
        BucketModelParameters(toml_dict, albedo = albedo, τc = FT(float(Δt)))
    bucket = BucketModel(
        parameters = bucket_parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = radiation,
    )
    return bucket
end

function ClimaCalibrate.forward_model(iteration, member, ::Type{BucketModel})
    output_dir = CALIBRATE_CONFIG.output_dir
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    nelements = CALIBRATE_CONFIG.nelements
    spinup = CALIBRATE_CONFIG.spinup
    extend = CALIBRATE_CONFIG.extend
    ensemble_member_path =
        ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, member)

    eki = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    minibatch = EKP.get_current_minibatch(eki)

    # Determine start date and end date from the sample date ranges
    start_date = first(sample_date_ranges[minimum(minibatch)]) - spinup
    stop_date = last(sample_date_ranges[maximum(minibatch)]) + extend
    Δt = 900.0

    # Convert to ITimes
    t0 = ITime(0, Dates.Second(1), start_date)
    tf = ITime(
        Dates.value(convert(Dates.Second, stop_date - start_date)),
        epoch = start_date,
    )
    Δt = ITime(Δt, epoch = start_date)

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    outdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)

    depth = FT(3.5)
    dz_tuple = FT.((1.0, 0.05))
    domain =
        ClimaLand.Domains.global_domain(FT; context, nelements, depth, dz_tuple)

    calibrate_params_path =
        ClimaCalibrate.parameter_path(output_dir, iteration, member)
    toml_dict =
        LP.create_toml_dict(FT, override_files = [calibrate_params_path])

    model = setup_model(
        FT,
        start_date,
        stop_date,
        Δt,
        domain,
        toml_dict,
        ClimaLand.Bucket.BucketModel,
    )

    # Set up diagnostics
    domain = ClimaLand.get_domain(model)
    diagnostic_domain =
        haskey(domain.space, :subsurface) ? domain.space.subsurface :
        domain.space.surface
    output_writer =
        ClimaDiagnostics.NetCDFWriter(diagnostic_domain, outdir; start_date)

    # Need to include "lhf", "shf", "lwu", "swu" because plotting the
    # leaderboard requires these diagnostics
    short_names = CALIBRATE_CONFIG.short_names
    short_names = unique!([short_names; ["lhf", "shf", "lwu", "swu"]])
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars = short_names,
        reduction_period = :monthly,
        reduction_type = :average,
    )

    # IC function
    function set_ic!(Y, p, t, bucket)
        temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4
        # Set temperature IC including anomaly, based on atmospheric setup
        T_sfc_0 = 271.0
        cds = ClimaCore.Fields.coordinate_field(Y.bucket.T)
        @. Y.bucket.T = T_sfc_0 + temp_anomaly_amip(cds)
        Y.bucket.W .= 0.15
        Y.bucket.Ws .= 0.0
        Y.bucket.σS .= 0.0
    end

    # Define timestepper and ODE algorithm
    timestepper = CTS.RK4()
    timestepper = CTS.ExplicitAlgorithm(timestepper)

    simulation = LandSimulation(
        t0,
        tf,
        Δt,
        model;
        set_ic!,
        timestepper,
        outdir,
        user_callbacks = (ClimaLand.ReportCallback(div((tf - t0), 10), t0),),
        diagnostics = diagnostics,
    )
    @info "Run: Global Bucket Model"
    @info "Resolution: $nelements"
    @info "Timestep: $Δt s"
    @info "Start Date: $start_date"
    @info "Stop Date: $stop_date"
    CP.log_parameter_information(
        toml_dict,
        joinpath(ensemble_member_path, "log_params_$member.toml"),
    )
    ClimaLand.Simulations.solve!(simulation)
    return nothing
end
