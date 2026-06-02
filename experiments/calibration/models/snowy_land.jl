# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 730 d
# Timestep: 900 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
    ::Type{ClimaLand.LandModel},
) where {FT}
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    # Forcing data - always use high resolution for calibration runs
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        context,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)

    # Build the canopy explicitly so SAI/RAI track spatially-varying max LAI
    # (MODIS 2000-2020), which sets the AR winter maintenance baseline, instead
    # of LandModel's default constant SAI/RAI. Other canopy defaults are unchanged.
    maxLAI = ClimaLand.Canopy.modis_max_lai(surface_space)
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        (;
            atmos,
            radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components,
        soil_moisture_stress = ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
            domain,
            toml_dict,
        ),
        biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
            surface_domain,
            LAI,
            maxLAI,
            toml_dict,
        ),
    )
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        canopy,
    )
    return land
end

function ClimaCalibrate.forward_model(
    model_interface::LandModelInterface,
    iteration,
    member,
    ::Type{ClimaLand.LandModel},
)
    (; config) = model_interface
    (; output_dir, sample_date_ranges, nelements, spinup, extend) = config
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

    domain = ClimaLand.Domains.global_box_domain(
        FT;
        context,
        nelements,
        mask_threshold = FT(0.99),
    )

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
        ClimaLand.LandModel,
    )

    # Set up diagnostics. lhf/shf/lwu/swu are needed for the leaderboard.
    # Observation-alias short_names (inv_nee, ...) map to their model diagnostic
    # names here; the aliases are reconstructed in data_sources.jl.
    obs_alias_to_diag = Dict(
        "inv_nee" => "nee",
        "sif_gpp" => "gpp",
        "res_er" => "er",
        "inv_hr" => "hr",
    )
    (; short_names) = config
    short_names = unique!(
        [
            [get(obs_alias_to_diag, sn, sn) for sn in short_names]
            ["lhf", "shf", "lwu", "swu", "gpp", "et", "hr", "ra", "er", "nee"]
        ],
    )
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date,
        outdir;
        output_vars = short_names,
        reduction_period = :monthly,
        reduction_type = :average,
    )

    simulation = LandSimulation(
        t0,
        tf,
        Δt,
        model;
        outdir,
        user_callbacks = (ClimaLand.ReportCallback(div((tf - t0), 10), t0),),
        diagnostics = diagnostics,
    )
    @info "Run: Global Soil-Canopy-Snow Model"
    @info "Resolution: $nelements"
    @info "Timestep: $Δt s"
    @info "Start Date: $start_date"
    @info "Stop Date: $stop_date"
    CP.log_parameter_information(
        toml_dict,
        joinpath(ensemble_member_path, "log_params_$member.toml"),
    )
    params = keys(TOML.parsefile(calibrate_params_path))
    strict = true
    CP.check_override_parameter_usage(toml_dict, params, strict)
    ClimaLand.Simulations.solve!(simulation)
    return nothing
end
