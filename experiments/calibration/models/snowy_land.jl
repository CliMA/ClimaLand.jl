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
    ::Type{ClimaLand.LandModel};
    use_rosetta = true,
    prognostic_lai = false,
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
    prognostic_land_components = (:canopy, :lake, :snow, :soil, :soilco2)

    # Build the soil model explicitly so we can choose the van Genuchten
    # retention parameters. use_rosetta = true (the default) uses the Rosetta
    # (Montzka et al. 2017) product, matching LandModel's default soil; false
    # uses the Gupta et al. (2020) product. The chosen parameters must be paired
    # with the matching spun-up initial conditions in `forward_model`.
    retention_parameters =
        use_rosetta ?
        ClimaLand.Soil.rosetta_soil_vangenuchten_parameters(
            domain.space.subsurface,
            FT,
        ) :
        ClimaLand.Soil.soil_vangenuchten_parameters(domain.space.subsurface, FT)
    soil = ClimaLand.Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        retention_parameters,
    )

    # Feed the soil's ν/θ_r into the moisture-stress thresholds so they stay
    # consistent with the (possibly Rosetta) retention parameters.
    soil_moisture_stress = ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
        domain,
        toml_dict;
        soil_params = (; ν = soil.parameters.ν, θ_r = soil.parameters.θ_r),
    )

    if prognostic_lai
        # Prognostic LAI: the canopy computes LAI with the optimal-LAI model
        # (Zhou et al. 2025), driven by the P-model potential GPP (PModel is the
        # CanopyModel default photosynthesis). The optimal-LAI callback that
        # advances LAI at local noon is added automatically by LandSimulation
        # via get_model_callbacks. This is what makes the optimal-LAI parameters
        # (optimal_lai_z/sigma/alpha/f0/k) affect the simulated `lai`.
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            surface_domain,
            (;
                atmos,
                radiation,
                ground = ClimaLand.PrognosticGroundConditions{FT}(),
            ),
            toml_dict;
            prognostic_land_components,
            soil_moisture_stress,
            biomass = ClimaLand.Canopy.ZhouOptimalLAIModel{FT}(
                surface_domain,
                toml_dict,
            ),
        )
        land = LandModel{FT}(
            forcing,
            toml_dict,
            domain,
            Δt;
            prognostic_land_components,
            soil,
            canopy,
        )
    else
        # Prescribed LAI (default): read LAI from MODIS data.
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            surface_space,
            start_date,
            stop_date,
        )
        # Build the canopy explicitly so SAI/RAI track spatially-varying max LAI
        # (MODIS 2000-2020), which sets the AR winter maintenance baseline,
        # instead of LandModel's default constant SAI/RAI. Other canopy defaults
        # are unchanged.
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
            soil_moisture_stress,
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
            soil,
            canopy,
        )
    end
    return land
end

function ClimaCalibrate.forward_model(
    model_interface::LandModelInterface,
    iteration,
    member,
    ::Type{ClimaLand.LandModel},
)
    (; config) = model_interface
    (;
        output_dir,
        sample_date_ranges,
        nelements,
        spinup,
        extend,
        use_rosetta,
        prognostic_lai,
    ) = config
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
        ClimaLand.LandModel;
        use_rosetta,
        prognostic_lai,
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
    # With prognostic LAI, also write the optimal-LAI annual potential GPP
    # (a0a_1M_average.nc, mol CO2 m^-2 yr^-1). This is the model-derived
    # A0_annual used to refresh the `optimal_lai_inputs` artifact; it only exists
    # for the ZhouOptimalLAIModel biomass. (a0d is the daily potential GPP.)
    prognostic_lai && push!(short_names, "a0a")
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date,
        outdir;
        output_vars = short_names,
        reduction_period = :monthly,
        reduction_type = :average,
    )

    # Initial conditions must match the soil retention parameters: the spun-up
    # Rosetta ICs go with the Rosetta soil, the Gupta-based saturated ICs with
    # the Gupta soil. LandSimulation defaults to the Rosetta IC, but we set it
    # explicitly so soil and IC stay paired regardless of use_rosetta.
    ic_path =
        use_rosetta ? ClimaLand.Artifacts.rosetta_spunup_ic_path(; context) :
        ClimaLand.Artifacts.saturated_land_ic_path(; context)

    simulation = LandSimulation(
        t0,
        tf,
        Δt,
        model;
        outdir,
        set_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
            ic_path,
            model,
        ),
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
