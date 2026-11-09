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

"""Observation short names whose model diagnostic has another name."""
const OBS_ALIAS_TO_DIAGNOSTIC =
    Dict("sif_gpp" => "gpp", "inv_nee" => "nee", "res_er" => "er")

"""Calibrated factor that multiplies the default PAR and NIR reflectances of
stems, so that the two change together."""
const STEM_REFLECTANCE_SCALE = "stem_reflectance_scale"

"""
    model_parameters_path(calibrate_params_path, member_path)

Return the path of a TOML file with the model parameters of an ensemble member.
The calibrated `STEM_REFLECTANCE_SCALE` multiplies the default values of
`stem_PAR_reflectance` and `stem_NIR_reflectance`; these are written, with the
other calibrated parameters, to `model_parameters.toml` in `member_path`.
Without that scale, the calibrated parameters are the model parameters.
"""
function model_parameters_path(calibrate_params_path, member_path)
    parameters = TOML.parsefile(calibrate_params_path)
    haskey(parameters, STEM_REFLECTANCE_SCALE) || return calibrate_params_path
    scale = pop!(parameters, STEM_REFLECTANCE_SCALE)["value"]
    default_toml = LP.create_toml_dict(FT)
    for name in ("stem_PAR_reflectance", "stem_NIR_reflectance")
        parameters[name] =
            Dict("value" => scale * default_toml[name], "type" => "float")
    end
    path = joinpath(member_path, "model_parameters.toml")
    open(io -> TOML.print(io, parameters), path, "w")
    return path
end

function setup_model(
    ::Type{FT},
    start_date,
    stop_date,
    Δt,
    domain,
    toml_dict,
    ::Type{ClimaLand.LandModel},
) where {FT}
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
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
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
    model_params_path =
        model_parameters_path(calibrate_params_path, ensemble_member_path)
    toml_dict = LP.create_toml_dict(FT, override_files = [model_params_path])

    model = setup_model(
        FT,
        start_date,
        stop_date,
        Δt,
        domain,
        toml_dict,
        ClimaLand.LandModel,
    )

    # Set up diagnostics
    # Need to include "lhf", "shf", "lwu", "swu" because plotting the
    # leaderboard requires these diagnostics
    (; short_names) = config
    short_names = unique!(
        [
            [get(OBS_ALIAS_TO_DIAGNOSTIC, name, name) for name in short_names]
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

    ic_path = something(
        config.ic_path,
        ClimaLand.Artifacts.rosetta_spunup_ic_path(; context),
    )
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
    @info "Initial conditions: $ic_path"
    CP.log_parameter_information(
        toml_dict,
        joinpath(ensemble_member_path, "log_params_$member.toml"),
    )
    params = keys(TOML.parsefile(model_params_path))
    strict = true
    CP.check_override_parameter_usage(toml_dict, params, strict)
    ClimaLand.Simulations.solve!(simulation)
    return nothing
end
