import TOML as toml
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON" # force climacomms to use singleton context and not MPI

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaAnalysis
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL

using Statistics
using Dates
import NCDatasets

FT = Float64;
time_interpolation_method = LinearInterpolation()
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

function setup_prob(t0, tf, Δt, params, start_date, outdir; nelements = (101, 15))
    earth_param_set = LP.LandParameters(FT)
    # Set up parameters
    p_names = collect(keys(params))
    p_values = [params[name]["value"] for name in p_names]
    params = (; zip(Symbol.(p_names), p_values)...)

    (; pc, sc, K_sat_plant, a, h_leaf, α_snow, α_soil_dry_scaler, τ_leaf_scaler, α_leaf_scaler, α_soil_scaler) = params
#    K_sat_plant = FT(5e-9) # m/s
#    a = FT(0.05 * 0.0098) # bulk modulus of elasticity
#    h_leaf = FT(1.0)
    #    α_snow = FT(0.67)
  #      α_leaf_scaler = FT(1.2) # below or above 1 OK
 #   τ_leaf_scaler = FT(1.1) # below or above 1 OK
  #  α_soil_dry_scaler = FT(1.0) # must be greater than or equal to 1 as implemented

    domain = ClimaLand.global_domain(FT; nelements = nelements)
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # Forcing data
#    era5_artifact_path =
#        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
#    era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
        era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
            date(t0),
            date(tf);
            context,
        )
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )


    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        f_max,
     ) = spatially_varying_soil_params

    # We need to ensure that the α < 1
    PAR_albedo_dry .= min.(PAR_albedo_dry .* α_soil_dry_scaler .* α_soil_scaler, FT(1))
    NIR_albedo_dry .= min.(NIR_albedo_dry .* α_soil_dry_scaler .* α_soil_scaler, FT(1))
    PAR_albedo_wet .= min.(PAR_albedo_wet .* α_soil_scaler, FT(1))
    NIR_albedo_wet .= min.(NIR_albedo_wet .* α_soil_scaler, FT(1))

    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
    )
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )

    # Spatially varying canopy parameters from CLM
    clm_parameters = ClimaLand.clm_canopy_parameters(surface_space)
    (;
        Ω,
        rooting_depth,
        is_c3,
        Vcmax25,
        g1,
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
    ) = clm_parameters

    #α <= 1
    α_PAR_leaf .= min.(α_leaf_scaler .* α_PAR_leaf, FT(1))
    α_NIR_leaf .= min.(α_leaf_scaler .* α_NIR_leaf, FT(1))

    # τ <= 1
    τ_PAR_leaf .= min.(τ_leaf_scaler .* τ_PAR_leaf, FT(1))
    τ_NIR_leaf .= min.(τ_leaf_scaler .* τ_NIR_leaf, FT(1))

    # We need to ensure that the scaling here does not violate α+τ < 1
    α_PAR_leaf .=
        ClimaLand.Canopy.enforce_albedo_constraint.(α_PAR_leaf, τ_PAR_leaf)
    α_NIR_leaf .=
        ClimaLand.Canopy.enforce_albedo_constraint.(α_NIR_leaf, τ_NIR_leaf)

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    conductivity_model =
        Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    zmax = FT(0.0)
    h_canopy = h_stem + h_leaf
    compartment_midpoints =
        n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces =
        n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    # soil microbes args
    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    # Set the soil CO2 BC to being atmospheric CO2
    soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),)

    soilco2_boundary_conditions =
        (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

    soilco2_args = (;
        boundary_conditions = soilco2_boundary_conditions,
        sources = soilco2_sources,
        domain = domain,
        parameters = soilco2_ps,
    )

    # Now we set up the canopy model, which we set up by component:
    # Component Types
    canopy_component_types = (;
        autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
        radiative_transfer = Canopy.TwoStreamModel{FT},
        photosynthesis = Canopy.FarquharModel{FT},
        conductance = Canopy.MedlynConductanceModel{FT},
        hydraulics = Canopy.PlantHydraulicsModel{FT},
        energy = Canopy.BigLeafEnergyModel{FT},
    )
    # Individual Component arguments
    # Set up autotrophic respiration
    autotrophic_respiration_args =
        (; parameters = Canopy.AutotrophicRespirationParameters(FT))
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = Canopy.TwoStreamParameters(
            FT;
            Ω,
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
            G_Function,
        )
    )
    # Set up conductance
    conductance_args =
        (; parameters = Canopy.MedlynConductanceParameters(FT; g1))
    # Set up photosynthesis
    photosynthesis_args = (;
        parameters = Canopy.FarquharParameters(
            FT,
            is_c3;
            Vcmax25 = Vcmax25,
            pc = pc,
            sc = sc,
        )
    )
    # Set up plant hydraulics
            modis_lai_ncdata_path = ClimaLand.Artifacts.find_modis_year_paths(
            date(t0),
            date(tf);
            context,
        )
#    modis_lai_artifact_path =
#    ClimaLand.Artifacts.modis_lai_forcing_data_path(; context)
#    modis_lai_ncdata_path =
#    joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")
    LAIfunction = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = time_interpolation_method,
    )
    ai_parameterization =
        Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    plant_hydraulics_args = (
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    )

    energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

    # Canopy component args
    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    # Other info needed
    shared_params = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    )

    canopy_model_args = (;
        parameters = shared_params,
        domain = ClimaLand.obtain_surface_domain(domain),
    )

    # Snow model
    snow_parameters = SnowParameters{FT}(Δt; earth_param_set = earth_param_set, α_snow = α_snow)
    snow_args = (;
        parameters = snow_parameters,
        domain = ClimaLand.obtain_surface_domain(domain),
    )
    snow_model_type = Snow.SnowModel

    land_input = (
        atmos = atmos,
        radiation = radiation,
        runoff = runoff_model,
        soil_organic_carbon = Csom,
    )
    land = LandModel{FT}(;
        soilco2_type = soilco2_type,
        soilco2_args = soilco2_args,
        land_args = land_input,
        soil_model_type = soil_model_type,
        soil_args = soil_args,
        canopy_component_types = canopy_component_types,
        canopy_component_args = canopy_component_args,
        canopy_model_args = canopy_model_args,
        snow_args = snow_args,
        snow_model_type = snow_model_type,
    )

    Y, p, cds = ClimaLand.initialize(land)

    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    evaluate!(p.snow.T, atmos.T, t0)
    ClimaLand.set_snow_initial_conditions!(
        Y,
        p,
        surface_space,
        ic_path,
        land.snow.parameters,
    )

    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
    T_bounds = extrema(Y.canopy.energy.T)

    ClimaLand.set_soil_initial_conditions!(
        Y,
        ν,
        θ_r,
        subsurface_space,
        ic_path,
        land.soil,
        T_bounds,
    )
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

    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics
    # num_points is the resolution of the output diagnostics
    # These are currently chosen to get a 1:1 ration with the number of
    # simulation points, ~101x101x4x4
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        subsurface_space,
        outdir;
        start_date,
        num_points = (4*nelements[1], 2*nelements[1], nelements[2]), # use default in `z`.
    )

    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = nc_writer,
        output_vars = :short,
        average_period = :monthly,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    walltime_info = WallTimeInfo()
    every1000steps(u, t, integrator) = mod(integrator.step, 1000) == 0
    report = let wt = walltime_info
        (integrator) -> report_walltime(wt, integrator)
    end
    report_cb = SciMLBase.DiscreteCallback(every1000steps, report)

    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(nancheck_freq, start_date, t0, Δt)

    return prob,
    SciMLBase.CallbackSet(driver_cb, diag_cb, report_cb, nancheck_cb),
    nc_writer
end

function CAL.forward_model(iteration, member)
    ensemble_member_path = path_to_ensemble_member(caldir, iteration, member)
    params_path = parameter_path(caldir, iteration, member)
    params = toml.parsefile(params_path)

    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes # hours in seconds
    days = 24hours # days in seconds
    years = 366days # years in seconds - 366 to make sure we capture at least full years
    t0 = 0.0 + iteration*1years
    tf = t0 + 2years
    Δt = 450.0

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)

    start_date = DateTime(2008)
    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)
    prob, cb, nc_writer = setup_prob(t0, tf, Δt, params, start_date, diagdir; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    close(nc_writer)
    return nothing
end
