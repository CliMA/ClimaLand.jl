module Simulations
using ClimaTimesteppers
using Dates
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP

struct LandSimulation
    context
    params
    model
    domain
    timestepper
    user_callbacks
    diagnostics
    required_callbacks
    callbacks
    problem
end

function GlobalLandSimulation(FT, context, start_date, t0, Δt;
                              params = LP.earth_param_set(FT),
                              domain = global_domain(FT; comms_ctx = context),
                              forcing = ClimaLand.prescribed_forcing_era5(ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;context,),
                                                                          domain.space.surface,
                                                                          start_date,
                                                                          params,
                                                                          FT;
                                                                          time_interpolation_method = LinearInterpolation(PeriodicCalendar())),
                              LAI = ClimaLand.prescribed_lai_modis(joinpath(ClimaLand.Artifacts.modis_lai_forcing_data_path(; context), "Yuan_et_al_2008_1x1.nc"),
                                                                   domain.space.surface,
                                                                   start_date;
                                                                   time_interpolation_method = LinearInterpolation(PeriodicCalendar())),
                              ic_file = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context),
                              timestepper = ClimaTimesteppers.IMEXAlgorithm(
                                  ClimaTimesteppers.ARS111(),
                                  ClimaTimesteppers.NewtonsMethod(
                                      max_iters = 3,
                                      update_j = ClimaTimesteppers.UpdateEvery(ClimaTimesteppers.NewNewtonIteration),
                                  ),
                              ),
                              user_callbacks = (
                                  ClimaLand.NaNCheckCallback(
                                      Dates.Month(1),
                                      start_date,
                                      t0,
                                      Δt;
                                      mask = ClimaLand.landsea_mask(domain),
                                  ),
                                  ClimaLand.ReportCallback(1000)
                              ),
                              diagnostics = (;output_vars = :short, average_period = :monthly, outdir = ""), # need to generalize
                              soil_model = (type = Soil.EnergyHydrology{FT},
                                            parameters = ClimaLand.default_spatially_varying_soil_parameters(
                                                domain.space.subsurface,
                                                domain.space.surface,
                                                FT,
                                            ),
                                            runoff_type = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}),
                              snow_model = (type = Snow.SnowModel, parameters = SnowParameters{FT}(Δt;earth_param_set = params)),
                              canopy_model = (; component_types = (;
                                                                   autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
                                                                   radiative_transfer = Canopy.TwoStreamModel{FT},
                                                                   photosynthesis = Canopy.FarquharModel{FT},
                                                                   conductance = Canopy.MedlynConductanceModel{FT},
                                                                   hydraulics = Canopy.PlantHydraulicsModel{FT},
                                                                   energy = Canopy.BigLeafEnergyModel{FT},
                                                                   ),
                                              parameters = ClimaLand.clm_canopy_parameters(domain.space.surface),
                                              LAI = ),
                              soilco2_model = (; type = Soil.Biogeochemistry.SoilCO2Model{FT},
                                               Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5)),
                                               parameters = Soil.Biogeochemistry.SoilCO2ModelParameters(FT))
                              )
    
    # convert times to Itime
    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)
    
    # Create land model - this is a sketch and can be improved separately
    soil_args = construct_soil_args(soil_model.type, soil_model.spatially_varying_parameters, param_set, domain, forcing)
    runoff_model = construct_soil_runoff_model(soil_model.runoff_type, soil_model.spatially_varying_parameters, param_set, domain)
    canopy_component_args = construct_canopy_component_args(canopy_model.component_types, canopy_model.spatially_varying_parameters, param_set, domain, forcing, LAI)
    canopy_model_args = construct_canopy_model_args(canopy_model.component_types, canopy_model.spatially_varying_parameters, param_set, domain)
    snow_args = construct_snow_args(snow_model.type, snow_model.parameters, param_set, domain, forcing)
    soilco2_args = construct_soilco2_args(soilco2_model.type, soilco2_model.parameters, param_set, domain, forcing)
    land_input = (
        atmos = forcing[1]
        radiation = forcing[2]
        runoff = runoff_model,
        soil_organic_carbon = soilco2_model.Csom,
    )
    land = LandModel{FT}(;
                         soilco2_type = soilco2_model.type,
                         soilco2_args = soilco2_args,
                         soil_model_type = soil_model.type,
                         soil_args = soil_args,
                         canopy_component_types = canopy_model.component_types,
                         canopy_component_args = canopy_component_args,
                         canopy_model_args = canopy_model_args,
                         snow_args = snow_args,
                         snow_model_type = snow_model.type,
                         land_args = land_input,
                         )

    # set initial conditions
    Y, p, cds = initialize(land)
    set_ic! = ClimaLand.set_ic_from_file(ic_file)
    set_ic!(Y, p, land, t0)
    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    
    # Create tendencies and jacobian update function
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_imp_tendency(land)
    jacobian! = ClimaLand.make_jacobian(land)
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
                  Wfact = jacobian!,
                  )
    
    # Create SciML ODE Problem
    problem = SciMLBase.ODEProblem(
        ClimaTimesteppers.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )
    
    # Required callbacks
    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    required_callbacks = (driver_cb,)
    
    # Diagnostics callbacks - can be generalized in the future
    if !(diagnostics isa Nothing)
        nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            domain.space.subsurface,
            diagnostics.outdir;
            start_date,
        )
        
        diags = ClimaLand.default_diagnostics(
            land,
            start_date;
            output_writer = nc_writer,
            output_vars = diagnostics.output_vars,
            average_period = diagnostics.average_period,
        )
        
        diagnostic_handler =
            ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)
        diag_cb = (ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler),)
    else
        diag_cb = (,)
    end
    
    
    # Collect all callbacks
    callbacks = SciMLBase.CallbackSet(user_callbacks..., diag_cb..., required_cb...)
    
    return LandSimulation(context,params, land, domain, timestepper, user_callbacks, diagnostics, required_callbacks, callbacks, problem)
end                 

end
