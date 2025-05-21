function EnergyHydrology(FT, domain, earth_param_set, forcing, prognostic_land_components;
                         runoff_model =  ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(f_over = FT(3.28), # extract from EPS
                                                                                  R_sb = FT(1.484e-4 / 1000),# extract from EPS
                                                                                  f_max = ClimaLand.topmodel_fmax(domain.space.surface,FT),
                                                                                  ),
                         retention_parameters = (; type = Soil.vanGenuchten{FT},
                                                 parameters = soil_vangenuchten_parameters(domain.space.subsurface, FT,),
                                                 ),
                         composition_parameters = soil_composition_parameters(domain.space.subsurface, FT,),
                         albedo = (; parameters = clm_soil_albedo_parameters(domain.space.surface, FT),
                                   ),
                         S_s = ClimaCore.Fields.zeros(domain.space.subsurface) .+ 1e-3,# extract from EPS or get from file
                         _top_bc = ClimaLand.AtmosDrivenFluxBC(forcing[1],
                                                              forcing[2],
                                                              runoff_model,
                                                              prognostic_land_components),
                         _bottom_bc = Soil.WaterHeatBC(; water = Soil.FreeDrainage(),
                                                       heat = Soil.HeatFluxBC((p, t) -> 0.0))
                         )
    if :canopy ∈ prognostic_land_components
        _sources = (RootExtraction{FT}(), Soil.PhaseChange{FT}())
    else
        _sources = (Soil.PhaseChange{FT}(),)
    end
    return EnergyHydrology()
end

function SnowModel(FT, domain, earth_param_set, forcing, prognostic_land_components, Δt;
                   boundary_conditions = Snow.AtmosDrivenSnowBC(forcing[1],
                                                                forcing[2],
                                                                prognostic_land_components,
                                                                ),
                   parameters = SnowParameters{FT}(Δt; earth_param_set = earth_param_set)
                   )
    return SnowModel{FT}()
end

function SoilCO2Model(FT, domain, earth_param_set, forcing, prognostic_land_components, soil_organic_carbon, soil_params;
                      top_bc = Soil.Biogeochemistry.AtmosCO2StateBC(),
                      bottom_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0),
                      parameters = Soil.Biogeochemistry.SoilCO2ModelParameters(FT),
                      sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),),
                      soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(Soil.Biogeochemistry.PrognosticMet(soil_params),
                                                                         soil_organic_carbon,
                                                                         forcing[1],
                                                                         ),
                      )
    return SoilCO2Model{FT}()
end
function CanopyModel(FT, domain, earth_param_set, scalar_params, forcing, prognostic_land_components;
                     autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT}(Canopy.AutotrophicRespirationParameters(FT)),
                     radiative_transfer = Canopy.TwoStreamModel{FT}(Canopy.TwoStreamParameters(FT, domain, scalar_params)),
                     photosynthesis = Canopy.FarquharModel{FT}(Canopy.FarquharParameters(FT, domain, scalar_params)),
                     conductance = Canopy.MedlynConductanceModel{FT}(Canopy.MedlynConductanceParameters(FT, domain, scalar_params)),
                     hydraulics = Canopy.PlantHydraulicsModel{FT}(Canopy.PlantHydraulicsParameters(scalar_params)),
                     energy = Canopy.BigLeafEnergyModel{FT}(Canopy.BigLeafEnergyParameters{FT}(scalar_params)),
                     parameters = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(scalar_params
                                                                                             earth_param_set,
                                                                                             )
                     )
    return CanopyModel()
end
# I need canopy component constructors that take in EPS and CLM path to spatially varying parameters, plus additional scalar params


function LandModel(FT, start_date, Δt, domain, earth_param_set;
                   prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
                   forcing = ClimaLand.prescribed_forcing_era5(
                       ClimaLand.Artifacts.era5_land_forcing_data2008_path(; ClimaComms.context(domain)),
                       domain.space.surface,
                       start_date,
                       earth_param_set,
                       FT;
                       time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                   ),
                   LAI = ClimaLand.prescribed_lai_modis(
                       ClimaLand.Artifacts.modis_lai_single_year_path(;
                                                                      context = ClimaComms.context(domain),
                                                                      year = Dates.year(start_date),
                                                                      ),
                       domain.space.surface,
                       start_date;
                       time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                   ),
                   soil = EnergyHydrology(FT, domain, earth_param_set, forcing, prognostic_land_components),
                   canopy = CanopyModel(FT, domain, earth_param_set, forcing, prognostic_land_components)
                   snow = SnowModel(FT, domain, earth_param_set, forcing, prognostic_land_components, Δt),
                   soilco2 = SoilCO2Model(FT, domain, earth_param_set, forcing, prognostic_land_components, soil_organic_carbon, soil_params),
                   )
    # Check that prognostic land components are the same for all components
    @assert snow.boundary_conditions.prognostic_land_components == prognostic_land_components
    @assert soil.boundary_conditions.top.prognostic_land_components == prognostic_land_components
    @assert canopy.boundary_conditions.prognostic_land_components == prognostic_land_components
    @assert progonostic_land_components == (:canopy, :snow, :soil, :soilco2)

    # Check that we are applying the correct boundary condition type, and that the forcings are the same
    @assert snow.boundary_conditions isa Snow.AtmosDrivenSnowBC
    @assert soil.boundary_conditions.top isa Soil.AtmosDrivenFluxBC
    @assert canopy.boundary_conditions isa Canopy.AtmosDrivenCanopyBC
    @assert soilco2.boundary_conditions.top isa AtmosCO2StateBC()

    @assert soil.boundary_conditions.top.atmos == forcing[1]
    @assert soil.boundary_conditions.top.radiation == forcing[2]
    @assert snow.boundary_conditions.atmos == forcing[1]
    @assert snow.boundary_conditions.radiation == forcing[2]
    @assert canopy.boundary_conditions.atmos == forcing[1]
    @assert canopy.boundary_conditions.radiation == forcing[2]
    @assert soilco2.soilco2_drivers.atmos == forcing[1]

    # Make sure the soilco2 model knows that the soil is prognostic
    @assert soilco2.soilco2_drivers.soil == Soil.Biogeochemistry.PrognosticMet(soil.parameters)

    # Make sure the canopy knows that the ground is prognostic
    @assert  canopy.boundary_conditions.ground == PrognosticGroundConditions()

    # Make sure that the soil knows that the canopy is present 
    @assert RootExtraction{FT}() ∈ soil.sources

    # Make sure all have the same earth_param_set
    @assert soil.parameters.earth_param_set == earth_param_set
    @assert snow.parameters.earth_param_set == earth_param_set
    @assert canopy.parameters.earth_param_set == earth_param_set
    @assert soilco2.parameters.earth_param_set == earth_param_set

    # Check for dt consistency:
    @assert FT(float(Δt)) == snow.parameters.Δt

    # Make sure that the LAI and other forcings have the same start date
    @assert forcing[1].start_date == start_date
    @assert forcing[2].start_date == start_date
#    @assert LAI.start_date == start_date
    
    # Make sure that the domains are consistent
    @assert soil.domain == domain
    @assert soilco2.domain == domain
    @assert snow.domain == ClimaLand.Domains.obtain_surface_domain(domain)
    @assert canopy.domain == ClimaLand.Domains.obtain_surface_domain(domain)

    land = LandModel{FT}(soilco2, soil, canopy, snow)
    return land
end
