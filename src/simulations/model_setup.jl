function LandModel(
    domain,
    start_date,
    params,
    FT;
    forcing = ClimaLand.prescribed_forcing_era5(
        ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context),# get context from domain
        domain.space.surface,
        start_date,
        params,
        FT;
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    ),
    LAI = ClimaLand.prescribed_lai_modis(
        joinpath(
            ClimaLand.Artifacts.modis_lai_forcing_data_path(; context),
            "Yuan_et_al_2008_1x1.nc",
        ),
        domain.space.surface,
        start_date;
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    ),
    soil_model = (
        type = Soil.EnergyHydrology{FT},
        parameters = ClimaLand.default_spatially_varying_soil_parameters(
            domain.space.subsurface,
            domain.space.surface,
            FT,
        ),
        runoff_type = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT},
    ),
    snow_model = (
        type = Snow.SnowModel,
        parameters = SnowParameters{FT}(Δt; earth_param_set = params),
    ),
    canopy_model = (;
        component_types = (;
            autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
            radiative_transfer = Canopy.TwoStreamModel{FT},
            photosynthesis = Canopy.FarquharModel{FT},
            conductance = Canopy.MedlynConductanceModel{FT},
            hydraulics = Canopy.PlantHydraulicsModel{FT},
            energy = Canopy.BigLeafEnergyModel{FT},
        ),
        parameters = ClimaLand.clm_canopy_parameters(domain.space.surface),
    ),
    soilco2_model = (;
        type = Soil.Biogeochemistry.SoilCO2Model{FT},
        Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput((t) -> 5),
        ),
        parameters = Soil.Biogeochemistry.SoilCO2ModelParameters(FT),
    ),
)
    # Create land model - this is a sketch
    soil_args = construct_soil_args(
        soil_model.type,
        soil_model.spatially_varying_parameters,
        params,
        domain,
        forcing,
    )
    runoff_model = construct_soil_runoff_model(
        soil_model.runoff_type,
        soil_model.spatially_varying_parameters,
        params,
        domain,
    )
    canopy_component_args = construct_canopy_component_args(
        canopy_model.component_types,
        canopy_model.spatially_varying_parameters,
        params,
        domain,
        forcing,
        LAI,
    )
    canopy_model_args = construct_canopy_model_args(
        canopy_model.component_types,
        canopy_model.spatially_varying_parameters,
        params,
        domain,
    )
    snow_args = construct_snow_args(
        snow_model.type,
        snow_model.parameters,
        params,
        domain,
        forcing,
    )
    soilco2_args = construct_soilco2_args(
        soilco2_model.type,
        soilco2_model.parameters,
        params,
        domain,
        forcing,
    )
    land_input = (
        atmos = forcing[1],
        radiation = forcing[2],
        runoff = runoff_model,
        soil_organic_carbon = soilco2_model.Csom,
    )
    land = LandModel{FT}(; # should this not be hardcoded?
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
    return land
end
