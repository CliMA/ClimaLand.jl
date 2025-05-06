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
        spatially_varying_parameters = ClimaLand.default_spatially_varying_soil_parameters(
            domain.space.subsurface,
            domain.space.surface,
            FT,
        ),
        runoff_scheme = (
            type = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT},
            parameters = (;
                f_over = FT(3.28),
                R_sb = FT(1.484e-4 / 1000),
                f_max = ClimaLand.default_spatially_varying_topmodel_fmax(
                    surface_space,
                    FT,
                ),
            ),
        ),
    ),
    soilco2_model = (
        type = Soil.Biogeochemistry.SoilCO2Model{FT},
        parameters = Soil.Biogeochemistry.SoilCO2ModelParameters(FT),
        Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput((t) -> 5),
        ),
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
        component_args = (; # needs lots of work
            autotrophic_respiration = Canopy.AutotrophicRespirationParameters(
                FT,
            ),
            radiative_transfer = Canopy.TwoStreamParameters(FT;),
            photosynthesis = Canopy.FarquharParameters(
                FT,
                is_c3;
                Vcmax25 = Vcmax25,
            ),
            conductance = Canopy.MedlynConductanceParameters(FT; g1),
            hydraulics = Canopy.PlantHydraulics.PlantHydraulicsParameters(;),
            energy = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),
        ),
        parameters = Canopy.SharedCanopyParameters{FT, typeof(params)}(
            z0_m,
            z0_b,
            params,
        ),
    ),
)
    # Create land model - Sketch

    # Soil Model
    runoff_model =
        soil_model.runoff_scheme.type(; soil_model.runoff_scheme.parameters...)
    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        soil_model.spatially_varying_parameters...,
    )
    soil_args = (parameters = soil_params, domain = domain)

    # Soil CO2  - easy
    soilco2_args = (parameters = soilco2_model.parameters, domain = domain)

    # Snow - also easy given defaults
    snow_args = (
        parameters = snow_model.parameters,
        domain = ClimaLand.obtain_surface_domain(domain),
    )

    # Canopy
    canopy_model_args = (;
        parameters = canopy_model.parameters,
        domain = ClimaLand.obtain_surface_domain(domain),
    )
    land_input = (
        atmos = forcing[1],
        radiation = forcing[2],
        runoff = runoff_model,
        soil_organic_carbon = soilco2_model.Csom,
    )
    land = LandModel{FT}(;
        soilco2_type = soilco2_model.type,
        soilco2_args = soilco2_args,
        soil_model_type = soil_model.type,
        soil_args = soil_args,
        canopy_component_types = canopy_model.component_types,
        canopy_component_args = canopy_model.component_args,
        canopy_model_args = canopy_model_args,
        snow_args = snow_args,
        snow_model_type = snow_model.type,
        land_args = land_input,
    )
    return land
end
