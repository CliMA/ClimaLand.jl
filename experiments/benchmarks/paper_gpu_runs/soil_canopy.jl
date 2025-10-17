function soil_canopy_integrator(;
    t0 = 0.0,
    tf = 86400.0,
    dt = 450.0,
    stepper = CTS.ARS111(),
    dlat_degrees = 1.0,
    FT = Float32,
    n_vertical_elements = 15,
    diagnostics = false,
    update_drivers = false,
    info = true,
)

    context = ClimaComms.context()
    ClimaComms.init(context)
    device = ClimaComms.device()
    info && @info "Running on $device"

    n_horizontal_elements, effective_resolution, num_columns =
        resolution(; dlat_degrees)
    info &&
        @info "Running with $(n_horizontal_elements) horizontal elements ($(round(effective_resolution, sigdigits = 2)) degrees, $num_columns columns)"

    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    prognostic_land_components = (:canopy, :soil, :soilco2)
    radius = FT(6378.1e3)
    depth = FT(3.5)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (n_horizontal_elements, n_vertical_elements),
        npolynomial = 0,
        dz_tuple = FT.((1.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    start_date = DateTime(2008)
    stop_date = start_date + Second(tf - t0)
    time_interpolation_method = LinearInterpolation(Throw())
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )

    soil_forcing = (; atmos, radiation)
    soil = Soil.EnergyHydrology{FT}(
        domain,
        soil_forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
    )

    # Soil microbes model
    soil_organic_carbon =
        ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    drivers = Soil.Biogeochemistry.SoilDrivers(
        co2_prognostic_soil,
        soil_organic_carbon,
        atmos,
    )
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(domain, drivers, toml_dict)

    # Now we set up the canopy model, which mostly use defaults for:
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_domain = ClimaLand.obtain_surface_domain(domain)
    canopy_forcing = (; atmos, radiation, ground)

    # Set up plant hydraulics
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
    )

    # Combine the soil and canopy models into a single prognostic land model
    land = SoilCanopyModel{FT}(soilco2, soil, canopy)


    function set_ic!(Y, p, t, land)
        plant_ν = FT(1.44e-4)
        soil_params = land.soil.parameters
        θ_r = soil_params.θ_r
        ν = soil_params.ν
        init_soil(ν, θ_r) = θ_r + (ν - θ_r) / 2
        Y.soil.ϑ_l .= init_soil.(ν, θ_r)
        Y.soil.θ_i .= FT(0.0)
        T = FT(276.85)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                soil_params.ρc_ds,
                soil_params.earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                soil_params.earth_param_set,
            )
        Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
        Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
        atmos = land.canopy.boundary_conditions.atmos
        evaluate!(Y.canopy.energy.T, atmos.T, t0)
    end


    sim = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        diagnostics = (),
        updateat = update_drivers ? Second(dt) : Second(2 * (tf - t0)), # we still want to update drivers on init
        user_callbacks = (),
        set_ic!,
    )
    return sim
end
