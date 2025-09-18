function snowy_land_integrator(;
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
    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (n_horizontal_elements, n_vertical_elements),
        npolynomial = 0,
        dz_tuple = FT.((10.0, 0.05)),
    )
    surface_space = domain.space.surface
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
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
        max_wind_speed = 25.0,
        time_interpolation_method,
    )
    forcing = (; atmos, radiation)

    # Read in LAI from MODIS data
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    # Overwrite some defaults for the canopy model

    # Plant hydraulics
    conductivity_model = Canopy.PlantHydraulics.Weibull(toml_dict)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve(toml_dict)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        LAI,
        toml_dict;
        conductivity_model,
        retention_model,
    )

    # Roughness lengths
    h_canopy = hydraulics.compartment_surfaces[end]
    z_0m = FT(0.13) * h_canopy
    z_0b = FT(0.1) * z_0m

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        hydraulics,
        z_0m,
        z_0b,
    )

    # Snow model setup
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict)
    horz_degree_res = FT(
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2,
    ) # mean of resolution in latitude and longitude, in degrees
    scf = Snow.WuWuSnowCoverFractionModel(toml_dict, horz_degree_res)
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        dt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
        α_snow,
        scf,
    )

    # Construct the land model with all default components except for snow
    land = LandModel{FT}(forcing, LAI, toml_dict, domain, dt; snow, canopy)
    function set_ic!(Y, p, t, land)
        plant_ν = FT(1.44e-4)
        soil_params = land.soil.parameters
        θ_r = soil_params.θ_r
        ν = soil_params.ν
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
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
        evaluate!(Y.canopy.energy.T, atmos.T, t)

        Y.snow.S .= 0.0
        Y.snow.U .= 0.0
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
