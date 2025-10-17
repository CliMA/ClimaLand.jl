import StaticArrays: SVector
using ClimaLand.Canopy.PlantHydraulics

function canopy_integrator(;
    t0 = 0.0,
    tf = 86400.0,
    dt = 450.0,
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
    depth = FT(3.5)
    domain = ClimaLand.Domains.SphericalSurface(;
        radius = radius,
        nelements = n_horizontal_elements,
        npolynomial = 0,
    )
    surface_space = domain.space.surface

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
    ground = PrescribedGroundConditions{FT}()
    forcing = (; atmos, radiation, ground)
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    canopy = CanopyModel{FT}(domain, forcing, LAI, toml_dict)
    sim = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        canopy;
        diagnostics = (),
        updateat = update_drivers ? Second(dt) : Second(2 * (tf - t0)), # we still want to update drivers on init
        user_callbacks = (),
        set_ic! = Returns(nothing),
    )
    return sim

    # # Provide initial conditions for the canopy hydraulics model

    # ψ_stem_0 = FT(-1e5 / 9800)
    # ψ_leaf_0 = FT(-2e5 / 9800)

    # S_l_ini =
    #     inverse_water_retention_curve.(
    #         retention_model,
    #         [ψ_stem_0, ψ_leaf_0],
    #         plant_ν,
    #         plant_S_s,
    #     )

    # for i in 1:2
    #     Y.canopy.hydraulics.ϑ_l.:($i) .= augmented_liquid_fraction.(plant_ν, S_l_ini[i])
    # end;

end
