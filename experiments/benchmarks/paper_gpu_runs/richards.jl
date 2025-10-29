function richards_integrator(;
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

    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)

    # Below, the preprocess_func argument is used to
    # 1. Convert precipitation to be negative (as it is downwards)
    # 2. Convert mass flux to equivalent liquid water flux
    # Precipitation:
    precip = TimeVaryingInput(
        era5_ncdata_path,
        "mtpr",
        surface_space;
        start_date,
        regridder_type = :InterpolationsRegridder,
        file_reader_kwargs = (; preprocess_func = (data) -> FT(-data / 1000)),
    )
    forcing = (;
        atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip),
    )
    model = ClimaLand.Soil.RichardsModel{FT}(domain, forcing)
    function set_ic!(Y, p, t, model)
        domain = ClimaLand.get_domain(model)
        z = domain.fields.z
        lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
        function hydrostatic_profile(
            lat::FT,
            z::FT,
            ν::FT,
            θ_r::FT,
            α::FT,
            n::FT,
            S_s::FT,
            fmax,
        ) where {FT}
            m = 1 - 1 / n
            zmin = FT(-50.0)
            zmax = FT(0.0)

            z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5 * (fmax - 0.35) / 0.7)
            if z > z_∇
                S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
                ϑ_l = S * (ν - θ_r) + θ_r
            else
                ϑ_l = -S_s * (z - z_∇) + ν
            end
            return FT(ϑ_l)
        end

        # Set initial state values
        hydrology_cm = model.parameters.hydrology_cm
        ν = model.parameters.ν
        θ_r = model.parameters.S_s
        S_s = model.parameters.θ_r
        # Set initial state values
        vg_α = hydrology_cm.α
        vg_n = hydrology_cm.n
        vg_α = hydrology_cm.α
        vg_n = hydrology_cm.n
        f_max = ClimaLand.Soil.topmodel_fmax(domain.space.surface, FT)
        Y.soil.ϑ_l .=
            hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s, f_max)
    end
    # Create model update functions
    sim = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        model;
        diagnostics = (),
        updateat = update_drivers ? Second(dt) : Second(2 * (tf - t0)), # we still want to update drivers on init
        user_callbacks = (),
        set_ic!,
    )
    return sim
end
