using Test
import ClimaParams
using ClimaLand
using SurfaceFluxes
import ClimaLand.Parameters as LP
using ClimaLand.InlandWater
using Dates
using ClimaLand.Domains: Point
using ClimaLand: prescribed_forcing_era5
using ClimaUtilities.TimeManager: ITime
using ClimaCore
@testset "SlabLake model and make_* functions test" begin
    FT = Float32
    toml_dict = LP.create_toml_dict(FT)
    default_params_from_toml = SlabLakeParameters(toml_dict)
    longlat = FT.((50.6689, 41.9350))
    z_sfc = FT(0)
    domain = Point(; z_sfc, longlat)
    surface_space = domain.space.surface
    start_date = DateTime(2008)
    stop_date = start_date + Day(1)
    forcing = prescribed_forcing_era5(
        start_date,
        stop_date,
        domain.space.surface,
        toml_dict,
        FT;
        use_lowres_forcing = true,
    )
    lake = InlandWater.SlabLakeModel(FT, domain, forcing, toml_dict)
    @test ClimaLand.name(lake) == :lake
    Y, p, cds = ClimaLand.initialize(lake)
    @test propertynames(Y.lake) == prognostic_vars(lake)
    @test propertynames(p.lake) == auxiliary_vars(lake)
    @test InlandWater.get_turb_fluxes_type(FT, lake.boundary_conditions) ==
          NamedTuple{
        (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T),
        Tuple{FT, FT, FT, FT, FT},
    }
    @test ClimaLand.surface_roughness_model(lake, Y, p) ==
          SurfaceFluxes.ConstantRoughnessParams{FT}(
        default_params_from_toml.z_0m,
        default_params_from_toml.z_0b,
    )
    @test ClimaLand.surface_height(lake, Y, p) == FT(0)
    @test ClimaLand.surface_displacement_height(lake, Y, p) == FT(0)
    @test lake.inland_water_mask == ClimaCore.Fields.ones(surface_space)
    t0 = ITime(0, Second(1), start_date)
    ClimaLand.Simulations.set_lake_initial_conditions!(Y, p, t0, lake)
    tmp = ClimaCore.Fields.zeros(surface_space)
    evaluate!(tmp, forcing.atmos.T, t0)
    expected_U = @. InlandWater.lake_energy_from_temperature(
        tmp,
        default_params_from_toml,
    )
    @test Y.lake.U == expected_U
    set_initial_cache! = make_set_initial_cache(lake)
    set_initial_cache!(p, Y, t0)
    @test p.lake.T == tmp
    @test p.lake.q_l == ClimaCore.Fields.zeros(surface_space)
    @test all(parent(p.lake.albedo) .== default_params_from_toml.ice_albedo)
    @test all(
        parent(
            (
                p.drivers.P_liq .+ p.drivers.P_snow .+
                p.lake.turbulent_fluxes.vapor_flux
            ) .+ p.lake.runoff,
        ) .≈ 0,
    )
    compute_exp_tendency! = make_compute_exp_tendency(lake)
    dY = similar(Y)
    compute_exp_tendency!(dY, Y, p, t0)
    @test all(.~isnan.(parent(dY.lake.U)))
    @test all(.~isinf.(parent(dY.lake.U)))
end


@testset "Lake parameterization functions test" begin
    FT = Float32
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    params = InlandWater.SlabLakeParameters(toml_dict)

    # Warm lake (liquid)
    T_warm = FT(290)
    U_warm = InlandWater.lake_energy_from_temperature(T_warm, params)
    q_l_warm = InlandWater.lake_liquid_fraction(U_warm, params)
    T_warm_back = InlandWater.lake_temperature(U_warm, q_l_warm, params)
    @test q_l_warm ≈ FT(1)
    @test T_warm_back ≈ T_warm

    # Cold lake (frozen)
    T_cold = FT(260)
    U_cold = InlandWater.lake_energy_from_temperature(T_cold, params)
    q_l_cold = InlandWater.lake_liquid_fraction(U_cold, params)
    T_cold_back = InlandWater.lake_temperature(U_cold, q_l_cold, params)
    @test q_l_cold ≈ FT(0)
    @test T_cold_back ≈ T_cold

    # At freezing — T == T_freeze takes the ice branch, so q_l = 0
    T_freeze = FT(LP.T_freeze(earth_param_set))
    U_freeze = InlandWater.lake_energy_from_temperature(T_freeze, params)
    q_l_freeze = InlandWater.lake_liquid_fraction(U_freeze, params)
    @test q_l_freeze ≈ FT(0)

    @test InlandWater.lake_surface_albedo(FT(0.5), params) ≈
          FT(0.5) * params.liquid_albedo + FT(0.5) * params.ice_albedo

    # Zero runoff gives zero energy flux
    @test InlandWater.lake_runoff_energy_flux(
        FT(0),
        FT(290),
        FT(1),
        earth_param_set,
    ) ≈ FT(0)

    # Non-zero runoff gives non-zero energy flux
    E = InlandWater.lake_runoff_energy_flux(
        FT(1e-3),
        FT(290),
        FT(1),
        earth_param_set,
    )
    @test isfinite(E)
    @test E != FT(0)

    # When precipitation falls on a warm lake with no sediment flux,
    # the energy tendency should equal -(surface flux + runoff energy flux).
    # Here we verify the runoff energy budget: precipitation adds water at T_air,
    # runoff removes it at T_lake, and the net energy is accounted for.

    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _cp_l = LP.cp_l(earth_param_set)
    _T_ref = LP.T_0(earth_param_set)

    T_lake = FT(290)
    q_l = FT(1)  # fully liquid
    P_liq = FT(1e-4)  # precipitation rate (m/s)

    # Runoff = -F_water = -P_liq (no evaporation) to maintain constant depth
    runoff = -P_liq

    # Precipitation enthalpy flux into lake (at air temperature = lake temperature here)
    T_air = T_lake
    Q_precip = P_liq * _ρ_l * _cp_l * (T_air - _T_ref)

    # Runoff energy flux out of lake (at lake temperature)
    Q_runoff = InlandWater.lake_runoff_energy_flux(
        runoff,
        T_lake,
        q_l,
        earth_param_set,
    )

    # ρe_int at lake temperature for fully liquid
    ρe_lake = InlandWater.lake_volumetric_internal_energy(
        T_lake,
        q_l,
        earth_param_set,
    )

    # Runoff energy should be -P_liq * ρe(T_lake)
    @test Q_runoff ≈ -P_liq * ρe_lake

    # When T_air == T_lake and no radiation/turbulent fluxes, the energy tendency
    # contribution from precip + runoff should be: -Q_precip - Q_runoff
    # = -P_liq*ρ*cp*(T-T0) - (-P_liq)*ρe(T) = -P_liq*ρ*cp*(T-T0) + P_liq*ρe(T)
    # For fully liquid: ρe(T) = ρ*cp*(T-T0), so these cancel exactly.
    net_energy = -Q_precip - Q_runoff
    @test abs(net_energy) < eps(FT) * abs(Q_precip)

    # When evaporation removes water, runoff is positive (water supplied to maintain depth).
    # Verify the runoff energy flux is consistent

    T_lake = FT(295)
    q_l = FT(1)

    # Evaporation rate (negative water flux, so E_liq < 0 in sign convention)
    E_liq = FT(-5e-5)  # m/s evaporated
    F_water = E_liq  # no precipitation
    runoff = -F_water  # positive: water must be supplied

    Q_runoff = InlandWater.lake_runoff_energy_flux(
        runoff,
        T_lake,
        q_l,
        earth_param_set,
    )
    ρe_lake = InlandWater.lake_volumetric_internal_energy(
        T_lake,
        q_l,
        earth_param_set,
    )

    # Runoff is positive and carries energy into the lake at lake temperature
    @test Q_runoff ≈ runoff * ρe_lake
end
