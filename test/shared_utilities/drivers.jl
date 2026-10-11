import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand.Parameters as LP
using ClimaCore
using Test
using StaticArrays
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
import Thermodynamics
import ClimaParams as CP
using Dates
import SurfaceFluxes
import SurfaceFluxes.Parameters as SFP

FT = Float32
@testset "Default model, FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    (pa, pr) = ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)

    domain = ClimaLand.Domains.Plane(;
        xlim = FT.((1.0, 2.0)),
        ylim = FT.((1.0, 2.0)),
        nelements = (1, 1),
        periodic = (true, true),
    )
    coords = ClimaLand.Domains.coordinates(domain)
    zero_instance = ClimaCore.Fields.zeros(axes(coords.surface))
    @test ClimaLand.initialize_drivers((pp,), coords) ==
          NamedTuple{(:P_liq,)}((zero_instance,))
    @test ClimaLand.initialize_drivers((), coords) == (;)
    pa_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    pa_vals = ([zero_instance for k in pa_keys]...,)
    all_pa_keys = (pa_keys...,)
    all_pa_vals = (pa_vals...,)
    @test ClimaLand.initialize_drivers((pa,), coords) ==
          NamedTuple{all_pa_keys}(all_pa_vals)
    pr_keys = (:SW_d, :LW_d, :cosθs, :frac_diff)
    pr_vals = ([zero_instance for k in pr_keys]...,)
    all_papr_keys = (pa_keys..., pr_keys...)
    all_papr_vals = (pa_vals..., pr_vals...)
    @test ClimaLand.initialize_drivers((pa, pr), coords) ==
          NamedTuple{all_papr_keys}(all_papr_vals)
end

@testset "Driver update functions" begin
    toml_dict = LP.create_toml_dict(FT)
    f = TimeVaryingInput((t) -> 10.0)
    pa = ClimaLand.PrescribedAtmosphere(f, f, f, f, f, f, f, FT(1), toml_dict)
    pr = ClimaLand.PrescribedRadiativeFluxes(FT, f, f, f)
    domain = ClimaLand.Domains.HybridBox(;
        xlim = FT.((1.0, 2.0)),
        ylim = FT.((1.0, 2.0)),
        zlim = FT.((1.0, 2.0)),
        nelements = (1, 1, 1),
    )
    coords = ClimaLand.Domains.coordinates(domain)
    sfc_instance = ClimaCore.Fields.zeros(axes(coords.surface)) .+ 10
    p = (; drivers = ClimaLand.initialize_drivers((), coords))
    nothing_update! = ClimaLand.make_update_drivers(())
    nothing_update!(p, 0.0)
    @test p.drivers == (;)
    p = (; drivers = ClimaLand.initialize_drivers((pa,), coords))
    atmos_only_update! = ClimaLand.make_update_drivers((pa,))
    atmos_only_update!(p, 0.0)
    @test p.drivers.P_liq == sfc_instance
    @test p.drivers.P_snow == sfc_instance
    @test p.drivers.P == sfc_instance
    @test p.drivers.T == sfc_instance
    @test p.drivers.q == sfc_instance
    @test p.drivers.u == sfc_instance
    @test p.drivers.c_co2 == (sfc_instance .* 0 .+ FT(4.2e-4))

    p = (; drivers = ClimaLand.initialize_drivers((pa, pr), coords))
    update! = ClimaLand.make_update_drivers((pa, pr))
    update!(p, 0.0)
    @test p.drivers.P_liq == sfc_instance
    @test p.drivers.P_snow == sfc_instance
    @test p.drivers.P == sfc_instance
    @test p.drivers.T == sfc_instance
    @test p.drivers.q == sfc_instance
    @test p.drivers.u == sfc_instance
    @test p.drivers.c_co2 == (sfc_instance .* 0 .+ FT(4.2e-4))
    @test p.drivers.SW_d == sfc_instance
    @test p.drivers.LW_d == sfc_instance
    @test all(isnan.(parent(p.drivers.cosθs)))
    @test all(isnan.(parent(p.drivers.frac_diff)))

    p = (; drivers = ClimaLand.initialize_drivers((pr,), coords))
    rad_only_update! = ClimaLand.make_update_drivers((pr,))
    rad_only_update!(p, 0.0)
    @test p.drivers.SW_d == sfc_instance
    @test p.drivers.LW_d == sfc_instance
    @test all(isnan.(parent(p.drivers.cosθs)))
    @test all(isnan.(parent(p.drivers.frac_diff)))

    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)
    precip_update! = ClimaLand.make_update_drivers((pp,))
    p = (; drivers = ClimaLand.initialize_drivers((pp,), coords))
    precip_update!(p, 0.0)
    @test p.drivers.P_liq == sfc_instance .* 0 .- FT(1)
end

@testset "PrescribedAtmosphere and PrescribedRadiativeFluxes show and summary" begin
    toml_dict = LP.create_toml_dict(FT)
    f = TimeVaryingInput((t) -> 10.0)
    start_date = DateTime(2005, 1, 1)
    for x in (
        ClimaLand.PrescribedAtmosphere(
            f,
            f,
            f,
            f,
            f,
            f,
            start_date,
            FT(1),
            toml_dict,
        ),
        ClimaLand.PrescribedRadiativeFluxes(FT, f, f, start_date),
    )
        typename = string(nameof(typeof(x)))

        out = sprint(show, MIME("text/plain"), x)
        @test occursin(typename, out)
        @test count(==('\n'), out) <= 10

        out2 = sprint(show, x)
        @test occursin(typename, out2)
        @test !occursin('\n', out2)
        out3 = sprint(show, MIME("text/plain"), x; context = :compact => true)
        @test out2 == out3

        out_summary = sprint(summary, x)
        @test occursin(typename, out_summary)
        @test !occursin('\n', out_summary)
    end
end

@testset "Gustiness models, FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    sf_params = LP.surface_fluxes_parameters(earth_param_set)
    β = SFP.gustiness_coeff(sf_params)
    z_i = SFP.gustiness_zi(sf_params)

    spec = SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1))
    @test spec isa SurfaceFluxes.AbstractGustinessSpec
    # The gustiness is a function of the stability parameter and the state,
    # so SurfaceFluxes uses its closed-form friction velocity
    @test !SurfaceFluxes.depends_on_ustar(spec)
    # Stable or neutral (non-positive buoyancy flux): only the floor applies
    for B in (FT(0), FT(-0.01))
        @test SurfaceFluxes.gustiness_value(spec, sf_params, B) == FT(1)
    end
    # Unstable: Deardorff value when it exceeds the floor
    B = FT(0.02)
    expected = β * cbrt(B * z_i)
    @test expected > 1
    @test SurfaceFluxes.gustiness_value(spec, sf_params, B) ≈ expected
    # A floor above the convective value wins
    @test SurfaceFluxes.gustiness_value(
        SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(10)),
        sf_params,
        B,
    ) == FT(10)

    # Floors
    c = SurfaceFluxes.ConstantGustinessSpec(FT(2))
    @test SurfaceFluxes.minimum_wind_speed(spec, sf_params) == FT(1)
    @test SurfaceFluxes.minimum_wind_speed(c, sf_params) == FT(2)
    @test SurfaceFluxes.minimum_wind_speed(
        SurfaceFluxes.DeardorffGustinessSpec(),
        sf_params,
    ) == 0

    # A flux solve with the model: in calm unstable conditions the effective
    # wind speed is at least the floor, so the fluxes exceed those of a solve
    # without gustiness; the solve is a function of its inputs alone
    roughness_model = SurfaceFluxes.ConstantRoughnessParams(FT(0.01), FT(0.001))
    solve(gustiness) = ClimaLand.surface_fluxes_at_a_point(
        FT(300), # T_sfc
        FT(0.015), # q_sfc
        nothing,
        nothing,
        FT(101325), # P_atmos
        FT(290), # T_atmos
        FT(0.005), # q_atmos
        FT(0.1), # u_atmos
        FT(10), # h_atmos
        FT(0), # h_sfc
        FT(0), # displ
        roughness_model,
        gustiness,
        earth_param_set,
    )
    with_gust = solve(spec)
    @test with_gust == solve(spec)
    no_gust = solve(SurfaceFluxes.ConstantGustinessSpec(FT(0)))
    floor_only = solve(SurfaceFluxes.ConstantGustinessSpec(FT(1)))
    @test with_gust.shf > no_gust.shf > 0 # upward sensible heat flux
    @test with_gust.shf >= floor_only.shf
    @test with_gust.ustar >= floor_only.ustar
    @test with_gust.ζ < 0 # unstable
    @test isfinite(with_gust.lhf) && isfinite(with_gust.L_MO)
    # Self-consistency: the effective wind speed of the solve, u* / sqrt(Cd),
    # is the gustiness β w* of its own buoyancy flux (the mean wind is
    # negligible here), and a solve with that constant gustiness agrees to
    # within the tolerance of the stability solve
    ΔU = with_gust.ustar / sqrt(with_gust.Cd)
    B = -with_gust.ustar^3 / (SFP.von_karman_const(sf_params) * with_gust.L_MO)
    @test B > 0
    @test ΔU ≈ β * cbrt(B * z_i) rtol = sqrt(eps(FT))
    @test ΔU > 1
    same = solve(SurfaceFluxes.ConstantGustinessSpec(ΔU))
    @test same.shf ≈ with_gust.shf rtol = 1e-3
    @test same.ustar ≈ with_gust.ustar rtol = 1e-3

    # Atmospheric drivers accept a number or a gustiness model
    f = TimeVaryingInput((t) -> 10.0)
    pa = ClimaLand.PrescribedAtmosphere(f, f, f, f, f, f, f, FT(1), toml_dict)
    @test pa.gustiness == SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1))
    pa2 = ClimaLand.PrescribedAtmosphere(
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        FT(1),
        toml_dict;
        gustiness = 2,
    )
    @test pa2.gustiness == SurfaceFluxes.ConstantGustinessSpec(FT(2))
    @test ClimaLand.gustiness_spec(pa2) === pa2.gustiness
    # A CoupledAtmosphere keeps the number the coupler reads; the flux solves
    # use it as a constant gustiness
    ca = ClimaLand.CoupledAtmosphere{FT, FT}(FT(1), FT(1))
    @test ca.gustiness == FT(1)
    @test ClimaLand.gustiness_spec(ca) ==
          SurfaceFluxes.ConstantGustinessSpec(FT(1))
end

@testset "Turbulent flux selections, FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    roughness_model = SurfaceFluxes.ConstantRoughnessParams(FT(0.01), FT(0.001))
    gustiness = SurfaceFluxes.ConstantGustinessSpec(FT(1))
    fluxes(stored) = ClimaLand.turbulent_fluxes_at_a_point(
        stored,
        FT(101325), # P_atmos
        FT(290), # T_atmos
        FT(0.005), # q_atmos
        FT(2), # u_atmos
        FT(10), # h_atmos
        FT(295), # T_sfc
        FT(0.015), # q_sfc
        roughness_model,
        nothing,
        nothing,
        FT(0), # h_sfc
        FT(0), # displ
        (args...) -> FT(1),
        (args...) -> FT(1),
        gustiness,
        earth_param_set,
    )
    # The Boolean selections are the contract of ClimaCoupler, which
    # evaluates them directly into its flux fields
    @test propertynames(fluxes(Val(false))) ==
          (:lhf, :shf, :vapor_flux, :∂lhf∂T, :∂shf∂T)
    @test propertynames(fluxes(Val(true))) == (
        :lhf,
        :shf,
        :vapor_flux,
        :∂lhf∂T,
        :∂shf∂T,
        :ρτxz,
        :ρτyz,
        :buoyancy_flux,
    )
    # A tuple of names selects from the full output, in the order given
    stored = Val((:ustar, :shf, :T_sfc, :q_sfc, :ζ, :Δz_eff, :buoyancy_flux))
    selected = fluxes(stored)
    @test propertynames(selected) ==
          (:ustar, :shf, :T_sfc, :q_sfc, :ζ, :Δz_eff, :buoyancy_flux)
    @test selected.shf == fluxes(Val(false)).shf
    @test selected.buoyancy_flux == fluxes(Val(true)).buoyancy_flux
    @test selected.T_sfc == FT(295)
    @test selected.q_sfc == FT(0.015)
    @test selected.Δz_eff == FT(10)
    @test selected.ustar > 0
    @test selected.ζ < 0 # unstable
    @test eltype(values(selected)) == FT
end

@testset "CoupledAtmosphere and CoupledRadiativeFluxes initialization" begin
    domain = ClimaLand.Domains.global_domain(FT)
    coords = ClimaLand.Domains.coordinates(domain)

    atmos = ClimaLand.CoupledAtmosphere{FT, FT}(FT(1), FT(1))
    radiation = ClimaLand.CoupledRadiativeFluxes{FT}()
    p = (; drivers = ClimaLand.initialize_drivers((atmos, radiation), coords))

    @test keys(p.drivers) == (
        :P_liq,
        :P_snow,
        :c_co2,
        :T,
        :P,
        :q,
        :u,
        :SW_d,
        :LW_d,
        :cosθs,
        :frac_diff,
    )
end

@testset "CoupledRadiativeFluxes" begin
    toml_dict = LP.create_toml_dict(FT)
    start_date = DateTime(Date(2020, 6, 15), Time(12, 0, 0))
    domain = ClimaLand.Domains.HybridBox(;
        xlim = FT.((0, 9)),
        ylim = FT.((0, 9)),
        zlim = FT.((1, 2)),
        nelements = (10, 10, 1),
        longlat = FT.((0, 0)),
    )
    coords = ClimaLand.Domains.coordinates(domain)

    # test CoupledRadiativeFluxes with no start_date provided (will not update cosθs)
    crf_no_zenith = ClimaLand.CoupledRadiativeFluxes{FT}()
    p = (; drivers = ClimaLand.initialize_drivers((crf_no_zenith,), coords))
    p.drivers.cosθs .= FT(0)
    no_update = ClimaLand.make_update_drivers((crf_no_zenith,))
    no_update(p, 0)
    @test all(isequal(FT(0)), ClimaCore.Fields.field2array(p.drivers.cosθs))
    crf = ClimaLand.CoupledRadiativeFluxes{FT}(
        start_date;
        latitude = coords.surface.lat,
        longitude = coords.surface.long,
        toml_dict,
    )
    p = (; drivers = ClimaLand.initialize_drivers((crf,), coords))
    update_cosθs_only = ClimaLand.make_update_drivers((crf,))
    update_cosθs_only(p, 0) # populate cosθs with cos(zenith) at noon mid-summer at equator
    @test all(
        x -> isapprox(x, 0.95; atol = 0.05),
        ClimaCore.Fields.field2array(p.drivers.cosθs),
    )
    update_cosθs_only(p, 60 * 60 * 12) # populate with cos(zenith) at night
    @test all((==)(0), ClimaCore.Fields.field2array(p.drivers.cosθs)) # cos(zenith angle) at nighttime should be 0
end

@testset "Ground Conditions" begin
    for FT in (Float32, Float64)
        soil_driver = PrescribedGroundConditions{FT}()
        prognostic_soil_driver = ClimaLand.PrognosticGroundConditions{FT}()
        @test ClimaLand.Canopy.ground_albedo_PAR(
            Val((:canopy,)),
            soil_driver,
            nothing,
            nothing,
            nothing,
        ) == FT(0.2)
        @test ClimaLand.Canopy.ground_albedo_NIR(
            Val((:canopy,)),
            soil_driver,
            nothing,
            nothing,
            nothing,
        ) == FT(0.4)
        dest = [-1.0]
        t = 2.0

        evaluate!(dest, soil_driver.θ, t)
        @test dest[1] == FT(0.4)
        evaluate!(dest, soil_driver.T, t)
        @test dest[1] == FT(298.0)

        domain = ClimaLand.Domains.Plane(;
            xlim = FT.((1.0, 2.0)),
            ylim = FT.((1.0, 2.0)),
            nelements = (1, 1),
            periodic = (true, true),
        )
        coords = ClimaLand.Domains.coordinates(domain)
        zero_instance = ClimaCore.Fields.zeros(axes(coords.surface))
        p_soil_driver = (;
            drivers = (;
                θ = copy(zero_instance),
                T_ground = copy(zero_instance),
            )
        )
        @test ClimaLand.initialize_drivers((soil_driver,), coords) ==
              p_soil_driver.drivers
        update_drivers! = make_update_drivers((soil_driver,))
        update_drivers!(p_soil_driver, 0.0)
        @test p_soil_driver.drivers.θ == zero_instance .+ FT(0.4)
        @test p_soil_driver.drivers.T_ground == zero_instance .+ 298

        @test ClimaLand.initialize_drivers((prognostic_soil_driver,), coords) ==
              (;)
        update_drivers! = make_update_drivers((prognostic_soil_driver,))
        p_soil_driver.drivers.θ .= FT(0.1)
        p_soil_driver.drivers.T_ground .= -1
        update_drivers!(p_soil_driver, 0.0)
        # no change
        @test p_soil_driver.drivers.θ == (zero_instance .+ FT(0.1))
        @test p_soil_driver.drivers.T_ground == (zero_instance .- 1)
    end
end

@testset "Screen-level reconstruction, FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    sf_params = LP.surface_fluxes_parameters(earth_param_set)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    κ = SFP.von_karman_const(sf_params)
    g = LP.grav(earth_param_set)
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    # The neutral heat profile carries the neutral turbulent Prandtl number
    Pr_0 = SFP.Pr_0(sf_params)

    z0m, z0h, Δz_eff = FT(0.1), FT(0.01), FT(20)
    T_sfc, T_air, q_sfc, q_air, ustar =
        FT(300), FT(290), FT(0.02), FT(0.01), FT(0.3)
    s = ClimaLand.screen_level_values(
        T_sfc,
        q_sfc,
        ustar,
        FT(0),
        Δz_eff,
        z0m,
        z0h,
        T_air,
        q_air,
        FT(2),
        FT(10),
        earth_param_set,
    )
    r = log((z0h + 2) / z0h) / log(Δz_eff / z0h)
    @test s.T ≈
          T_sfc + (T_air - T_sfc) * r + g / cp_d * (r * Δz_eff - (z0h + 2))
    @test s.q ≈ q_sfc + (q_air - q_sfc) * r
    @test s.u ≈ ustar / κ * log((z0m + 10) / z0m)
    @test s.g_h ≈ κ * ustar / (Pr_0 * log(Δz_eff / z0h))
    # Forcing at or below the screen height: forcing values are returned
    s_low = ClimaLand.screen_level_values(
        T_sfc,
        q_sfc,
        ustar,
        FT(0),
        FT(1.5),
        z0m,
        z0h,
        T_air,
        q_air,
        FT(2),
        FT(10),
        earth_param_set,
    )
    @test s_low.T ≈ T_air
    @test s_low.q ≈ q_air
    @test s_low.u ≈ ustar / κ * log(FT(1.5) / z0m)
    # Screen height at or below the roughness length is clamped to z0h
    s_z0 = ClimaLand.screen_level_values(
        T_sfc,
        q_sfc,
        ustar,
        FT(0),
        Δz_eff,
        z0m,
        z0h,
        T_air,
        q_air,
        FT(-1),
        FT(-1),
        earth_param_set,
    )
    @test s_z0.T ≈ T_sfc - g / cp_d * z0h
    @test s_z0.q ≈ q_sfc
    @test s_z0.u == 0

    # Weighted mean over surfaces: area fraction times conductance
    s1 = (; T = FT(1), q = FT(1), u = FT(1), g_h = FT(2))
    s2 = (; T = FT(3), q = FT(3), u = FT(3), g_h = FT(1))
    @test ClimaLand.screen_level_mean(Val(:T), (FT(1), s1)) == 1
    @test ClimaLand.screen_level_mean(Val(:T), (FT(1), s1), (FT(1), s2)) ≈
          (2 * 1 + 1 * 3) / 3
    @test ClimaLand.screen_level_mean(Val(:u), (FT(0.5), s1), (FT(1), s2)) ≈
          (1 * 1 + 1 * 3) / 2
    s0 = (; T = FT(7), q = FT(0), u = FT(0), g_h = FT(0))
    @test ClimaLand.screen_level_mean(Val(:T), (FT(1), s0), (FT(1), s0)) == 7
end
