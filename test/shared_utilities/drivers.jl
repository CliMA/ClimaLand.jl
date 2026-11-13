import ClimaComms
ClimaComms.@import_required_backends
import ClimaLand.Parameters as LP
using ClimaCore
using Test
using StaticArrays
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities
using ClimaLand
import Thermodynamics
import ClimaParams as CP
using Dates
import SurfaceFluxes
import SurfaceFluxes.Parameters as SFP

# A time-varying input defined only between two dates, like one read from files
struct BoundedTestInput{F} <:
       ClimaUtilities.TimeVaryingInputs.AbstractTimeVaryingInput
    data_handler::Tuple{DateTime, DateTime}
    f::F
end
Base.in(date::DateTime, input::BoundedTestInput) =
    input.data_handler[1] <= date <= input.data_handler[2]
function ClimaUtilities.TimeVaryingInputs.evaluate!(
    dest,
    input::BoundedTestInput,
    date,
)
    date in input || error("$date is outside the data")
    dest .= input.f(date)
    return nothing
end

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
    reference = SurfaceFluxes.ReferenceAboveSurface()
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
        reference,
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
    # use it as the floor of a FlooredDeardorffGustinessSpec
    ca = ClimaLand.CoupledAtmosphere{FT, FT}(FT(1), FT(1))
    @test ca.gustiness == FT(1)
    @test ca.reference === SurfaceFluxes.ReferenceAboveApparentSink()
    @test ClimaLand.gustiness_spec(ca) ==
          SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(1))
    @test SurfaceFluxes.without_floor(ClimaLand.gustiness_spec(ca)) ==
          SurfaceFluxes.FlooredDeardorffGustinessSpec(FT(0))
end

@testset "Turbulent flux selections, FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)
    z_0m = FT(0.01)
    z_0b = FT(0.001)
    roughness_model = SurfaceFluxes.ConstantRoughnessParams(z_0m, z_0b)
    gustiness = SurfaceFluxes.ConstantGustinessSpec(FT(1))
    reference = SurfaceFluxes.ReferenceAboveSurface()
    displ = FT(0)
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
        displ,
        (args...) -> FT(1),
        (args...) -> FT(1),
        gustiness,
        reference,
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
    @test selected.Δz_eff == FT(10) - displ
    @test selected.ustar > 0
    @test selected.ζ < 0 # unstable
    @test eltype(values(selected)) == FT
    # Above the apparent sink, the forcing lies z_0m above the displacement
    # height in addition
    above_sink = ClimaLand.turbulent_fluxes_at_a_point(
        Val((:Δz_eff,)),
        FT(101325),
        FT(290),
        FT(0.005),
        FT(2),
        FT(10),
        FT(295),
        FT(0.015),
        roughness_model,
        nothing,
        nothing,
        FT(0),
        displ,
        (args...) -> FT(1),
        (args...) -> FT(1),
        gustiness,
        SurfaceFluxes.ReferenceAboveApparentSink(),
        earth_param_set,
    )
    @test above_sink.Δz_eff ≈ FT(10) + z_0m
    # The deprecated form without a reference, which the ClimaCoupler bucket
    # calls, measures the forcing height from the ground
    no_reference(return_extra_fluxes) = ClimaLand.turbulent_fluxes_at_a_point(
        return_extra_fluxes,
        FT(101325),
        FT(290),
        FT(0.005),
        FT(2),
        FT(10),
        FT(295),
        FT(0.015),
        roughness_model,
        nothing,
        nothing,
        FT(0),
        displ,
        (args...) -> FT(1),
        (args...) -> FT(1),
        gustiness,
        earth_param_set,
    )
    @test no_reference(Val(false)) == fluxes(Val(false))
    @test no_reference(Val(true)) == fluxes(Val(true))
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
    # Arguments of different types, as with dual numbers, promote
    s_mixed = ClimaLand.screen_level_values(
        Float64(T_sfc),
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
    @test s_mixed.T isa Float64
    @test s_mixed.T ≈ s.T rtol = 10 * eps(FT)

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

@testset "CRU-JRA versions and cycled years" begin
    import ClimaUtilities.TimeVaryingInputs: evaluate!
    import ClimaUtilities.TimeManager: ITime
    @test ClimaLand.Artifacts.crujra_years() == 1901:2023
    @test ClimaLand.Artifacts.crujra_years("4.10") == 1901:2025
    @test_throws ErrorException ClimaLand.Artifacts.crujra_years("3.0")
    @test_throws ErrorException ClimaLand.Artifacts.crujra_forcing_data_folder_path(;
        version = "4.10",
    )

    crujra_path = mktempdir()
    for year in 2024:2025
        touch(joinpath(crujra_path, "crujra_forcing_data_$(year)_0.5x0.5.nc"))
    end
    find_paths(start_date, stop_date; version = "4.10") =
        ClimaLand.Artifacts.find_crujra_year_paths(
            start_date,
            stop_date;
            version,
            crujra_path,
        )
    @test basename.(find_paths(DateTime(2024), DateTime(2025, 6))) == [
        "crujra_forcing_data_2024_0.5x0.5.nc",
        "crujra_forcing_data_2025_0.5x0.5.nc",
    ]
    @test length(find_paths(DateTime(2025), DateTime(2026))) == 1
    @test_throws ErrorException find_paths(DateTime(2025), DateTime(2026, 2))
    @test_throws ErrorException find_paths(DateTime(2023), DateTime(2025))
    @test_throws ErrorException find_paths(
        DateTime(2024),
        DateTime(2025);
        version = "2.5",
    )

    years = 1901:1920
    @test ClimaLand.cycled_date(DateTime(1700, 3, 5, 6), years) ==
          DateTime(1920, 3, 5, 6)
    @test ClimaLand.cycled_date(DateTime(1901, 1, 1), years) == DateTime(1901)
    @test ClimaLand.cycled_date(DateTime(1921, 1, 1), years) == DateTime(1901)
    @test ClimaLand.cycled_date(DateTime(2025, 12, 31, 18), years) ==
          DateTime(1905, 12, 31, 18)
    @test ClimaLand.cycled_date(DateTime(1720, 2, 29), years) ==
          DateTime(1920, 2, 29)
    @test ClimaLand.cycled_date(DateTime(1712, 2, 29), years) ==
          DateTime(1912, 2, 29)
    @test ClimaLand.cycled_date(DateTime(1716, 2, 29), years) ==
          DateTime(1916, 2, 29)
    @test ClimaLand.cycled_date(DateTime(1704, 2, 29), years) ==
          DateTime(1904, 2, 29)
    @test ClimaLand.cycled_date(DateTime(1724, 2, 29), years) ==
          DateTime(1904, 2, 29)
    @test ClimaLand.cycled_date(DateTime(1800, 3, 1), years) ==
          DateTime(1920, 3, 1)
    @test ClimaLand.cycled_date(DateTime(1802, 2, 28), years) ==
          DateTime(1902, 2, 28)

    year_and_day = TimeVaryingInput(
        (date) -> Dates.year(date) + Dates.dayofyear(date) / 1000,
    )
    start_date = DateTime(1700)
    cycled = ClimaLand.CycledYearsInput(year_and_day, start_date, years)
    dest = zeros(1)
    seconds =
        Dates.value(Dates.Millisecond(DateTime(1745, 2, 1) - start_date)) / 1000
    evaluate!(dest, cycled, seconds)
    @test dest[1] ≈ 1905.032
    evaluate!(dest, cycled, ITime(seconds; epoch = start_date))
    @test dest[1] ≈ 1905.032

    # The hold-window interpolation of the state variables runs on the cycled
    # records: after the last cycled year, its samples come from the first one
    year_and_hour =
        TimeVaryingInput((date) -> Dates.year(date) + Dates.hour(date) / 100)
    held = ClimaLand.HoldWindowStateInput(
        ClimaLand.CycledYearsInput(year_and_hour, start_date, years),
        start_date,
        zeros(1);
        updateat = Dates.Hour(3),
        method = :linear,
    )
    sec_of(date) = Dates.value(Dates.Millisecond(date - start_date)) / 1000
    # Window 21:00-24:00 of 1720-12-31 (cycled 1920): midpoint 22:30, between
    # the samples of 1920-12-31 18:00 and of 1721-01-01 00:00 (cycled 1901)
    evaluate!(dest, held, sec_of(DateTime(1720, 12, 31, 21)))
    @test dest[1] ≈ 0.25 * 1920.18 + 0.75 * 1901
    @test ClimaLand.in_data_range(
        ClimaLand.CycledYearsInput(year_and_hour, start_date, years),
        DateTime(1721),
    )
end

@testset "Cos-zenith distribution of interval-mean shortwave" begin
    import ClimaUtilities.TimeVaryingInputs: evaluate!
    import ClimaUtilities.TimeManager: ITime
    import Insolation
    FT = Float64
    insol_params = LP.LandParameters(LP.create_toml_dict(FT)).insol_params
    # Equator, mid-latitudes, polar day (89.9°N in June), polar night (89.9°S)
    # with a nonzero and with a zero interval mean
    latitude = FT[0, 45, -30, 89.9, -89.9, -89.9]
    longitude = FT[0, -120, 150, 10, 10, 10]
    polar_day, polar_night, dark_zero = 4, 5, 6
    cosθs(date) = [
        Insolation.insolation(date, φ, λ, insol_params).μ for
        (φ, λ) in zip(latitude, longitude)
    ]
    # Interval means that differ between intervals, stamped at the interval start
    interval_mean(date) =
        (FT(100) + Dates.hour(floor(date, Dates.Hour(6)))) .*
        FT[1, 1.5, 0.5, 2, 1, 0]
    start_date = DateTime(2003, 6, 21)
    step = Dates.Minute(15)
    half_step = Dates.Millisecond(step) ÷ 2
    sw_input(updateat; start_date = start_date) =
        ClimaLand.CosZenithWeightedInput(
            TimeVaryingInput(interval_mean),
            start_date,
            latitude,
            longitude,
            insol_params;
            updateat,
        )
    @test_throws ErrorException ClimaLand.CosZenithWeightedInput(
        TimeVaryingInput(interval_mean),
        start_date,
        latitude,
        longitude,
        insol_params;
        quadrature_step = Dates.Minute(7),
    )
    @test_throws ErrorException sw_input(Dates.Hour(0))
    @test_throws ErrorException sw_input([0.0, 3600.0])
    dest = zeros(FT, length(latitude))
    ndark = 0
    # Mean over each 6-hour interval of the values held between updates
    for updateat in (Dates.Minute(15), Dates.Hour(1), Dates.Hour(3))
        sw = sw_input(updateat)
        nwindows = div(Dates.Hour(6), updateat)
        for interval_start in start_date .+ Dates.Hour(6) .* (0:7)
            window_starts = interval_start .+ updateat .* (0:(nwindows - 1))
            values = map(window_starts) do date
                evaluate!(dest, sw, date)
                copy(dest)
            end
            @test sum(values) / nwindows ≈ interval_mean(interval_start) rtol =
                10eps(FT)
            @test all(v -> all(>=(0), v), values)
            # Polar night: the zero mean stays zero; a nonzero mean is spread
            # evenly over the interval
            @test all(v[dark_zero] == 0 for v in values)
            @test all(
                v[polar_night] ≈ interval_mean(interval_start)[polar_night] for
                v in values
            )
            # Polar day: nearly constant sun, nearly constant shortwave
            @test all(
                isapprox(
                    v[polar_day],
                    interval_mean(interval_start)[polar_day];
                    rtol = 0.01,
                ) for v in values
            )
            # Sun down over the whole window, but up at some time in the
            # interval: no shortwave
            interval_cosθs =
                sum(o -> cosθs(interval_start + o), step .* (0:23) .+ half_step)
            offsets = step .* (0:(div(updateat, step) - 1)) .+ half_step
            for (date, value) in zip(window_starts, values)
                window_cosθs = sum(o -> cosθs(date + o), offsets)
                dark = (window_cosθs .== 0) .& (interval_cosθs .> 0)
                ndark += count(dark)
                @test all(value[dark] .== 0)
            end
        end
    end
    @test ndark > 0
    # At 15-minute updates, 45°N 120°W is dark before its sunrise in the
    # 12–18 UTC interval (≈ 04:25 local solar time), and the equator gets more
    # than the interval mean at noon
    sw = sw_input(Dates.Minute(15))
    evaluate!(dest, sw, DateTime(2003, 6, 21, 12))
    @test dest[2] == 0
    @test sw.window_cosθs[2] == 0
    @test dest[1] > interval_mean(DateTime(2003, 6, 21, 12))[1]
    # The held zenith angle is the shortwave-weighted mean over the window
    sw = sw_input(Dates.Hour(3))
    date = DateTime(2003, 6, 21, 6)
    evaluate!(dest, sw, date)
    nodes = date .+ step .* (0:11) .+ half_step
    μ = [cosθs(d)[1] for d in nodes]
    @test sw.window_cosθs[1] ≈ sum(μ .^ 2) / sum(μ)
    @test 0 < sw.window_cosθs[1] < 1
    # Windows that cross interval boundaries and midnight UTC: the mean over the
    # window is the mean of its parts
    sw2 = sw_input(Dates.Hour(2))
    sw1 = sw_input(Dates.Hour(1))
    window = DateTime(2003, 6, 21, 23)
    evaluate!(dest, sw2, window)
    crossing = copy(dest)
    evaluate!(dest, sw1, window)
    first_hour = copy(dest)
    evaluate!(dest, sw1, window + Dates.Hour(1))
    @test crossing ≈ (first_hour .+ dest) ./ 2 rtol = 10eps(FT)
    # and energy is conserved over whole days with misaligned 2-hour windows
    sw2 = sw_input(Dates.Hour(2))
    windows = DateTime(2003, 6, 21, 1) .+ Dates.Hour(2) .* (0:11)
    energy = sum(windows) do date
        evaluate!(dest, sw2, date)
        copy(dest)
    end
    intervals = DateTime(2003, 6, 21) .+ Dates.Hour(6) .* (0:4)
    sw1 = sw_input(Dates.Hour(1))
    expected = sum(DateTime(2003, 6, 21, 1) .+ Dates.Hour(1) .* (0:23)) do date
        evaluate!(dest, sw1, date)
        copy(dest)
    end
    @test energy .* 2 ≈ expected rtol = 10eps(FT)
    # Simulation times in seconds and as ITime
    sw = sw_input(Dates.Hour(3))
    date = DateTime(2003, 6, 22, 15)
    evaluate!(dest, sw, date)
    expected = copy(dest)
    seconds = Dates.value(Dates.Millisecond(date - start_date)) / 1000
    evaluate!(dest, sw, seconds)
    @test dest == expected
    evaluate!(dest, sw, ITime(seconds; epoch = start_date))
    @test dest == expected
    # The diffuse fraction is constant over each interval
    diffuse(date) = fill(FT(Dates.hour(date) + Dates.minute(date) / 60), 2)
    frac_diff =
        ClimaLand.IntervalStartInput(TimeVaryingInput(diffuse), start_date)
    fd = zeros(FT, 2)
    for date in start_date .+ Dates.Minute(15) .* (0:47)
        evaluate!(fd, frac_diff, date)
        @test all(fd .== Dates.hour(floor(date, Dates.Hour(6))))
    end
end

@testset "Cos-zenith shortwave in PrescribedRadiativeFluxes" begin
    import ClimaUtilities.TimeManager: ITime
    FT = Float32
    toml_dict = LP.create_toml_dict(FT)
    insol_params = LP.LandParameters(toml_dict).insol_params
    start_date = DateTime(2003, 3, 20)
    domain = ClimaLand.Domains.HybridBox(;
        xlim = FT.((0, 9)),
        ylim = FT.((0, 9)),
        zlim = FT.((1, 2)),
        nelements = (2, 2, 1),
        longlat = FT.((30, 10)),
    )
    coords = ClimaLand.Domains.coordinates(domain)
    sw = ClimaLand.CosZenithWeightedInput(
        TimeVaryingInput((t) -> FT(400)),
        start_date,
        coords.surface.lat,
        coords.surface.long,
        insol_params,
    )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        sw,
        TimeVaryingInput((t) -> FT(300)),
        start_date;
        cosθs = ClimaLand.WindowCosZenith(sw),
        toml_dict,
        frac_diff = ClimaLand.IntervalStartInput(
            TimeVaryingInput((t) -> FT(0.3)),
            start_date,
        ),
    )
    p = (; drivers = ClimaLand.initialize_drivers((radiation,), coords))
    update! = ClimaLand.make_update_drivers((radiation,))
    sw_sum = zero(FT)
    for k in 0:7
        update!(p, ITime(3 * 3600 * k; epoch = start_date))
        sw_sum += sum(ClimaCore.Fields.field2array(p.drivers.SW_d)) / 4
        @test ClimaCore.Fields.field2array(p.drivers.cosθs) ==
              ClimaCore.Fields.field2array(sw.window_cosθs)
        if k == 4 # 12 UTC, near local noon: sun up, cosθs between 0 and 1
            @test all(>(0), ClimaCore.Fields.field2array(p.drivers.cosθs))
        end
    end
    @test sw_sum / 8 ≈ 400 rtol = 1e-5
    @test all(==(FT(0.3)), ClimaCore.Fields.field2array(p.drivers.frac_diff))
end

@testset "Interval-mean fluxes held between driver updates" begin
    import ClimaUtilities.TimeVaryingInputs: evaluate!
    import ClimaUtilities.TimeManager: ITime
    FT = Float64
    start_date = DateTime(2003, 1, 1)
    # Interval means stamped at the interval start, different in every interval
    interval_mean(date) =
        FT[1, -2, 0] .* (1 + Dates.value(date - start_date) / 3_600_000)
    flux(updateat) = ClimaLand.IntervalMeanInput(
        TimeVaryingInput(interval_mean),
        start_date,
        zeros(FT, 3);
        updateat,
    )
    @test_throws ErrorException flux(Dates.Hour(0))
    dest = zeros(FT, 3)
    for updateat in (
            Dates.Minute(15),
            Dates.Hour(1),
            Dates.Hour(2),
            Dates.Hour(3),
            Dates.Hour(4),
            Dates.Hour(5),
            Dates.Hour(6),
            Dates.Hour(8),
        ),
        offset in (Dates.Hour(0), Dates.Hour(1), Dates.Minute(15))

        input = flux(updateat)
        h = Dates.value(Dates.Millisecond(updateat)) / 3_600_000
        # Windows that tile more than a day, against the integral of the
        # piecewise-constant flux over the same span (minute by minute)
        windows =
            (start_date + offset) .+ updateat .* (0:div(Dates.Day(1), updateat))
        total = sum(windows) do date
            evaluate!(dest, input, date)
            dest .* h
        end
        minutes =
            windows[1]:Dates.Minute(1):(windows[end] + updateat - Dates.Minute(
                1,
            ))
        expected =
            sum(m -> interval_mean(floor(m, Dates.Hour(6))) ./ 60, minutes)
        @test total ≈ expected rtol = 1e-12
    end
    # A window inside an interval returns its mean; a window across a boundary
    # is weighted by the overlaps
    input = flux(Dates.Hour(2))
    evaluate!(dest, input, start_date + Dates.Hour(3))
    @test dest == interval_mean(start_date)
    evaluate!(dest, input, start_date + Dates.Minute(270))
    @test dest ≈
          0.75 .* interval_mean(start_date) .+
          0.25 .* interval_mean(start_date + Dates.Hour(6))
    # Seconds, ITime and DateTime times agree
    input = flux(Dates.Hour(3))
    evaluate!(dest, input, 5 * 3600.0)
    seconds = copy(dest)
    evaluate!(dest, input, ITime(5 * 3600; epoch = start_date))
    @test dest == seconds
    evaluate!(dest, input, start_date + Dates.Hour(5))
    @test dest == seconds
end

@testset "Sampled states at the middle of the driver hold window" begin
    import ClimaUtilities.TimeVaryingInputs: evaluate!
    import ClimaUtilities.TimeManager: ITime
    FT = Float64
    start_date = DateTime(2003, 1, 1)
    hours(date) = Dates.value(date - start_date) / 3_600_000
    state(f; kwargs...) = ClimaLand.HoldWindowStateInput(
        TimeVaryingInput(f),
        start_date,
        zeros(FT, 2);
        kwargs...,
    )
    @test_throws ErrorException state(identity; method = :spline)
    @test_throws ErrorException state(identity; updateat = Dates.Hour(0))
    dest = zeros(FT, 2)
    # The cubic is exact for quadratics, the linear interpolation for lines,
    # both evaluated at the middle of the window
    quadratic(date) =
        FT[1, -1] .* (280 + 0.5 * hours(date) - 0.01 * hours(date)^2)
    line(date) = FT[1, 2] .* (3 - 0.1 * hours(date))
    for updateat in (Dates.Minute(15), Dates.Hour(1), Dates.Hour(3)),
        date in start_date .+ Dates.Minute(45) .* (0:40)

        midpoint = date + Dates.Millisecond(updateat) ÷ 2
        evaluate!(dest, state(quadratic; updateat), date)
        @test dest ≈ quadratic(midpoint) rtol = 1e-12
        evaluate!(dest, state(line; updateat, method = :linear), date)
        @test dest ≈ line(midpoint) rtol = 1e-12
    end
    # The harmonic interpolation is exact for a mean, a diurnal and a
    # semidiurnal cosine
    harmonics(date) =
        FT[1, 2] .* (
            3 +
            2 * cospi(hours(date) / 12) +
            sinpi(hours(date) / 12) +
            0.5 * cospi(hours(date) / 6)
        )
    for updateat in (Dates.Hour(1), Dates.Hour(3)),
        date in start_date .+ Dates.Minute(45) .* (0:40)

        evaluate!(dest, state(harmonics; updateat, method = :harmonic), date)
        @test dest ≈ harmonics(date + Dates.Millisecond(updateat) ÷ 2) rtol =
            1e-12
    end
    # All pass through the samples
    sample(date) = FT[sin(hours(date)), cos(hours(date))]
    for method in (:cubic, :harmonic, :linear)
        input = state(sample; updateat = Dates.Hour(2), method)
        evaluate!(dest, input, start_date + Dates.Hour(5))
        @test dest ≈ sample(start_date + Dates.Hour(6))
    end
    # A diurnal cycle sampled every 6 hours: the cubic follows it more closely
    # and keeps more of its range than the linear interpolation
    diurnal(date) = FT[1, 1] .* cospi((hours(date) - 15) / 12)
    errors_and_ranges = map((:cubic, :linear)) do method
        input = state(diurnal; updateat = Dates.Hour(1), method)
        dates = start_date .+ Dates.Hour.(0:23)
        values = map(dates) do date
            evaluate!(dest, input, date)
            dest[1]
        end
        truth = [diurnal(date + Dates.Minute(30))[1] for date in dates]
        (sum(abs, values .- truth) / 24, maximum(values) - minimum(values))
    end
    @test errors_and_ranges[1][1] < errors_and_ranges[2][1] / 2
    @test errors_and_ranges[1][2] > errors_and_ranges[2][2]
    # Cubic overshoots below zero are removed if requested
    spike(date) = fill(FT(Dates.hour(date) in (0, 18) ? 10 : 0), 2)
    for nonnegative in (false, true)
        input = state(spike; updateat = Dates.Hour(1), nonnegative)
        evaluate!(dest, input, start_date + Dates.Hour(8))
        @test nonnegative ? all(dest .== 0) : all(dest .< 0)
    end
    # Samples outside the data fall back on the nearest sample of the interval
    data_end = start_date + Dates.Hour(12)
    bounded = BoundedTestInput((start_date, data_end), quadratic)
    input = ClimaLand.HoldWindowStateInput(
        bounded,
        start_date,
        zeros(FT, 2);
        updateat = Dates.Hour(1),
    )
    evaluate!(dest, input, start_date)
    @test all(isfinite, dest)
    evaluate!(dest, input, data_end - Dates.Hour(1))
    @test all(isfinite, dest)
    @test minimum(quadratic(data_end)) <= minimum(dest)
    # Seconds, ITime and DateTime times agree
    input = state(quadratic; updateat = Dates.Hour(3))
    evaluate!(dest, input, 7 * 3600.0)
    seconds = copy(dest)
    evaluate!(dest, input, ITime(7 * 3600; epoch = start_date))
    @test dest == seconds
    evaluate!(dest, input, start_date + Dates.Hour(7))
    @test dest == seconds
end

@testset "Interval-mean and hold-window inputs in PrescribedAtmosphere" begin
    import ClimaUtilities.TimeManager: ITime
    FT = Float32
    toml_dict = LP.create_toml_dict(FT)
    start_date = DateTime(2003, 3, 20)
    domain = ClimaLand.Domains.HybridBox(;
        xlim = FT.((0, 9)),
        ylim = FT.((0, 9)),
        zlim = FT.((1, 2)),
        nelements = (2, 2, 1),
        longlat = FT.((30, 10)),
    )
    coords = ClimaLand.Domains.coordinates(domain)
    scratch = zero(coords.surface.lat)
    updateat = Dates.Hour(3)
    hours(date) = Dates.value(date - start_date) / 3_600_000
    precip_mean(date) = FT(-1e-7 * (1 + Dates.hour(date) / 6))
    temperature(date) = FT(280 + 0.5 * hours(date) - 0.01 * hours(date)^2)
    atmos = ClimaLand.PrescribedAtmosphere(
        ClimaLand.IntervalMeanInput(
            TimeVaryingInput(precip_mean),
            start_date,
            scratch;
            updateat,
        ),
        TimeVaryingInput((t) -> FT(0)),
        ClimaLand.HoldWindowStateInput(
            TimeVaryingInput(temperature),
            start_date,
            scratch;
            updateat,
        ),
        TimeVaryingInput((t) -> FT(2)),
        TimeVaryingInput((t) -> FT(0.005)),
        TimeVaryingInput((t) -> FT(1e5)),
        start_date,
        FT(10),
        toml_dict,
    )
    p = (; drivers = ClimaLand.initialize_drivers((atmos,), coords))
    update! = ClimaLand.make_update_drivers((atmos,))
    precip_sum = zero(FT)
    for k in 0:7
        date = start_date + updateat * k
        update!(p, ITime(3 * 3600 * k; epoch = start_date))
        precip_sum += sum(ClimaCore.Fields.field2array(p.drivers.P_liq)) / 4
        @test all(
            ClimaCore.Fields.field2array(p.drivers.T) .≈
            temperature(date + Dates.Minute(90)),
        )
    end
    expected = sum(precip_mean.(start_date .+ Dates.Hour(6) .* (0:3))) * 2
    @test precip_sum ≈ expected rtol = 1e-5
end
