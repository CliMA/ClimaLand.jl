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
    zero_thermal_state = ClimaCore.Fields.zeros(
        Thermodynamics.PhaseEquil{FT},
        axes(coords.surface),
    )
    pa_vals = ([zero_instance for k in pa_keys]...,)
    all_pa_keys = (pa_keys..., :thermal_state)
    all_pa_vals = (pa_vals..., zero_thermal_state)
    @test ClimaLand.initialize_drivers((pa,), coords) ==
          NamedTuple{all_pa_keys}(all_pa_vals)
    pr_keys = (:SW_d, :LW_d, :cosθs, :frac_diff)
    pr_vals = ([zero_instance for k in pr_keys]...,)
    all_papr_keys = (pa_keys..., :thermal_state, pr_keys...)
    all_papr_vals = (pa_vals..., zero_thermal_state, pr_vals...)
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


    soc = TimeVaryingInput((t) -> 1.0)
    subsfc_instance = ClimaCore.Fields.zeros(axes(coords.subsurface)) .+ 1
    soc_driver = ClimaLand.PrescribedSoilOrganicCarbon{FT}(soc)
    soc_update! = ClimaLand.make_update_drivers((soc_driver,))
    soc_p = (; drivers = ClimaLand.initialize_drivers((soc_driver,), coords))
    soc_update!(soc_p, 0.0)
    @test soc_p.drivers.soc == subsfc_instance
end

@testset "Dewpoint to RH" begin
    Td = 288.0:0.5:293.0
    Ta = 290.0:1.0:300.0
    rh = ClimaLand.rh_from_dewpoint.(Td, Ta)
    soln = [
        0.9701877799488919,
        0.9629770622625089,
        0.9558600702031558,
        0.9488352550428512,
        0.9419010988316328,
        0.9350561136897647,
        0.928298841118352,
        0.9216278513277926,
        0.9150417425835539,
        0.9085391405688034,
        0.9021186977633818,
    ]
    @test all(rh .- soln .≈ 0)
end

@testset "CoupledAtmosphere and CoupledRadiativeFluxes initialization" begin
    domain = ClimaLand.Domains.global_domain(FT)
    coords = ClimaLand.Domains.coordinates(domain)

    atmos = ClimaLand.CoupledAtmosphere{FT}()
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
        :thermal_state,
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
    @test all((<)(0), ClimaCore.Fields.field2array(p.drivers.cosθs)) # zenith angle at nighttime should be > 90 degrees
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
