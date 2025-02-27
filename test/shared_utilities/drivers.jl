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

FT = Float32
@testset "Default model, FT = $FT" begin
    earth_param_set = LP.LandParameters(FT)
    pa = ClimaLand.PrescribedAtmosphere(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        FT(1),
        earth_param_set,
    )
    pr = ClimaLand.PrescribedRadiativeFluxes(FT, nothing, nothing, nothing)
    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)
    coords = (; surface = [1, 2, 3])
    @test ClimaLand.initialize_drivers((pp,), coords) ==
          NamedTuple{(:P_liq,)}((zeros(FT, 3),))
    @test ClimaLand.initialize_drivers((), coords) == (;)
    pa_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    zero_thermal_state = map(
        _ -> ClimaCore.RecursiveApply.rzero(Thermodynamics.PhaseEquil{FT}),
        coords.surface,
    )
    pa_vals = ([zeros(FT, 3) for k in pa_keys]...,)
    all_pa_keys = (pa_keys..., :thermal_state)
    all_pa_vals = (pa_vals..., zero_thermal_state)
    @test ClimaLand.initialize_drivers((pa,), coords) ==
          NamedTuple{all_pa_keys}(all_pa_vals)
    pr_keys = (:SW_d, :LW_d, :cosθs, :frac_diff)
    pr_vals = ([zeros(FT, 3) for k in pr_keys]...,)
    all_papr_keys = (pa_keys..., :thermal_state, pr_keys...)
    all_papr_vals = (pa_vals..., zero_thermal_state, pr_vals...)
    @test ClimaLand.initialize_drivers((pa, pr), coords) ==
          NamedTuple{all_papr_keys}(all_papr_vals)
end

## Driver update callback tests
@testset "driver callback cond, affect! test, FT = $FT" begin
    mutable struct Integrator{FT}
        t::Any
        dt::Any
        p::NamedTuple
    end

    function upd_integrator(integrator::Integrator{FT}, dt) where {FT}
        integrator.t += dt
    end

    t0 = Float64(0)
    dt = Float64(1)
    tf = Float64(100)
    t_range = collect(t0:dt:tf)
    nsteps = Int((tf - t0) / dt)
    # specify when to update in the callback
    updateat = collect((t0 + dt * 0.5):(5.5 * dt):tf)

    # set up components of callback
    cond = ClimaLand.update_condition(updateat)
    @test cond(nothing, t0 - dt, nothing) == false
    @test cond(nothing, tf, nothing) == false
    @test cond(nothing, t0 + 0.5 * dt, nothing) == true
    @test cond(nothing, updateat[end], nothing) == true
    updatefunc = (p, t) -> p.drivers.q .= [t]
    affect! = ClimaLand.DriverAffect(updateat, updatefunc)
    @test affect!.updateat == updateat
    @test affect!.updatefunc == updatefunc

    p_init = (; drivers = (; q = [FT(-1)]))
    integrator = Integrator{FT}(t0, dt, p_init)
    cb = (; affect! = affect!)
    ClimaLand.driver_initialize(cb, nothing, t0, integrator)
    @test integrator.p.drivers.q == [t0]
    # simulate callback behavior
    for i in 1:nsteps
        t = integrator.t
        if cond(0, t, 0)
            next_t = first(affect!.updateat)
            affect!(integrator)
            @test integrator.p.drivers.q == [next_t]
        end
        upd_integrator(integrator, dt)
    end
end

@testset "Driver update functions" begin
    earth_param_set = LP.LandParameters(FT)
    f = TimeVaryingInput((t) -> 10.0)
    pa = ClimaLand.PrescribedAtmosphere(
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        FT(1),
        earth_param_set,
    )
    pr = ClimaLand.PrescribedRadiativeFluxes(FT, f, f, f)
    coords = (; surface = [1], subsurface = [1, 2])
    p = (; drivers = ClimaLand.initialize_drivers((), coords))
    nothing_update! = ClimaLand.make_update_drivers(())
    nothing_update!(p, 0.0)
    @test p.drivers == (;)
    p = (; drivers = ClimaLand.initialize_drivers((pa,), coords))
    atmos_only_update! = ClimaLand.make_update_drivers((pa,))
    atmos_only_update!(p, 0.0)
    @test p.drivers.P_liq == [FT(10)]
    @test p.drivers.P_snow == [FT(10)]
    @test p.drivers.P == [FT(10)]
    @test p.drivers.T == [FT(10)]
    @test p.drivers.q == [FT(10)]
    @test p.drivers.u == [FT(10)]
    @test p.drivers.c_co2 == [FT(4.2e-4)]

    p = (; drivers = ClimaLand.initialize_drivers((pa, pr), coords))
    update! = ClimaLand.make_update_drivers((pa, pr))
    update!(p, 0.0)
    @test p.drivers.P_liq == [FT(10)]
    @test p.drivers.P_snow == [FT(10)]
    @test p.drivers.P == [FT(10)]
    @test p.drivers.T == [FT(10)]
    @test p.drivers.q == [FT(10)]
    @test p.drivers.u == [FT(10)]
    @test p.drivers.c_co2 == [FT(4.2e-4)]
    @test p.drivers.SW_d == [FT(10)]
    @test p.drivers.LW_d == [FT(10)]
    @test all(isnan.(p.drivers.cosθs))
    @test all(isnan.(p.drivers.frac_diff))

    p = (; drivers = ClimaLand.initialize_drivers((pr,), coords))
    rad_only_update! = ClimaLand.make_update_drivers((pr,))
    rad_only_update!(p, 0.0)
    @test p.drivers.SW_d == [FT(10)]
    @test p.drivers.LW_d == [FT(10)]
    @test all(isnan.(p.drivers.cosθs))
    @test all(isnan.(p.drivers.frac_diff))

    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)
    precip_update! = ClimaLand.make_update_drivers((pp,))
    p = (; drivers = ClimaLand.initialize_drivers((pp,), coords))
    precip_update!(p, 0.0)
    @test p.drivers.P_liq == [FT(-1.0)]


    soc = TimeVaryingInput((t) -> 1.0)
    soc_driver = ClimaLand.PrescribedSoilOrganicCarbon{FT}(soc)
    soc_update! = ClimaLand.make_update_drivers((soc_driver,))
    soc_p = (; drivers = ClimaLand.initialize_drivers((soc_driver,), coords))
    soc_update!(soc_p, 0.0)
    @test soc_p.drivers.soc == [FT(1.0), FT(1.0)]
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
