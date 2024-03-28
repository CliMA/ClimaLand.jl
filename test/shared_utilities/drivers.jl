using ClimaCore
using Test
using StaticArrays
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand

FT = Float32
@testset "Default model, FT = $FT" begin
    pa = ClimaLand.PrescribedAtmosphere(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        FT(1);
    )
    pr = ClimaLand.PrescribedRadiativeFluxes(FT, nothing, nothing, nothing)
    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)
    coords = (; surface = [1, 2, 3])
    @test ClimaLand.initialize_drivers(pp, nothing, coords) ==
          NamedTuple{(:P_liq,)}((zeros(FT, 3),))
    @test ClimaLand.initialize_drivers(nothing, nothing, coords) == (;)
    pa_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    pa_vals = ([zeros(FT, 3) for k in pa_keys]...,)
    @test ClimaLand.initialize_drivers(pa, nothing, coords) ==
          NamedTuple{pa_keys}(pa_vals)
    papr_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2, :SW_d, :LW_d, :θs)
    papr_vals = ([zeros(FT, 3) for k in papr_keys]...,)
    @test ClimaLand.initialize_drivers(pa, pr, coords) ==
          NamedTuple{papr_keys}(papr_vals)
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
    f = TimeVaryingInput((t) -> 10.0)
    pa = ClimaLand.PrescribedAtmosphere(f, f, f, f, f, f, f, FT(1);)
    pr = ClimaLand.PrescribedRadiativeFluxes(FT, f, f, f)
    coords = (; surface = [1])
    p = (; drivers = ClimaLand.initialize_drivers(nothing, nothing, coords))
    nothing_update! = ClimaLand.make_update_drivers(nothing, nothing)
    nothing_update!(p, 0.0)
    @test p.drivers == (;)
    p = (; drivers = ClimaLand.initialize_drivers(pa, nothing, coords))
    atmos_only_update! = ClimaLand.make_update_drivers(pa, nothing)
    atmos_only_update!(p, 0.0)
    @test p.drivers.P_liq == [FT(10)]
    @test p.drivers.P_snow == [FT(10)]
    @test p.drivers.P == [FT(10)]
    @test p.drivers.T == [FT(10)]
    @test p.drivers.q == [FT(10)]
    @test p.drivers.u == [FT(10)]
    @test p.drivers.c_co2 == [FT(4.2e-4)]

    p = (; drivers = ClimaLand.initialize_drivers(pa, pr, coords))
    update! = ClimaLand.make_update_drivers(pa, pr)
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
    @test p.drivers.θs == [FT(0)]

    p = (; drivers = ClimaLand.initialize_drivers(nothing, pr, coords))
    rad_only_update! = ClimaLand.make_update_drivers(nothing, pr)
    rad_only_update!(p, 0.0)
    @test p.drivers.SW_d == [FT(10)]
    @test p.drivers.LW_d == [FT(10)]
    @test p.drivers.θs == [FT(0)]

    liquid_precip = TimeVaryingInput((t) -> -1.0)
    pp = ClimaLand.PrescribedPrecipitation{FT}(liquid_precip)
    precip_update! = ClimaLand.make_update_drivers(pp, nothing)
    p = (; drivers = ClimaLand.initialize_drivers(pp, nothing, coords))
    precip_update!(p, 0.0)
    @test p.drivers.P_liq == [FT(-1.0)]
end
