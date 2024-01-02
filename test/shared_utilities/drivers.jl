using ClimaCore
using Test
using StaticArrays
using ClimaLSM

FT = Float32
@testset "Default model, FT = $FT" begin
    pa = ClimaLSM.PrescribedAtmosphere(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        FT(1);
    )
    pr = ClimaLSM.PrescribedRadiativeFluxes(FT, nothing, nothing, nothing)
    coords = (; surface = [1, 2, 3])
    @test ClimaLSM.initialize_drivers(nothing, nothing, coords) == (;)
    pa_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2)
    pa_vals = ([zeros(FT, 3) for k in pa_keys]...,)
    @test ClimaLSM.initialize_drivers(pa, nothing, coords) ==
          NamedTuple{pa_keys}(pa_vals)
    papr_keys = (:P_liq, :P_snow, :T, :P, :u, :q, :c_co2, :SW_d, :LW_d, :θs)
    papr_vals = ([zeros(FT, 3) for k in papr_keys]...,)
    @test ClimaLSM.initialize_drivers(pa, pr, coords) ==
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
    cond = ClimaLSM.update_condition(updateat)
    @test cond(nothing, t0 - dt, nothing) == false
    @test cond(nothing, tf, nothing) == false
    @test cond(nothing, t0 + 0.5 * dt, nothing) == true
    @test cond(nothing, updateat[end], nothing) == true
    updatefunc = (drivers, t) -> drivers.q .= [t]
    affect! = ClimaLSM.DriverAffect(updateat, updatefunc)
    @test affect!.updateat == updateat
    @test affect!.updatefunc == updatefunc

    p_init = (; drivers = (; q = [FT(-1)]))
    integrator = Integrator{FT}(t0, dt, p_init)
    cb = (; affect! = affect!)
    ClimaLSM.driver_initialize(cb, nothing, t0, integrator)
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
    f = (t) -> 10.0
    pa = ClimaLSM.PrescribedAtmosphere(f, f, f, f, f, f, f, FT(1);)
    pr = ClimaLSM.PrescribedRadiativeFluxes(FT, f, f, f)
    coords = (; surface = [1])
    drivers = ClimaLSM.initialize_drivers(nothing, nothing, coords)
    nothing_update! = ClimaLSM.make_update_drivers(nothing, nothing)
    nothing_update!(drivers, 0.0)
    @test drivers == (;)
    drivers = ClimaLSM.initialize_drivers(pa, nothing, coords)
    atmos_only_update! = ClimaLSM.make_update_drivers(pa, nothing)
    atmos_only_update!(drivers, 0.0)
    @test drivers.P_liq == [FT(10)]
    @test drivers.P_snow == [FT(10)]
    @test drivers.P == [FT(10)]
    @test drivers.T == [FT(10)]
    @test drivers.q == [FT(10)]
    @test drivers.u == [FT(10)]
    @test drivers.c_co2 == [FT(4.2e-4)]

    drivers = ClimaLSM.initialize_drivers(pa, pr, coords)
    update! = ClimaLSM.make_update_drivers(pa, pr)
    update!(drivers, 0.0)
    @test drivers.P_liq == [FT(10)]
    @test drivers.P_snow == [FT(10)]
    @test drivers.P == [FT(10)]
    @test drivers.T == [FT(10)]
    @test drivers.q == [FT(10)]
    @test drivers.u == [FT(10)]
    @test drivers.c_co2 == [FT(4.2e-4)]
    @test drivers.SW_d == [FT(10)]
    @test drivers.LW_d == [FT(10)]
    @test drivers.θs == [FT(0)]

    drivers = ClimaLSM.initialize_drivers(nothing, pr, coords)
    rad_only_update! = ClimaLSM.make_update_drivers(nothing, pr)
    rad_only_update!(drivers, 0.0)
    @test drivers.SW_d == [FT(10)]
    @test drivers.LW_d == [FT(10)]
    @test drivers.θs == [FT(0)]
end
