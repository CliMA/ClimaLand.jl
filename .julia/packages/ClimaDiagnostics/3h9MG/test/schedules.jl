using Test
import SciMLBase

import Dates

import ClimaDiagnostics.Schedules:
    DivisorSchedule, EveryDtSchedule, EveryCalendarDtSchedule

import ClimaDiagnostics: time_to_date

import ClimaUtilities.TimeManager: ITime

include("TestTools.jl")

@testset "Schedules" begin
    # Test a callback that is called at every other iteration with DivisorSchedule

    called = Ref(0)
    function callback_func0(integrator)
        called[] += 1
    end

    divisor = 2
    scheduled_func = DivisorSchedule(divisor)
    @test "$scheduled_func" == "2it"

    callback_everystep = SciMLBase.DiscreteCallback(
        (_, _, integrator) -> scheduled_func(integrator),
        callback_func0,
    )

    t0 = 0.0
    tf = 1.0
    dt = 1e-3

    space = ColumnCenterFiniteDifferenceSpace()
    args, kwargs = create_problem(space; t0, tf, dt)

    expected_called = convert(Int, (tf - t0) / (divisor * dt))

    SciMLBase.solve(args...; kwargs..., callback = callback_everystep)

    @test called[] == expected_called

    # EveryDtSchedule

    called = Ref(0)
    function callback_func(integrator)
        called[] += 1
    end

    called2 = Ref(0)
    function callback_func2(integrator)
        called2[] += 1
    end

    dt_callback = 0.2
    scheduled_func = EveryDtSchedule(dt_callback)
    @test "$scheduled_func" == "0.2s"

    scheduled_func_test1 = EveryDtSchedule(dt_callback; t_last = 0.1)
    scheduled_func_test2 = EveryDtSchedule(dt_callback; t_last = 0.1)
    @test scheduled_func_test1 == scheduled_func_test2
    @test !(scheduled_func_test1 === scheduled_func_test2)

    dt_callback2 = 0.3
    t_last2 = 0.1
    scheduled_func2 = EveryDtSchedule(dt_callback2; t_last = t_last2)

    callback_dt = SciMLBase.DiscreteCallback(
        (_, _, integrator) -> scheduled_func(integrator),
        callback_func,
    )
    callback_dt2 = SciMLBase.DiscreteCallback(
        (_, _, integrator) -> scheduled_func2(integrator),
        callback_func2,
    )

    args, kwargs = create_problem(space; t0, tf, dt)

    expected_called = convert(Int, (tf - t0) / dt_callback)
    expected_called2 = convert(Int, floor((tf - t0 - t_last2) / dt_callback2))

    SciMLBase.solve(
        args...;
        kwargs...,
        callback = SciMLBase.CallbackSet(callback_dt, callback_dt2),
    )

    @test called[] == expected_called
    @test called2[] == expected_called2

    # Test EveryCalendarDtSchedule
    dt_callback_month = Dates.Month(1)
    dummyintegrator = (; t = 1.0)

    schedule0 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
    )
    @test !schedule0(dummyintegrator)

    schedule1 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 7),
    )
    @test schedule1(dummyintegrator)

    schedule2 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 12),
    )
    @test !schedule2(dummyintegrator)

    # Check precisely one month
    schedule3 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 10),
    )
    # one month from the reference date in seconds
    one_month_from_reference = Dates.value(
        Dates.Second(
            Dates.DateTime(2024, 10) + Dates.Month(1) -
            Dates.DateTime(2024, 10),
        ),
    )

    dummyintegrator_one_month = (; t = 1.0 * one_month_from_reference)
    @test schedule3(dummyintegrator_one_month)

    # Test EveryCalendarDtSchedule with ITime
    # Schedules should be the same internally
    schedule1 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 7),
    )
    schedule2 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 12),
    )

    itime_schedule1 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = ITime(0, epoch = Dates.DateTime(2024, 10)),
        date_last = ITime(0, epoch = Dates.DateTime(2024, 7)),
    )
    @test schedule1 == itime_schedule1

    itime_schedule = EveryCalendarDtSchedule(
        # DateTime(2024, 7) - DateTime(2024, 10) |> Seconds is DateTime(2024, 10)
        dt_callback_month,
        ITime(-7948800, epoch = Dates.DateTime(2024, 10)),
    )
    @test schedule1 == itime_schedule

    # Schedules should work when integrator.t is an ITime
    dummyitimeintegrator = (; t = ITime(1.0, epoch = Dates.DateTime(2024, 10)))
    schedule1 = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2024, 10),
        date_last = Dates.DateTime(2024, 7),
    )
    @test schedule1(dummyitimeintegrator)
    @test !schedule2(dummyitimeintegrator)

    called_month = Ref(0)
    function callback_func_month(integrator)
        called_month[] += 1
    end

    # Test with start_date read from integrator.p
    scheduled_func_month = EveryCalendarDtSchedule(
        dt_callback_month,
        start_date = Dates.DateTime(2000),
    )
    @test "$scheduled_func_month" == "1M"

    callback_dt_month = SciMLBase.DiscreteCallback(
        (_, _, integrator) -> scheduled_func_month(integrator),
        callback_func_month,
    )
    t0 = 0.0
    tf = 2 * 32 * 86400  # At least two months
    dt = tf / 100
    args, kwargs = create_problem(space; t0, tf, dt)

    # See https://github.com/JuliaLang/julia/issues/55406
    dt_callback_month_in_s = 86400 * Dates.days(dt_callback_month)
    expected_called_month =
        convert(Int, floor((tf - t0) / dt_callback_month_in_s))

    SciMLBase.solve(
        args...;
        kwargs...,
        callback = SciMLBase.CallbackSet(callback_dt_month),
    )

    @test called_month[] == expected_called_month
end

@testset "time_to_date" begin
    start_date = Dates.DateTime(2024, 1, 1, 0, 0, 0)

    # Test with whole seconds
    @test time_to_date(0.0, start_date) == start_date
    @test time_to_date(60.0, start_date) == start_date + Dates.Minute(1)

    # Test with fractional seconds
    @test time_to_date(0.5, start_date) == start_date + Dates.Millisecond(500)
    @test time_to_date(1.25, start_date) ==
          start_date + Dates.Second(1) + Dates.Millisecond(250)

    # Test with negative times
    @test time_to_date(-1.0, start_date) == start_date - Dates.Second(1)
    @test time_to_date(-0.5, start_date) == start_date - Dates.Millisecond(500)

    # Test with ITime
    @test time_to_date(
        ITime(1, epoch = start_date),
        start_date + Dates.Second(42),
    ) == start_date + Dates.Second(1)
    @test time_to_date(
        ITime(-1, epoch = start_date),
        start_date + Dates.Second(42),
    ) == start_date - Dates.Second(1)
end
