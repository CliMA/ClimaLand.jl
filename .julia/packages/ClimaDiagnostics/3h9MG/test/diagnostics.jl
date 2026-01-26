using Test
import ClimaDiagnostics
import ClimaDiagnostics.Schedules
import ClimaDiagnostics.ScheduledDiagnostics
import ClimaDiagnostics.Writers

import SciMLBase

include("TestTools.jl")

@testset "Diagnostics" begin
    t0 = 0.0
    tf = 1.0
    dt = 0.1

    space = ColumnCenterFiniteDifferenceSpace()
    Y = ClimaCore.Fields.FieldVector(; my_var = ones(space))
    p = (; tau = -0.1)

    function exp_tendency!(dY, Y, p, t)
        @. dY.my_var = p.tau * Y.my_var
    end

    prob = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(T_exp! = exp_tendency!),
        Y,
        (t0, tf),
        p,
    )
    algo = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

    function computet!(out, u, p, t)
        if isnothing(out)
            return [t]
        else
            out .= t
            return nothing
        end
    end

    # A simple diagnostic that just saves the time at every timestep
    simple_var = ClimaDiagnostics.DiagnosticVariable(;
        compute! = computet!,
        short_name = "YO",
        long_name = "YO YO",
    )

    dict_writer = Writers.DictWriter()

    # Let's start with a simple diagnostic that is computed and output
    # at every time step (the default)
    diagnostic_every_step = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
    )
    short_name = ScheduledDiagnostics.output_short_name(diagnostic_every_step)

    diagnostic_handler = ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic_every_step],
        Y,
        p,
        t0;
        dt,
    )

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    prob = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(T_exp! = exp_tendency!),
        Y,
        (t0, tf),
        p,
    )
    algo = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

    SciMLBase.solve(prob, algo, dt = dt, callback = diag_cb)

    @test length(keys(dict_writer.dict[short_name])) ==
          convert(Int, 1 + (tf - t0) / dt)
    @test dict_writer[short_name][t0] == [t0]
    @test dict_writer[short_name][tf] == [tf]

    # Now test accumulation and average
    dict_writer = Writers.DictWriter()

    diagnostic_accumulate_every_step = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
        reduction_time_func = +,
        output_schedule_func = Schedules.DivisorSchedule(5),
    )
    short_name =
        ScheduledDiagnostics.output_short_name(diagnostic_accumulate_every_step)

    diagnostic_handler = ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic_accumulate_every_step],
        Y,
        p,
        t0;
        dt,
    )

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    prob = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(T_exp! = exp_tendency!),
        Y,
        (t0, tf),
        p,
    )
    algo = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

    SciMLBase.solve(prob, algo, dt = dt, callback = diag_cb)

    @test length(keys(dict_writer.dict[short_name])) ==
          convert(Int, (tf - t0) / 5dt)

    @test dict_writer[short_name][0.5][] ≈ sum(t0:dt:(tf / 2))
    @test dict_writer[short_name][1.0][] ≈ sum((tf / 2 + dt):dt:tf)

    # Incompatible timestep in output
    diagnostic_incompatible_timestep = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
        reduction_time_func = +,
        output_schedule_func = Schedules.EveryDtSchedule(1.5dt),
    )
    @test_throws ErrorException ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic_incompatible_timestep],
        Y,
        p,
        t0;
        dt,
    )
    # Incompatible timestep in compute
    diagnostic_incompatible_timestep = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
        reduction_time_func = +,
        compute_schedule_func = Schedules.EveryDtSchedule(1.5dt),
    )
    @test_throws ErrorException ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic_incompatible_timestep],
        Y,
        p,
        t0;
        dt,
    )

    # Test duplicated diagnostics
    diagnostic1 = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
    )
    diagnostic2 = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
    )
    diagnostic3 = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dict_writer,
        reduction_time_func = +,
    )

    @test_logs (
        :warn,
        "Given list of diagnostics contains duplicates, removing them",
    ) handler_dup = ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic1, diagnostic2, diagnostic3],
        Y,
        p,
        t0;
        dt,
    )

    handler_dup = ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic1, diagnostic2, diagnostic3],
        Y,
        p,
        t0;
        dt,
    )

    @test length(handler_dup.scheduled_diagnostics) == 2

end
