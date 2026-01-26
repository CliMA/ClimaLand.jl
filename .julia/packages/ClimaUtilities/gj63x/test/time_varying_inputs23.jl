using Test
using Dates
using Artifacts

import ClimaUtilities
import ClimaUtilities: DataHandling
import ClimaUtilities: TimeVaryingInputs
import ClimaUtilities.TimeManager: ITime

import ClimaCore: Domains, Geometry, Fields, Meshes, Topologies, Spaces
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

import Interpolations
import NCDatasets
import ClimaCoreTempestRemap

const context = ClimaComms.context()
ClimaComms.init(context)
const singleton_cpu_context =
    ClimaComms.SingletonCommsContext(ClimaComms.device())

include("TestTools.jl")

@testset "InterpolatingTimeVaryingInput23D" begin
    PATH = joinpath(artifact"era5_example", "era5_t2m_sp_u10n_20210101.nc")
    regridder_types = (:InterpolationsRegridder, :TempestRegridder)

    # Tempestremap is single threaded and CPU only
    if context isa ClimaComms.MPICommsContext ||
       ClimaComms.device() isa ClimaComms.CUDADevice ||
       Sys.iswindows()
        regridder_types = (:InterpolationsRegridder,)
    end
    for FT in (Float32, Float64)
        for regridder_type in regridder_types
            if regridder_type == :TempestRegridder
                target_spaces = (make_spherical_space(FT; context).horizontal,)
            else
                target_spaces = (
                    make_spherical_space(FT; context).horizontal,
                    make_regional_space(FT; context).horizontal,
                )
            end
            for target_space in target_spaces
                data_handler = DataHandling.DataHandler(
                    PATH,
                    "u10n",
                    target_space;
                    start_date = Dates.DateTime(2021, 1, 1, 1),
                    regridder_type,
                )

                # Test constructor with multiple variables
                regridder_type_interp = :InterpolationsRegridder
                compose_function = (x, y) -> x + y
                input_multiple_vars = TimeVaryingInputs.TimeVaryingInput(
                    PATH,
                    ["u10n", "t2m"],
                    target_space;
                    start_date = Dates.DateTime(2021, 1, 1, 1),
                    t_start = 0.0,
                    regridder_type = regridder_type_interp,
                    compose_function,
                )

                input_nearest = TimeVaryingInputs.TimeVaryingInput(
                    PATH,
                    "u10n",
                    target_space;
                    start_date = Dates.DateTime(2021, 1, 1, 1),
                    regridder_type,
                    method = TimeVaryingInputs.NearestNeighbor(),
                )

                input_nearest_flat = TimeVaryingInputs.TimeVaryingInput(
                    PATH,
                    "u10n",
                    target_space;
                    start_date = Dates.DateTime(2021, 1, 1, 1),
                    regridder_type,
                    method = TimeVaryingInputs.NearestNeighbor(
                        TimeVaryingInputs.Flat(),
                    ),
                )

                input_nearest_periodic_calendar =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        start_date = Dates.DateTime(2021, 1, 1, 1),
                        regridder_type,
                        method = TimeVaryingInputs.NearestNeighbor(
                            TimeVaryingInputs.PeriodicCalendar(),
                        ),
                    )

                input_linear_periodic_calendar =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        start_date = Dates.DateTime(2021, 1, 1, 1),
                        regridder_type,
                        method = TimeVaryingInputs.LinearInterpolation(
                            TimeVaryingInputs.PeriodicCalendar(),
                        ),
                    )

                # Repeat one day over and over
                period = Dates.Day(1)
                repeat_date = Dates.DateTime(2021, 1, 1, 1)

                # Note, the data is define on the day of 2021/1/1, so Dates.DateTime(2021, 1, 1,
                # 1) is one hour in. This allows us to test misaligned periods
                input_nearest_periodic_calendar_date =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        start_date = Dates.DateTime(2021, 1, 1, 1),
                        regridder_type,
                        method = TimeVaryingInputs.NearestNeighbor(
                            TimeVaryingInputs.PeriodicCalendar(
                                period,
                                repeat_date,
                            ),
                        ),
                    )

                input_linear_periodic_calendar_date =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        start_date = Dates.DateTime(2021, 1, 1, 1),
                        regridder_type,
                        method = TimeVaryingInputs.LinearInterpolation(
                            TimeVaryingInputs.PeriodicCalendar(
                                period,
                                repeat_date,
                            ),
                        ),
                    )

                # Since TimeVaryingInputs internally use dates, we don't need a
                # start date
                input_nearest_no_start_date =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        regridder_type,
                        method = TimeVaryingInputs.NearestNeighbor(),
                    )

                input_linear_periodic_calendar_no_start_date =
                    TimeVaryingInputs.TimeVaryingInput(
                        PATH,
                        "u10n",
                        target_space;
                        regridder_type,
                        method = TimeVaryingInputs.LinearInterpolation(
                            TimeVaryingInputs.PeriodicCalendar(),
                        ),
                    )

                @test 0.0 in input_nearest
                @test !(82801.0 in input_nearest)
                @test Dates.DateTime(2021, 1, 1, 0) in input_nearest
                @test !(Dates.DateTime(2021, 1, 2, 0) in input_nearest)

                available_times = DataHandling.available_times(data_handler)
                available_dates = DataHandling.available_dates(data_handler)
                dest = Fields.zeros(target_space)

                # Time outside of range
                for t in (
                    FT(-40000),
                    Dates.DateTime(2020),
                    ITime(0, epoch = Dates.DateTime(2020)),
                )
                    @test_throws ErrorException TimeVaryingInputs.evaluate!(
                        dest,
                        input_nearest,
                        t,
                    )
                    if !(t isa AbstractFloat)
                        @test_throws ErrorException TimeVaryingInputs.evaluate!(
                            dest,
                            input_nearest_no_start_date,
                            t,
                        )
                    end
                end

                # We are testing NearestNeighbor, so we can just have to check if the fields agree

                # Left nearest point
                target_time = available_times[10] + 1
                target_date = available_dates[10] + Second(1)
                for t in (
                    target_time,
                    target_date,
                    ITime(1, epoch = available_dates[10]),
                )
                    TimeVaryingInputs.evaluate!(dest, input_nearest, t)

                    # We use isequal to handle NaNs
                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[10],
                                ),
                            ),
                        ),
                    )

                    if !(t isa AbstractFloat)
                        TimeVaryingInputs.evaluate!(
                            dest,
                            input_nearest_no_start_date,
                            t,
                        )
                        @test isequal(
                            Array(parent(dest)),
                            Array(
                                parent(
                                    DataHandling.regridded_snapshot(
                                        data_handler,
                                        available_times[10],
                                    ),
                                ),
                            ),
                        )
                    end
                end

                # Right nearest point
                target_time = available_times[9] - 1
                target_date = available_dates[9] - Second(1)
                for t in (
                    target_time,
                    target_date,
                    ITime(-1, epoch = available_dates[9]),
                )
                    TimeVaryingInputs.evaluate!(dest, input_nearest, t)

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[9],
                                ),
                            ),
                        ),
                    )

                    if !(t isa AbstractFloat)
                        TimeVaryingInputs.evaluate!(
                            dest,
                            input_nearest_no_start_date,
                            t,
                        )
                        @test isequal(
                            Array(parent(dest)),
                            Array(
                                parent(
                                    DataHandling.regridded_snapshot(
                                        data_handler,
                                        available_times[9],
                                    ),
                                ),
                            ),
                        )
                    end
                end

                # On node
                target_time = available_times[11]
                target_date = available_dates[11]
                for t in (
                    target_time,
                    target_date,
                    ITime(0, epoch = available_dates[11]),
                )
                    TimeVaryingInputs.evaluate!(dest, input_nearest, t)

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[11],
                                ),
                            ),
                        ),
                    )

                    if !(t isa AbstractFloat)
                        TimeVaryingInputs.evaluate!(
                            dest,
                            input_nearest_no_start_date,
                            t,
                        )
                        @test isequal(
                            Array(parent(dest)),
                            Array(
                                parent(
                                    DataHandling.regridded_snapshot(
                                        data_handler,
                                        available_times[11],
                                    ),
                                ),
                            ),
                        )
                    end
                end

                # Flat left
                target_time = available_times[begin] - 1
                target_date = available_dates[begin] - Second(1)
                for t in (
                    target_time,
                    target_date,
                    ITime(-1, epoch = available_dates[begin]),
                )
                    TimeVaryingInputs.evaluate!(dest, input_nearest_flat, t)

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[begin],
                                ),
                            ),
                        ),
                    )
                end

                # Flat right
                target_time = available_times[end] + 1
                target_date = available_dates[end] + Second(1)
                for t in (
                    target_time,
                    target_date,
                    ITime(1, epoch = available_dates[end]),
                )
                    TimeVaryingInputs.evaluate!(dest, input_nearest_flat, t)

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[end],
                                ),
                            ),
                        ),
                    )
                end

                # Nearest periodic calendar
                dt = available_times[2] - available_times[1]
                target_time = available_times[end] + 0.1dt
                d_date = available_dates[2] - available_dates[1]
                target_date = available_dates[end] + 0.1d_date
                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(
                        dest,
                        input_nearest_periodic_calendar,
                        t,
                    )

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[end],
                                ),
                            ),
                        ),
                    )
                end

                # With date
                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(
                        dest,
                        input_nearest_periodic_calendar_date,
                        t,
                    )

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[end],
                                ),
                            ),
                        ),
                    )
                end

                dt = available_times[2] - available_times[1]
                target_time = available_times[end] + 0.6dt
                d_date = available_dates[2] - available_dates[1]
                target_date = available_dates[end] + 0.6d_date
                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(
                        dest,
                        input_nearest_periodic_calendar,
                        t,
                    )

                    @test isequal(
                        Array(parent(dest)),
                        Array(
                            parent(
                                DataHandling.regridded_snapshot(
                                    data_handler,
                                    available_times[begin],
                                ),
                            ),
                        ),
                    )
                end

                # Now testing LinearInterpolation
                input_linear = TimeVaryingInputs.TimeVaryingInput(data_handler)

                # Time outside of range
                @test_throws ErrorException TimeVaryingInputs.evaluate!(
                    dest,
                    input_linear,
                    FT(-40000),
                )
                @test_throws ErrorException TimeVaryingInputs.evaluate!(
                    dest,
                    input_linear,
                    Dates.DateTime(2021, 1, 1) + Second(-40000),
                )

                left_value = DataHandling.regridded_snapshot(
                    data_handler,
                    available_times[10],
                )
                right_value = DataHandling.regridded_snapshot(
                    data_handler,
                    available_times[11],
                )

                target_time = available_times[10] + 30
                left_time = available_times[10]
                right_time = available_times[11]

                target_date = available_dates[10] + Second(30)

                expected = Fields.zeros(target_space)
                expected .=
                    left_value .+
                    (target_time - left_time) / (right_time - left_time) .*
                    (right_value .- left_value)

                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(dest, input_linear, t)
                    @test parent(dest) ≈ parent(expected)
                end

                # LinearInterpolation with PeriodicCalendar
                time_delta = 0.1dt
                target_time = available_times[end] + time_delta

                left_value = DataHandling.regridded_snapshot(
                    data_handler,
                    available_times[end],
                )
                right_value = DataHandling.regridded_snapshot(
                    data_handler,
                    available_times[begin],
                )

                time_delta = 0.1dt
                target_time = available_times[end] + time_delta
                left_time = available_times[10]
                right_time = available_times[11]

                expected = Fields.zeros(target_space)
                expected .=
                    left_value .+ time_delta / dt .* (right_value .- left_value)

                target_date = available_dates[end] + 0.1d_date
                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(
                        dest,
                        input_linear_periodic_calendar,
                        t,
                    )
                    @test parent(dest) ≈ parent(expected)
                    if !(t isa AbstractFloat)
                        TimeVaryingInputs.evaluate!(
                            dest,
                            input_linear_periodic_calendar_no_start_date,
                            t,
                        )
                        @test parent(dest) ≈ parent(expected)
                    end
                end

                # With offset of one period
                target_time += 86400
                target_date += Second(86400)
                for t in
                    (target_time, target_date, ITime(0, epoch = target_date))
                    TimeVaryingInputs.evaluate!(
                        dest,
                        input_linear_periodic_calendar_date,
                        t,
                    )
                    @test parent(dest) ≈ parent(expected)
                end

                close(input_multiple_vars)
                close(input_nearest)
                close(input_linear)
                close(input_nearest_flat)
                close(input_nearest_periodic_calendar)
                close(input_linear_periodic_calendar)
                close(input_nearest_no_start_date)
                close(input_linear_periodic_calendar_no_start_date)
            end
        end
    end
end
