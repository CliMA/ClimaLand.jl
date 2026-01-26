using Test
using Dates
using Artifacts

import ClimaUtilities
import ClimaUtilities: DataHandling
import ClimaUtilities: TimeVaryingInputs
using ClimaUtilities.TimeManager

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

@testset "Analytic TimeVaryingInput" begin
    fun = (x) -> 2x
    input = TimeVaryingInputs.TimeVaryingInput(fun)

    FT = Float32

    # Prepare a field
    domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(0),
        Geometry.ZPoint{FT}(5),
        boundary_names = (:bottom, :top),
    )
    mesh = Meshes.IntervalMesh(domain; nelems = 10)
    topology = Topologies.IntervalTopology(singleton_cpu_context, mesh)

    column_space = Spaces.CenterFiniteDifferenceSpace(topology)
    column_field = Fields.zeros(column_space)

    TimeVaryingInputs.evaluate!(column_field, input, 10.0)
    @test Array(parent(column_field))[1] == fun(10.0)

    # Check with args and kwargs
    fun2 = (x, y; z) -> 2x * y * z
    input2 = TimeVaryingInputs.TimeVaryingInput(fun2)

    TimeVaryingInputs.evaluate!(column_field, input2, 10.0, 20.0; z = 30.0)
    @test Array(parent(column_field))[1] == fun2(10.0, 20.0; z = 30.0)
end

@testset "InterpolatingTimeVaryingInput0D" begin
    # Check times not sorted
    xs = [1.0, 0.0]
    ys = [1.0, 2.0]

    @test_throws ErrorException TimeVaryingInputs.TimeVaryingInput(xs, ys)

    # Test with PeriodicCalendar with a given period
    # This is not allowed for 1D data because we don't have dates
    @test_throws ErrorException TimeVaryingInputs.TimeVaryingInput(
        sort(xs),
        ys,
        method = TimeVaryingInputs.NearestNeighbor(
            TimeVaryingInputs.PeriodicCalendar(Month(1), Date(2024)),
        ),
    )

    # Test with LinearPeriodFillingInterpolation with a given period
    # This is not allowed for 1D data because we don't have dates
    @test_throws ErrorException TimeVaryingInputs.TimeVaryingInput(
        sort(xs),
        ys,
        method = TimeVaryingInputs.LinearPeriodFillingInterpolation(),
    )

    # Test PeriodicCalendar with non simple duration
    @test_throws ErrorException TimeVaryingInputs.PeriodicCalendar(
        Month(2),
        Date(2024),
    )

    # Test LinearPeriodFillingInterpolation with non simple duration
    @test_throws ErrorException TimeVaryingInputs.LinearPeriodFillingInterpolation(
        Month(2),
        TimeVaryingInputs.Throw(),
    )

    # Test LinearPeriodFillingInterpolation with PeriodicCalendar (they are incompatible)
    @test_throws ErrorException TimeVaryingInputs.LinearPeriodFillingInterpolation(
        Year(1),
        TimeVaryingInputs.PeriodicCalendar(),
    )

    for FT in (Float32, Float64)
        # Prepare spaces/fields

        domain = Domains.IntervalDomain(
            Geometry.ZPoint{FT}(0),
            Geometry.ZPoint{FT}(5),
            boundary_names = (:bottom, :top),
        )
        mesh = Meshes.IntervalMesh(domain; nelems = 10)
        topology = Topologies.IntervalTopology(singleton_cpu_context, mesh)

        column_space = Spaces.CenterFiniteDifferenceSpace(topology)
        point_space = Spaces.level(column_space, 1)
        column_field = Fields.zeros(column_space)
        point_field = Fields.zeros(point_space)

        start_date = Dates.DateTime(2014)
        times_ft = collect(FT(0):FT(0.5):FT(100))
        vals = sin.(times_ft)
        itime_no_epoch = [promote(ITime.(times_ft)...)...]
        add_test_epoch = t -> promote(t, ITime(0; epoch = start_date))[1]
        itime_with_epoch = [promote(add_test_epoch.(itime_no_epoch)...)...]

        # tuple where first element is vector used to create test TimeVaryingInputs
        # second element converts floats into the desired test inputs
        # third element converts the elements of the times vector to the desired test inputs
        for (times, ft_to_input, time_to_input) in (
            (times_ft, identity, identity), # FT TVI, inputs as FT
            (times_ft, ITime, ITime), # FT TVI, inputs as ITime no date
            (
                times_ft,
                t -> add_test_epoch(ITime(t)),
                t -> add_test_epoch(ITime(t)),
            ), # FT TVI, inputs as ITime with date
            (itime_no_epoch, ITime, identity), # ITime TVI no dates, inputs as ITime (no dates)
            (itime_no_epoch, identity, float), # ITime TVI no dates, inputs as floats
            (itime_no_epoch, t -> add_test_epoch(ITime(t)), add_test_epoch), # ITime TVI no dates, inputs as ITime with dates
            (itime_with_epoch, t -> add_test_epoch(ITime(t)), identity), # ITime TVI with dates, inputs as ITime with dates
            (
                itime_with_epoch,
                t -> ITime(t),
                t -> ITime(t.counter; period = t.period),
            ), # ITime TVI with dates, inputs as ITime with no dates
            (itime_with_epoch, identity, float), # ITime TVI with dates, inputs as floats
            (itime_with_epoch, t -> start_date + Second(t), DateTime), # ITime TVI with dates, inputs as datetimes
        )
            dt = times[2] - times[1]
            # Nearest neighbor interpolation
            input = TimeVaryingInputs.TimeVaryingInput(
                times,
                vals;
                method = TimeVaryingInputs.NearestNeighbor(),
            )

            # Nearest neighbor interpolation with Flat
            input_clamp = TimeVaryingInputs.TimeVaryingInput(
                times,
                vals;
                method = TimeVaryingInputs.NearestNeighbor(
                    TimeVaryingInputs.Flat(),
                ),
            )

            # Nearest neighbor interpolation with PeriodicCalendar
            input_periodic_calendar = TimeVaryingInputs.TimeVaryingInput(
                times,
                vals;
                method = TimeVaryingInputs.NearestNeighbor(
                    TimeVaryingInputs.PeriodicCalendar(),
                ),
            )

            # Linear interpolation with PeriodicCalendar
            input_periodic_calendar_linear = TimeVaryingInputs.TimeVaryingInput(
                times,
                vals;
                method = TimeVaryingInputs.LinearInterpolation(
                    TimeVaryingInputs.PeriodicCalendar(),
                ),
            )

            # Test extrapolation_bc
            @test TimeVaryingInputs.extrapolation_bc(
                TimeVaryingInputs.NearestNeighbor(),
            ) == TimeVaryingInputs.Throw()

            # Test in
            if ft_to_input(FT(3.0)) isa eltype(times)
                @test FT(3.0) in input
                @test !(FT(-3.0) in input)
            end

            # Test with different types of spaces
            for dest in (point_field, column_field)
                # Time outside of range
                @test_throws ErrorException TimeVaryingInputs.evaluate!(
                    dest,
                    input,
                    ft_to_input(FT(-4)),
                )

                # Time outside of range with Flat left
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_clamp,
                    ft_to_input(FT(-4)),
                )
                @test Array(parent(dest))[1] == vals[begin]

                # Time outside of range with Flat right
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_clamp,
                    ft_to_input(FT(400)),
                )
                @test Array(parent(dest))[1] == vals[end]

                # Time outside of range with PeriodicCalendar
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar,
                    time_to_input(times[begin]),
                )
                @test Array(parent(dest))[1] == vals[begin]

                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar,
                    time_to_input(times[end]),
                )
                @test Array(parent(dest))[1] == vals[end]

                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar,
                    time_to_input(times[end] + 10dt),
                )
                @test Array(parent(dest))[1] == vals[begin + 9]

                # Check after times[end] but before times[end] + 0.5dt, should lead be
                # equivalent to times[end]
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar,
                    time_to_input(times[end] + 0.3dt),
                )
                @test Array(parent(dest))[1] == vals[end]

                # Check after times[end] and after times[end] + 0.5dt, should lead be equivalent
                # to times[begin]
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar,
                    time_to_input(times[end] + 0.8dt),
                )
                @test Array(parent(dest))[1] == vals[begin]

                # Linear interpolation

                TimeVaryingInputs.evaluate!(
                    dest,
                    input,
                    time_to_input(times[10]),
                )

                @test Array(parent(dest))[1] == vals[10]

                # Linear interpolation
                input = TimeVaryingInputs.TimeVaryingInput(times, vals)

                TimeVaryingInputs.evaluate!(dest, input, ft_to_input(FT(10)))

                if eltype(time) <: Number
                    # searchsortedfirst only needs to work with floats, as TVI0d will
                    # internally handle conversions
                    index = searchsortedfirst(times, 0.1)
                    @test times[index - 1] <= 0.1 <= times[index]
                    expected =
                        vals[index - 1] +
                        (vals[index] - vals[index - 1]) /
                        (times[index] - times[index - 1]) *
                        (0.1 - times[index - 1])

                    @test Array(parent(dest))[1] ≈ expected
                end

                # Check edge case
                TimeVaryingInputs.evaluate!(dest, input, ft_to_input(FT(0.0)))

                @test Array(parent(dest))[1] ≈ 0.0

                # Linear interpolation with PeriodicCalendar
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[10]),
                )
                @test Array(parent(dest))[1] == vals[10]

                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[1]),
                )
                @test Array(parent(dest))[1] == vals[1]

                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[end]),
                )
                @test Array(parent(dest))[1] == vals[end]

                # t_end + dt is equivalent to t_init
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[end] + dt),
                )
                @test Array(parent(dest))[1] ≈ vals[1]

                # t_end + 2dt is equivalent to t_init + dt
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[end] + 2dt),
                )
                @test Array(parent(dest))[1] ≈ vals[2]

                # In between t_end and t_init
                expected =
                    vals[end] +
                    (vals[begin] - vals[end]) / FT(float(dt)) * FT(0.1float(dt))
                TimeVaryingInputs.evaluate!(
                    dest,
                    input_periodic_calendar_linear,
                    time_to_input(times[end] + 0.1dt),
                )
                eltype(time) <: Number &&
                    @test Array(parent(dest))[1] ≈ expected
            end
        end

        # test errors when using datetime with TVI of floats or ITime without epoch
        for dest in (point_field, column_field),
            times in (times_ft, itime_no_epoch)

            input = TimeVaryingInputs.TimeVaryingInput(
                times,
                vals;
                method = TimeVaryingInputs.NearestNeighbor(),
            )
            @test_throws ErrorException TimeVaryingInputs.evaluate!(
                dest,
                input,
                start_date,
            )
        end
    end
end
