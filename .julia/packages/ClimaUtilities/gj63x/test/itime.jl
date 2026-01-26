import ClimaUtilities
using ClimaUtilities.TimeManager

using Test, Dates

@testset "ITime" begin
    @testset "Constructors" begin
        # Constructor with just an integer counter
        t1 = ITime(10)
        @test t1.counter == 10
        @test t1.period == Dates.Second(1)
        @test isnothing(t1.epoch)

        # Constructor with just an integer counter in Int32
        t1_int32 = ITime(Int32(10))
        @test t1_int32.counter isa Int32
        @test t1_int32.counter == 10

        # Constructor with start date
        epoch = Dates.DateTime(2024, 1, 1)
        t3 = ITime(10, epoch = epoch)
        @test t3.epoch == epoch

        # Start date as Date (not DateTime)
        t4 = ITime(10, epoch = Dates.Date(2024, 1, 1))
        @test t4.epoch == epoch

        # Explicit period
        t5 = ITime(10, period = Dates.Millisecond(100))
        @test t5.period == Dates.Millisecond(100)

        # From float
        t6 = ITime(0.0)
        @test t6.counter == 0
        @test t6.period == Dates.Second(1)

        t7 = ITime(1.0)
        @test t7.counter == 1
        @test t7.period == Dates.Second(1)

        t8 = ITime(0.001)
        @test t8.counter == 1
        @test t8.period == Dates.Millisecond(1)

        t9 = ITime(1.5)
        @test t9.counter == 1500
        @test t9.period == Dates.Millisecond(1)


        t10 = ITime(1.5; epoch = Dates.DateTime(2020, 1, 1))
        @test t10.epoch == Dates.DateTime(2020, 1, 1)

        # Cannot be represented exactly
        @test_throws ErrorException ITime(1e-20)

        # Compare ITime(::Int) to ITime(::Float)
        function compare_itime_constructors(counter, t)
            a = ITime(counter)
            b = ITime(t)
            @test a.counter == b.counter
            @test a.period == b.period
            @test a.epoch == b.epoch
        end
        for nums in (0.0, 1.0, 60.0, 120.0, 3600.0)
            compare_itime_constructors(Int(nums), nums)
        end
    end

    @testset "Accessors" begin
        t = ITime(10, Dates.Millisecond(50), Dates.DateTime(2024, 1, 1))
        @test counter(t) == 10
        @test period(t) == Dates.Millisecond(50)
        @test epoch(t) == Dates.DateTime(2024, 1, 1)
    end

    @testset "date" begin
        t1 = ITime(10, epoch = Dates.DateTime(2024, 1, 1))
        @test date(t1) == Dates.DateTime(2024, 1, 1) + Dates.Second(10)
        @test Dates.DateTime(t1) ==
              Dates.DateTime(2024, 1, 1) + Dates.Second(10)

        # Cannot convert to date without a start date
        t2 = ITime(10)
        @test_throws ErrorException date(t2) # No start date
    end

    @testset "Promote" begin
        t1 = ITime(10, period = Dates.Millisecond(100))
        t2 = ITime(100, period = Dates.Millisecond(10))
        t1_promoted, t2_promoted = promote(t1, t2)
        @test t1_promoted.counter == 100
        @test t2_promoted.counter == 100
        @test t1_promoted.period == Dates.Millisecond(10)

        t3 = ITime(10, epoch = Dates.DateTime(2024, 1, 1))
        t4 = ITime(20)

        t3_promoted, t4_promoted = promote(t3, t4)
        @test t3_promoted.epoch == Dates.DateTime(2024, 1, 1)
        @test t4_promoted.epoch == Dates.DateTime(2024, 1, 1)

        @test_throws ErrorException promote(
            t3,
            ITime(10, epoch = Dates.DateTime(2024, 1, 2)),
        )
    end

    @testset "Arithmetic Operations" begin
        t1 = ITime(10)
        t2 = ITime(5)
        @test t1 + t2 == ITime(15)
        @test t1 - t2 == ITime(5)
        @test -t1 == ITime(-10)
        @test abs(t1) == ITime(10)

        t3 = ITime(10, period = Dates.Millisecond(100))
        t4 = ITime(100, period = Dates.Millisecond(10))
        @test t3 + t4 == ITime(200, period = Dates.Millisecond(10)) # 10*100 + 100*10 = 2000 ms = 2s

        @test !(t1 < t2)
        @test t1 > t2
        @test t1 == ITime(10)
        @test t1 >= ITime(10)
        @test t1 <= ITime(10)
        @test t1 != ITime(5)
        @test isapprox(t1, ITime(10))

        @test t1 * 2 == ITime(20)
        @test 2 * t1 == ITime(20)
        @test div(t1, 2) == ITime(5)
        @test t1 / t2 == 2
        @test t2 / t1 == 1 // 2

        @test one(t1) == 1
        @test oneunit(t1) == ITime(1)
        @test zero(t1) == ITime(0)


        t5 = ITime(10, epoch = Dates.DateTime(2024, 1, 1))
        t6 = ITime(5, epoch = Dates.DateTime(2024, 1, 1))
        @test (t5 + t6).epoch == t5.epoch

        t7 = ITime(5, epoch = Dates.DateTime(2024, 10, 1))
        # Arithmetic operations between ITime with different epochs are disallowed
        @test_throws ErrorException t6 + t7
    end

    @testset "Float Conversion and Broadcasting" begin
        t1 = ITime(10, period = Dates.Millisecond(100))
        @test float(t1) == 1.0
        @test Float64(t1) === Float64(1)
        @test Float32(t1) === Float32(1)

        # Test broadcasting (simple example)
        @test float.(t1) == 1.0
    end

    @testset "Show method" begin
        t1 = ITime(10)
        @test sprint(show, t1) ==
              "10.0 seconds [counter = 10, period = 1 second]"

        t2 = ITime(10, epoch = Dates.DateTime(2024, 1, 1))
        @test sprint(show, t2) ==
              "10.0 seconds (2024-01-01T00:00:10) [counter = 10, period = 1 second, epoch = 2024-01-01T00:00:00]"

        t3 = ITime(
            10,
            period = Dates.Hour(1),
            epoch = Dates.DateTime(2024, 1, 1),
        )
        @test sprint(show, t3) ==
              "10.0 hours (2024-01-01T10:00:00) [counter = 10, period = 1 hour, epoch = 2024-01-01T00:00:00]"
    end

    @testset "Find common epoch" begin
        t1 = ITime(0, period = Second(1), epoch = Date(2011))
        t2 = ITime(0, period = Minute(1), epoch = Date(2011))
        t3 = ITime(0, period = Hour(1))
        t4 = ITime(0, period = Day(1), epoch = Date(2012))

        @test ClimaUtilities.TimeManager.find_common_epoch(t1, t2) == Date(2011)
        @test ClimaUtilities.TimeManager.find_common_epoch(t1, t3) == Date(2011)
        @test_throws ErrorException ClimaUtilities.TimeManager.find_common_epoch(
            t1,
            t4,
        )
        @test_throws ErrorException ClimaUtilities.TimeManager.find_common_epoch(
            t1,
            t2,
            t3,
            t4,
        )
    end

    @testset "Range" begin
        start = ITime(0, period = Hour(1))
        step = ITime(1, period = Second(1), epoch = Date(2011))
        stop = ITime(1, period = Minute(1))
        stop1 = ITime(0, period = Day(1), epoch = Date(2012))

        start2stop = collect(start:step:stop)
        @test length(start2stop) == 61
        @test start2stop[begin] ==
              ITime(0, period = Second(1), epoch = Date(2011))
        @test start2stop[2] == ITime(1, period = Second(1), epoch = Date(2011))
        @test start2stop[end] ==
              ITime(60, period = Second(1), epoch = Date(2011))

        @test_throws ErrorException start:step:stop1
    end

    @testset "% and mod operators" begin
        t1 = ITime(0, period = Hour(1))
        t2 = ITime(1, period = Second(1))
        t3 = ITime(10, period = Day(1))
        t4 = ITime(7, period = Second(1))
        t5 = ITime(2, period = Second(1))
        t6 = ITime(9, period = Second(2))
        t7 = ITime(2, period = Second(2))
        @test t1 % t2 == ITime(0, period = Second(1))
        @test t3 % t2 == ITime(0, period = Second(1))
        @test t4 % t5 == ITime(1, period = Second(1))
        @test t6 % t7 == ITime(1, period = Second(2))

        @test mod(t1, t2) == ITime(0, period = Second(1))
        @test mod(t3, t2) == ITime(0, period = Second(1))
        @test mod(t4, t5) == ITime(1, period = Second(1))
        @test t6 % t7 == ITime(1, period = Second(2))
    end

    @testset "Base.:(==)" begin
        @test ITime(0; period = Day(1), epoch = DateTime(2010)) ==
              DateTime(2010)
        @test DateTime(2011) ==
              ITime(365; period = Day(1), epoch = DateTime(2010))
        @test ITime(365; period = Day(1), epoch = DateTime(2010)) !=
              DateTime(2010)
        @test_throws "Cannot compare ITime with a Number" ITime(
            1;
            period = Day(1),
        ) == 60 * 60 * 24
        @test_throws "Cannot compare ITime with a Number" Inf == ITime(1.0;)
    end

    @testset "iszero" begin
        t1 = ITime(0, period = Hour(1))
        t2 = ITime(1, period = Second(1))
        @test iszero(t1) == true
        @test iszero(t2) == false
    end

    @testset "length" begin
        t1 = ITime(0, period = Hour(1))
        @test length(t1) == 1
    end

    @testset "Multiply a float by an ITime" begin
        t1 = ITime(1, Minute(1), DateTime(2010))
        t2 = ITime(1, Second(1), DateTime(2010))
        @test_throws ErrorException 4.2 * t1
        @test_throws ErrorException t1 * 4.2
        @test_throws ErrorException -4.2 * t1
        @test_throws ErrorException t1 * -4.2

        t1_no_round = 1.0 * t1
        t1_round_down = 0.3 * t1
        t1_round_up = 0.6 * t1
        t2_no_round = 1.0 * t2
        t2_round_down = 0.3 * t2
        t2_round_up = 0.6 * t2

        @test t1_no_round == ITime(1, Minute(1), DateTime(2010))
        @test t1_round_down == ITime(0, Minute(1), DateTime(2010))
        @test t1_round_up == ITime(1, Minute(1), DateTime(2010))
        @test t2_no_round == ITime(1, Second(1), DateTime(2010))
        @test t2_round_down == ITime(0, Second(1), DateTime(2010))
        @test t2_round_up == ITime(1, Second(1), DateTime(2010))
    end

    @testset "Int32 test" begin
        # Check if any operation with Int32 in ITime can lead to an Int64
        t1 = ITime(Int32(0), period = Second(1))
        t2 = ITime(Int32(1), period = Minute(1))
        t3 = ITime(Int32(2), period = Hour(1), epoch = DateTime(2010))

        # Test promote
        tt1, tt2, tt3 = promote(t1, t2, t3)
        @test typeof(tt1.counter) == Int32
        @test typeof(tt2.counter) == Int32
        @test typeof(tt3.counter) == Int32

        # Test addition
        t7 = t1 + t2
        t8 = t2 + t3
        @test typeof(t7.counter) == Int32
        @test typeof(t8.counter) == Int32

        # Test subtraction
        t7 = t1 - t2
        t8 = t2 - t3
        @test typeof(t7.counter) == Int32
        @test typeof(t8.counter) == Int32

        # Test multiplication by an int
        t7 = t1 * Int32(1)
        t9 = t3 * Int32(3)
        @test typeof(t7.counter) == Int32
        @test typeof(t9.counter) == Int32

        # Test multiplication by float
        ft1 = 0.6 * t2
        @test typeof(ft1.counter) == Int32
        ft2 = 0.3 * t2
        @test typeof(ft2.counter) == Int32

        # Test oneunit, zero, and mod
        @test typeof(oneunit(t1).counter) == Int32
        @test typeof(zero(t1).counter) == Int32
        @test typeof(mod(t3, t2).counter) == Int32
    end

    @testset "Overflow test" begin
        t1 = ITime(typemax(Int64), period = Second(2))
        # Dates.seconds(t.period) returns an integer if t.period is in seconds
        # Multiplying this by t.counter, an integer, results in an integer that
        # could overflow
        @test isapprox(float(t1), float(typemax(Int64)) * 2.0)

        # The computation epoch(t) + counter(t) * period(t) could overflow since
        # counter(t) * period(t) could overflow. In this case, the maximum value
        # of period(t) is 2^63 - 1 nanoseconds. Multiplying 2 nanoseconds by the
        # counter overflow and give -2 nanoseconds. Here, we decide to throw an
        # error instead of circumventing the overflow. One way to circumvent the
        # overflow is to add counter(t) * oneunit(period(t)) n times, where n is
        # the value of period.
        t2 = ITime(
            typemax(Int64),
            period = Nanosecond(2),
            epoch = Dates.DateTime(2010),
        )
        @test_throws ErrorException date(t2)

        # Same tests as above, but with Int32
        t3 = ITime(typemax(Int32), period = Second(2))
        @test isapprox(float(t3), float(typemax(Int32)) * 2.0)

        t4 = ITime(
            typemax(Int32),
            period = Nanosecond(2),
            epoch = Dates.DateTime(2010),
        )
        # We do not overflow because Period uses Int64 internally
        @test date(t4) == Dates.DateTime(2010) + typemax(Int32) * Nanosecond(2)

        # Should not overflow
        t5 = ClimaUtilities.TimeManager.ITime(
            typemax(Int64),
            period = Nanosecond(1),
            epoch = DateTime(2010),
        )
        @test date(t5) ==
              Dates.DateTime(2010) + (typemax(Int64) * Dates.Nanosecond(1))
    end
end
