using CFTime
import Dates
using Dates: DateTime, Day, @dateformat_str
using Test
import CFTime: datetuple
datetuple(dt::DateTime) = (
    Dates.year(dt), Dates.month(dt), Dates.day(dt),
    Dates.hour(dt), Dates.minute(dt), Dates.second(dt),
    Dates.millisecond(dt),
)

@testset "Check with reference values" begin
    # reference value from Meeus, Jean (1998)
    # launch of Sputnik 1
    @test CFTime.datetuple_ymd(DateTimeStandard, 2_436_116 - 2_400_001) == (1957, 10, 4)
    @test CFTime.datenum_gregjulian(1957, 10, 4, true) == 36115
    @test CFTime.datenum_gregjulian(333, 1, 27, false) == -557288


    # reference values from python's cftime
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(1) == DateTimeJulian(2000, 01, 02)
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(12) == DateTimeJulian(2000, 01, 13)
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(123) == DateTimeJulian(2000, 05, 03)
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(1234) == DateTimeJulian(2003, 05, 19)
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(12345) == DateTimeJulian(2033, 10, 19)
    @test DateTimeJulian(2000, 1, 1) + Dates.Day(12346) == DateTimeJulian(2033, 10, 20)
    @test DateTimeJulian(1, 1, 1) + Dates.Day(1234678) == DateTimeJulian(3381, 05, 14)
end

function datenum_datetuple_all_calendars(::Type{T}) where {T}
    #dayincrement = 11
    dayincrement = 11000

    for Z in (-2_400_000 + CFTime.DATENUM_OFFSET):dayincrement:(600_000 + CFTime.DATENUM_OFFSET)
        y, m, d = CFTime.datetuple_ymd(T, Z)
        Z2 = CFTime.datenum(T, y, m, d)
        if Z2 !== Z
            @show Z2, (y, m, d), Z
        end
        @test Z2 == Z
    end
    return
end

@testset "Internal consistency" begin
    for T in [
            DateTimeStandard, DateTimeJulian, DateTimeProlepticGregorian,
            DateTimeAllLeap, DateTimeNoLeap, DateTime360Day,
        ]
        datenum_datetuple_all_calendars(T)
    end
end

@testset "Leap day" begin
    # leap day
    @test DateTimeAllLeap(2001, 2, 28) + Dates.Day(1) == DateTimeAllLeap(2001, 2, 29)
    @test DateTimeNoLeap(2001, 2, 28) + Dates.Day(1) == DateTimeNoLeap(2001, 3, 1)
    @test DateTimeJulian(2001, 2, 28) + Dates.Day(1) == DateTimeJulian(2001, 3, 1)
    @test DateTimeJulian(1900, 2, 28) + Dates.Day(1) == DateTimeJulian(1900, 2, 29)
    @test DateTime360Day(2001, 2, 28) + Dates.Day(1) == DateTime360Day(2001, 2, 29)
    @test DateTime360Day(2001, 2, 29) + Dates.Day(1) == DateTime360Day(2001, 2, 30)


    @test DateTimeAllLeap(2001, 2, 29) - DateTimeAllLeap(2001, 2, 28) == Dates.Day(1)
    @test DateTimeNoLeap(2001, 3, 1) - DateTimeNoLeap(2001, 2, 28) == Dates.Day(1)
    @test DateTimeJulian(2001, 3, 1) - DateTimeJulian(2001, 2, 28) == Dates.Day(1)
    @test DateTimeJulian(1900, 2, 29) - DateTimeJulian(1900, 2, 28) == Dates.Day(1)
    @test DateTime360Day(2001, 2, 29) - DateTime360Day(2001, 2, 28) == Dates.Day(1)
    @test DateTime360Day(2001, 2, 30) - DateTime360Day(2001, 2, 29) == Dates.Day(1)
end

# generic tests
function stresstest_DateTime(::Type{DT}) where {DT}
    t0 = DT(1, 1, 1)
    for n in -800000:11:800000
        #@show n
        t = t0 + Dates.Day(n)
        y, m, d, h, mi, s, ms = CFTime.datetuple(t)
        @test DT(y, m, d, h, mi, s, ms) == t
    end
    return
end

@testset "Stresstest" begin
    for DT in [
            DateTimeStandard,
            DateTimeJulian,
            DateTimeProlepticGregorian,
            DateTimeAllLeap,
            DateTimeNoLeap,
            DateTime360Day,
        ]

        dtime = DT(1959, 12, 30, 23, 39, 59, 123)
        @test Dates.year(dtime) == 1959
        @test Dates.month(dtime) == 12
        @test Dates.day(dtime) == 30
        @test Dates.hour(dtime) == 23
        @test Dates.minute(dtime) == 39
        @test Dates.second(dtime) == 59
        @test Dates.millisecond(dtime) == 123

        @test string(DT(2001, 2, 20)) == "2001-02-20T00:00:00"
        @test CFTime.datetuple(DT(1959, 12, 30, 23, 39, 59, 123)) == (1959, 12, 30, 23, 39, 59, 123)

        stresstest_DateTime(DT)
    end
end

@testset "Error handling" begin
    @test_throws ErrorException DateTime360Day(2010, 0, 1)
    @test_throws ErrorException DateTime360Day(2010, 1, 0)
end


# test show
@testset "show" begin
    io = IOBuffer()
    show(io, DateTimeJulian(-1000, 1, 1))
    @test isempty(findfirst("Julian", String(take!(io)))) == false
end


@testset "Conversions" begin
    dt = CFTime.reinterpret(DateTimeStandard, DateTimeJulian(1900, 2, 28))
    @test typeof(dt) <: DateTimeStandard
    @test CFTime.datetuple(dt) == (1900, 2, 28, 0, 0, 0, 0)

    dt = CFTime.reinterpret(DateTime, DateTimeNoLeap(1900, 2, 28))
    @test typeof(dt) == DateTime
    @test Dates.year(dt) == 1900
    @test Dates.month(dt) == 2
    @test Dates.day(dt) == 28

    dt = CFTime.reinterpret(DateTimeNoLeap, DateTime(1900, 2, 28))
    @test typeof(dt) <: DateTimeNoLeap
    @test Dates.year(dt) == 1900
    @test Dates.month(dt) == 2
    @test Dates.day(dt) == 28

    # check conversion

    for T1 in [DateTimeProlepticGregorian, DateTimeStandard, DateTime]
        for T2 in [DateTimeProlepticGregorian, DateTimeStandard, DateTime]
            local dt1, dt2
            # datetuple should not change after 1582-10-15
            # for Gregorian Calendars
            dt1 = T1(2000, 01, 03)
            dt2 = convert(T2, dt1)

            @test CFTime.datetuple(dt1) == CFTime.datetuple(dt2)
        end
    end


    for T1 in [DateTimeStandard, DateTimeJulian]
        for T2 in [DateTimeStandard, DateTimeJulian]
            local dt1, dt2
            # datetuple should not change before 1582-10-15
            # for Julian Calendars
            dt1 = T1(200, 01, 03)
            dt2 = convert(T2, dt1)

            @test CFTime.datetuple(dt1) == CFTime.datetuple(dt2)
        end
    end

    for T1 in [DateTimeProlepticGregorian, DateTimeJulian, DateTimeStandard, DateTime]
        for T2 in [DateTimeProlepticGregorian, DateTimeJulian, DateTimeStandard, DateTime]
            local dt1, dt2
            # verify that durations (even accross 1582-10-15) are maintained
            # after convert
            dt1 = [T1(2000, 01, 03), T1(-100, 2, 20)]
            dt2 = convert.(T2, dt1)
            @test dt1[2] - dt1[1] == dt2[2] - dt2[1]
        end
    end
end  # testset "Conversions"
