using CFTime
import CFTime:
    AbstractCFDateTime,
    DateTimeStandard,
    datenum,
    datetuple,
    datetuple_ymd,
    parseDT,
    timetuplefrac,
    timeunits
using Dates
using Test

function same_tuple(t1, t2)
    len = min(length(t1), length(t2))
    return (t1[1:len] == t2[1:len]) &&
        all(==(0), t1[(len + 1):end]) &&
        all(==(0), t2[(len + 1):end])
end

@test timetuplefrac(CFTime.Period((2 * 24 * 60 * 60 + 3 * 60 * 60 + 4 * 60 + 5) * 1000, 1))[1:4] == (2, 3, 4, 5)
@test timetuplefrac(CFTime.Period((2 * 24 * 60 * 60 + 3 * 60 * 60 + 4 * 60 + 5), 1000))[1:4] == (2, 3, 4, 5)

@testset "Base methods on CFTime.Period" begin
    p = CFTime.Period(10, 1)
    @test abs(-p) == p
    @test zero(p) == CFTime.Period(0, 1)
    @test one(p) == CFTime.Period(1, 1)
end

tuf = (2, 3, 4, 5, 6, 7, 8) # day, hour, minute, seconds, ...
factor = 1
exponent = -9
p = CFTime.Period(tuf, factor, exponent)
@test timetuplefrac(p)[1:length(tuf)] == tuf

@test CFTime._type(p) == CFTime._type(typeof(p))

@test one(p) == CFTime.Period(1, factor, exponent)
# test promotion rules

p1 = CFTime.Period(1, :second)
p2 = CFTime.Period(1000, :millisecond)

@test promote_type(typeof(p1), typeof(p2)) == typeof(p2)

pp1, pp2 = promote(p1, p2)
@test pp1 === pp2


pp1, pp2 = promote(p1, Dates.Day(1))
@test pp1 == p1
@test pp2 == Dates.Day(1)
@test typeof(pp2) == typeof(pp1)


pp1, pp2 = promote(p1, Dates.Hour(1))
@test pp1 == p1
@test pp2 == Dates.Hour(1)
@test typeof(pp2) == typeof(pp1)

# missing

p = CFTime.Period(1, :second)
@test ismissing(p == missing)
@test ismissing(missing == p)
@test ismissing(DateTimeStandard(2000, 1, 1) + missing)
@test ismissing(DateTimeStandard(2000, 1, 1) - missing)
@test ismissing(missing + DateTimeStandard(2000, 1, 1))

# arithmetic

p1 = CFTime.Period(1, :microsecond)
p2 = CFTime.Period(10, :microsecond)
@test p1 + p2 == CFTime.Period(11, :microsecond)


p1 = CFTime.Period(1, :microsecond)
p2 = CFTime.Period(10, :nanosecond)
@test p1 + p2 == CFTime.Period(1010, :nanosecond)


Δt = DateTimeStandard(Int64, 2000, 1, 2) - DateTimeStandard(Int64, 2000, 1, 1)
@test Dates.value(typemax(typeof(Δt))) == typemax(Int64)
@test Dates.value(typemax(Δt)) == typemax(Int64)

@test Dates.value(typemin(typeof(Δt))) == typemin(Int64)
@test Dates.value(typemin(Δt)) == typemin(Int64)

# show and string
p1 = CFTime.Period(1, :second)
io = IOBuffer()
show(io, p1)
str = String(take!(io))
@test occursin("1 second", str)
@test occursin("1 second", string(p1))
@test occursin("second", CFTime.units(p1))
@test Dates.Second(p1) === Dates.Second(1)
@test Dates.Millisecond(p1) === Dates.Millisecond(1000)

# conversion
@test Dates.Hour(CFTime.Period(60, :minute)) == Dates.Hour(1)
@test_throws InexactError Dates.Hour(CFTime.Period(61, :minute))


# unknown CFTime.Period
p1 = CFTime.Period(1000, 1, -21) # 1000 zeptosecond = 1 attosecond
@test p1 == CFTime.Period(1, :attosecond)
@test occursin("-21", string(p1))
@test occursin("-21", CFTime.units(p1))

# decode dates

dt = DateTimeStandard(1000, "milliseconds since 2000-01-01");
@test same_tuple((2000, 1, 1, 0, 0, 1), datetuple(dt))

dt = DateTimeStandard(1, "seconds since 2000-01-01")
@test same_tuple((2000, 1, 1, 0, 0, 1), datetuple(dt))

dt = DateTimeStandard(1, "seconds since 2000-01-01")
@test same_tuple((2000, 1, 1, 0, 0, 1), datetuple(dt))

dt = DateTimeStandard(10^9, "nanoseconds since 2000-01-01");
@test same_tuple((2000, 1, 1, 0, 0, 1), datetuple(dt))

dt = DateTimeStandard(10^9, "nanoseconds since 2000-01-01T23:59:59")
@test same_tuple((2000, 1, 2), datetuple(dt))
@test Dates.day(dt) == 2
@test Dates.second(dt) == 0
@test Dates.millisecond(dt) == 0
@test Dates.microsecond(dt) == 0

dt = DateTimeStandard(1, "microseconds since 2000-01-01")
@test same_tuple((2000, 1, 1, 0, 0, 0, 0, 1), datetuple(dt))

# arithmetic with datetime and periods
dt = DateTimeStandard(1, "microseconds since 2000-01-01")
@test Dates.microsecond(dt + Dates.Microsecond(1)) == 2

@test Dates.nanosecond(dt) == 0

@test Dates.nanosecond(dt + Dates.Nanosecond(1)) == 1
@test Dates.nanosecond(dt + Dates.Nanosecond(1000)) == 0

dt = DateTimeStandard(0, "microseconds since 2000-01-01")
@test Dates.microsecond(dt + Dates.Nanosecond(1000)) == 1

dt = DateTimeStandard(1, "milliseconds since 2000-01-01T23:59:59.999")
@test same_tuple((2000, 1, 2), datetuple(dt))

dt = DateTimeStandard(1, "microseconds since 2000-01-01T23:59:59.999999")
@test same_tuple((2000, 1, 2), datetuple(dt))

dt = DateTimeStandard(1, "microseconds since 2000-01-01T23:59:59.999999")
@test same_tuple((2000, 1, 2), datetuple(dt))

dt = DateTimeStandard(1, "nanoseconds since 2000-01-01T23:59:59.999999999")
@test same_tuple((2000, 1, 2), datetuple(dt))

dt = DateTimeStandard(2001, 1, 1)
@test same_tuple((2001, 1, 1), datetuple(dt))

dt = DateTimeStandard(2001, 1, 1, 1, 2, 3, 100, 200, 300, units = :nanosecond)
@test same_tuple((2001, 1, 1, 1, 2, 3, 100, 200, 300), datetuple(dt))

dt = DateTimeStandard(0, "microseconds since 2000-01-01T23:59:59.999999")
@test string(dt) == "2000-01-01T23:59:59.999999"

# loss of precission
dt = DateTimeStandard(Float32(24 * 60 * 60 * 1000), "milliseconds since 2000-01-01")
@test Dates.hour(dt) < 24
@test Dates.minute(dt) < 60
@test Dates.second(dt) < 60

dt = DateTimeStandard(Float64(24 * 60 * 60 * 1000), "milliseconds since 2000-01-01")
@test same_tuple(datetuple(dt), (2000, 1, 2))


dt = DateTimeStandard(0, "milliseconds since 2000-01-01")
dt2 = dt + Dates.Millisecond(10) + Dates.Microsecond(20) + Dates.Nanosecond(30);
@test same_tuple(datetuple(dt2), (2000, 1, 1, 0, 0, 0, 10, 20, 30))


dt1 = DateTimeStandard(0, "microseconds since 2000-01-01")
dt2 = DateTimeStandard(10, "microseconds since 2000-01-01")
dt3 = DateTimeStandard(2, "days since 2000-01-01")
@test (dt2 - dt1) == Dates.Microsecond(10)
@test (dt2 - dt1) == Dates.Nanosecond(10_000)

@test dt1 < dt2
@test dt1 < dt3
@test dt3 > dt1

@test !(dt1 == dt2)

# ranges
dr = dt1:Dates.Microsecond(2):dt2;

@test dr[1] == dt1
@test dr[2] - dr[1] == Dates.Microsecond(2)
@test length(dr) == 6

@test_throws ArgumentError dt1:Dates.Microsecond(0):dt2
@test_throws ArgumentError Dates.len(dt1, dt2, CFTime.Period(0, :nanosecond))


#(dt2 - dt1) % Dates.Microsecond(2)
Delta = convert(CFTime.Period, Dates.Nanosecond(10_000))
@test Delta == CFTime.Period(10_000, :nanosecond)


@test convert(DateTime, DateTimeStandard(2001, 2, 3)) == DateTime(2001, 2, 3)
@test convert(DateTime, DateTimeProlepticGregorian(2001, 2, 3)) == DateTime(2001, 2, 3)
@test convert(DateTime, DateTimeStandard(2024)) == DateTime(2024)
@test convert(DateTime, DateTimeStandard(2024, 2)) == DateTime(2024, 2)

# https://en.wikipedia.org/w/index.php?title=Conversion_between_Julian_and_Gregorian_calendars&oldid=1194423852

@test convert(DateTime, DateTimeJulian(500, 2, 28)) == DateTime(500, 3, 1)
@test convert(DateTime, DateTimeJulian(1900, 3, 1)) == DateTime(1900, 3, 14)
@test convert(DateTime, DateTimeJulian(2024, 2, 13)) == DateTime(2024, 2, 26)


@test DateTimeJulian(500, 2, 28) == convert(DateTimeJulian, DateTimeProlepticGregorian(500, 3, 1))

@test DateTimeJulian(500, 2, 28) == convert(DateTimeJulian, DateTime(500, 3, 1))
@test DateTimeJulian(1900, 3, 1) == convert(DateTimeJulian, DateTimeProlepticGregorian(1900, 3, 14))
@test DateTimeJulian(2024, 2, 13) == convert(DateTimeJulian, DateTimeProlepticGregorian(2024, 2, 26))


@test DateTimeJulian(2024, 2, 13) == convert(DateTimeJulian, DateTime(2024, 2, 26))


@test DateTimeJulian(2024, 2, 13) == DateTimeProlepticGregorian(2024, 2, 26)


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


dt1 = DateTimeStandard(1, "day since 2000-01-01")
dt2 = DateTimeStandard(24, "hours since 2001-01-01")
T = typeof(dt1)
@test dt2 == convert(T, dt2)

@test daysinmonth(DateTimeAllLeap, 2001, 2) == 29
@test daysinmonth(DateTimeStandard, 1582, 10) == 21


Delta = Dates.Year(1) + Dates.Day(1)
@test DateTimeStandard(2000, 1, 1) + Delta == DateTimeStandard(2001, 1, 2)
@test DateTimeStandard(2001, 1, 2) - Delta == DateTimeStandard(2000, 1, 1)

@test daysinmonth(DateTimeStandard(1582, 10, 1)) == 21

@test parse(DateTimeNoLeap, "1999-12-05", dateformat"yyyy-mm-dd") == DateTimeNoLeap(1999, 12, 05)

dt = DateTimeStandard(1, "nanoseconds since 1999-12-31T23:59:59.999999999")
@test same_tuple((2000, 1, 1), datetuple(dt))
@test CFTime.year(dt) == 2000
@test CFTime.month(dt) == 1
@test CFTime.day(dt) == 1

@test CFTime.Year(dt) == Dates.Year(2000)
@test CFTime.Day(dt) == Dates.Day(1)

tt = (
    2000, # year
    1,   # month
    2,   # day
    3,   # hour
    4,   # minute
    5,   # second
    6,   # millisecond
    7,   # microsecond
    8,   # nanosecond
    9,   # picosecond
    10,   # femtosecond
    11,   # attosecond
    12,   # zeptosecond
    13,   # yoctosecond
    14,   # rontosecond
    15,   # quectosecond
)


dt = DateTimeStandard(Int128, tt[1:12]...)
@test CFTime.attosecond(dt) == 11

if :quectosecond in CFTime.TIME_NAMES
    dt = DateTimeStandard(Int128, tt[1:15]...)
    @test CFTime.rontosecond(dt) == 14
    @test CFTime.quectosecond(dt) == 15


    dt = DateTimeStandard(BigInt, tt...)
    @test CFTime.rontosecond(dt) == 14
    @test CFTime.quectosecond(dt) == 15

    CFTime.nanosecond(dt) == 8

    dt = DateTimeStandard(Int64, tt[1:9]...)

    @show datetuple(dt)


    tt = (
        1900, # year
        1,  # month
        2,   # day
        0,   # hour
        0,   # minute
        0,   # second
        0,   # millisecond
        0,   # microsecond
        0,   # nanosecond
        0,   # picosecond
        0,   # femtosecond
        0,   # attosecond
        0,   # zeptosecond
        0,   # yoctosecond
        0,   # rontosecond
        0,   # quectosecond
    )

    dt = DateTimeStandard(Int128, tt[1:10]...)
end


@test CFTime.datetuple(CFTime.timedecode(0, "days since -4713-01-01T12:00:00", "julian", prefer_datetime = false)) ==
    (-4713, 1, 1, 12, 0, 0, 0)

dt = CFTime.timedecode([0, 1], "years since 2000-01-01T00:00:00", prefer_datetime = false)
@test Dates.value(Dates.Millisecond(dt[2] - dt[1])) == CFTime.SOLAR_YEAR

dt = CFTime.timedecode([0, 1], "months since 2000-01-01T00:00:00", prefer_datetime = false)
@test Dates.value(Dates.Millisecond(dt[2] - dt[1])) == CFTime.SOLAR_YEAR ÷ 12

dt = CFTime.timedecode(1, "days since 2000-01-01T00:00:00", prefer_datetime = false)
@test Dates.value(dt) == dt.instant.duration


# issue 52

t1 = DateTimeStandard(Float64, 2001, 1, 1, 12; units = :day)
t2 = DateTimeStandard(Float64, 2001, 1, 1, 12)
@test Dates.value(t1 - t2) ≈ 0

# issue #47

@test Dates.Second(1) == CFTime.Picosecond(10^12)
@test CFTime.Picosecond(1) == CFTime.Femtosecond(10^3)
@test CFTime.Picosecond(1) == CFTime.Attosecond(10^6)
@test Dates.Second(1) == CFTime.Attosecond(10^18)
@test Dates.Day(1) == CFTime.Attosecond(24 * 60 * 60 * BigInt(10)^18)

dt = CFTime.timedecode(Int128[0, 1], "picoseconds since 2000-01-01T00:00:00", prefer_datetime = false)
@test dt[2] - dt[1] == CFTime.Picosecond(1)


r = DateTimeStandard(2000, 1, 1):CFTime.Picosecond(1):DateTimeStandard(2000, 1, 2);
@test length(r) == 24 * 60 * 60 * 10^12 + 1

r = DateTimeStandard(Int128, 2000, 1, 1):CFTime.Picosecond(1):DateTimeStandard(Int128, 2000, 12, 31);
@test length(r) == 365 * 24 * 60 * 60 * Int128(10)^12 + 1
