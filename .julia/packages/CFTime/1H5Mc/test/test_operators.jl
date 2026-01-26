using CFTime
using Test
using Dates

# check adding durations

dt = DateTimeNoLeap(1959, 12, 31, 23, 39, 59, 123)
@test dt + Dates.Millisecond(7) == DateTimeNoLeap(1959, 12, 31, 23, 39, 59, 130)
@test dt + Dates.Second(7) == DateTimeNoLeap(1959, 12, 31, 23, 40, 6, 123)
@test dt + Dates.Minute(7) == DateTimeNoLeap(1959, 12, 31, 23, 46, 59, 123)
@test dt + Dates.Hour(7) == DateTimeNoLeap(1960, 1, 1, 6, 39, 59, 123)
@test dt + Dates.Day(7) == DateTimeNoLeap(1960, 1, 7, 23, 39, 59, 123)
@test dt + Dates.Month(7) == DateTimeNoLeap(1960, 7, 31, 23, 39, 59, 123)
@test dt + Dates.Year(7) == DateTimeNoLeap(1966, 12, 31, 23, 39, 59, 123)
@test dt + Dates.Month(24) == DateTimeNoLeap(1961, 12, 31, 23, 39, 59, 123)

@test dt - Dates.Month(0) == DateTimeNoLeap(1959, 12, 31, 23, 39, 59, 123)
@test dt - Dates.Month(24) == DateTimeNoLeap(1957, 12, 31, 23, 39, 59, 123)
@test dt - Dates.Year(7) == DateTimeNoLeap(1952, 12, 31, 23, 39, 59, 123)


# difference between time instances
@test DateTimeStandard(2000, 1, 2) - DateTimeStandard(2000, 1, 1) == Dates.Day(1)

@test DateTimeStandard(2000, 1, 2) - DateTimeStandard(2000, 1, 1, units = :day) == Dates.Day(1)

@test (
    DateTimeStandard(2000, 1, 2) -
        DateTimeStandard(2000, 1, 1, origin = (1970, 1, 1))
) == Dates.Day(1)

@test CFTime.units(DateTimeStandard(2000, 1, 1, units = :day, origin = (1970, 1, 1))) == "days since 1970-01-01"

# check ordering

@test DateTimeStandard(2000, 1, 1) < DateTimeStandard(2000, 1, 2)
@test DateTimeStandard(2000, 1, 1) ≤ DateTimeStandard(2000, 1, 1)

@test DateTimeStandard(2000, 1, 03) > DateTimeStandard(2000, 1, 2)
@test DateTimeStandard(2000, 1, 03) ≥ DateTimeStandard(2000, 1, 1)

# issue #55

dt0 = DateTimeProlepticGregorian(2015, 1, 1)
dt1 = DateTimeProlepticGregorian(2015, 1, 2)
D = dt1 - dt0

@test D / Dates.Second(1) == 24 * 60 * 60
@test D ÷ Dates.Second(1) == 24 * 60 * 60
@test D ÷ CFTime.Picosecond(1) == 24 * 60 * 60 * 10^12
@test D / CFTime.Picosecond(1) ≈ 24 * 60 * 60 * 10^12
@test D / (24 * 60 * 60) == Dates.Second(1)

# issue #58

@test CFTime.Picosecond(1) * 2 == CFTime.Picosecond(2)
@test 2 * CFTime.Picosecond(1) == CFTime.Picosecond(2)
