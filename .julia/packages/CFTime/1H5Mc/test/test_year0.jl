using CFTime
using Test
using Dates: DateTime

# handling of year zero, reference values
# from pythons cftime 1.6.0
# issue #17

@test DateTimeProlepticGregorian(1, 1, 1) - Day(1) == DateTimeProlepticGregorian(0, 12, 31)
@test DateTimeJulian(1, 1, 1) - Day(1) == DateTimeJulian(-1, 12, 31)
@test DateTimeStandard(1, 1, 1) - Day(1) == DateTimeStandard(-1, 12, 31)
@test DateTimeAllLeap(1, 1, 1) - Day(1) == DateTimeAllLeap(0, 12, 31)
@test DateTimeNoLeap(1, 1, 1) - Day(1) == DateTimeNoLeap(0, 12, 31)
@test DateTime360Day(1, 1, 1) - Day(1) == DateTime360Day(0, 12, 30)


dt = CFTime.timedecode(0, "days since 0001-01-01", prefer_datetime = false)
@test dt == DateTimeStandard(1, 1, 1)

# https://www.fourmilab.ch/documents/calendar/
dt = CFTime.timedecode(0, "days since 0001-01-01") # julian date
@test dt == DateTime(0, 12, 30) # converted to gregorian date


dt = CFTime.timedecode(0, "days since 0001-01-01", "proleptic_gregorian")
@test dt == DateTime(1, 1, 1)

#=
import netCDF4

dates = netCDF4.num2date([-1],units="days since 0001-01-01",calendar="proleptic_gregorian")
print("dates",dates)

=#
dt = CFTime.timedecode(-1, "days since 0001-01-01", "proleptic_gregorian")
@test dt == DateTime(0, 12, 31)


dt = CFTime.timedecode(0, "days since 0001-01-01", prefer_datetime = false)
@test dt == DateTimeStandard(1, 1, 1)

# no year 0
dt = CFTime.timedecode(-1, "days since 0001-01-01", prefer_datetime = false)
@test dt == DateTimeStandard(-1, 12, 31)
@test CFTime.year(dt) == -1
@test CFTime.month(dt) == 12

# test a dummy calendar with no year 0 and regular month length
struct DummyDataTime{T, Torigintuple} <: CFTime.AbstractCFDateTime{T, Torigintuple}
end
import CFTime: _cum_month_length, _hasyear0
_cum_month_length(::Type{DummyDataTime}) = (0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360)
_hasyear0(::Type{DummyDataTime}) = false

Z = datenum(DummyDataTime, -1, 1, 1)
@test datetuple_ymd(DummyDataTime, Z) == (-1, 1, 1)
