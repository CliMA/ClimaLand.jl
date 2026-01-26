# solar year in ms (the interval between 2 successive passages of the sun
# through vernal equinox)
const SOLAR_YEAR = round(Int64, 365.242198781 * 24 * 60 * 60 * 1000)

const DEFAULT_TIME_UNITS = "days since 1900-01-01 00:00:00"

# Introduction of the Gregorian Calendar 1582-10-15
const GREGORIAN_CALENDAR = (1582, 10, 15)

# Time offset in days for the time origin
# if DATENUM_OFFSET = 0, then datenum_gregjulian
# corresponds to  Modified Julian Days (MJD).
# MJD is the number of days since midnight on 1858-11-17)

#const DATENUM_OFFSET = 2_400_000.5 # for Julian Days
const DATENUM_OFFSET = 0 # for Modified Julian Days

# Introduction of the Gregorian Calendar 1582-10-15
# expressed in MJD (if DATENUM_OFFSET = 0)

const DN_GREGORIAN_CALENDAR = -100840 + DATENUM_OFFSET

# DateTime(Dates.UTInstant{Millisecond}(Dates.Millisecond(0)))
# returns 0000-12-31T00:00:00
# 678576 is the output of -datenum(DateTimeProlepticGregorian,-1,12,31)

const DATETIME_OFFSET = Dates.Millisecond(678576 * (24 * 60 * 60 * Int64(1000)))

# all supported time units, e.g.
# 1 day =  24*60*60 × 10⁰ s
# 1 nanosecond = 1 × 10⁻⁹ s

# Note: the time of first call for datetuple does increases significantly
# if we use zeptoseconds and beyond.
# Currently we support only resolutions down to attosecond resolution (as
# does numpy).

const TIME_DIVISION = (
    # name           factor, exponent
    (:day, 24 * 60 * 60, 0),
    (:hour, 60 * 60, 0),
    (:minute, 60, 0),
    (:second, 1, 0),
    (:millisecond, 1, -3),
    (:microsecond, 1, -6), # time of first call for datetuple
    (:nanosecond, 1, -9), # 1.27 seconds
    (:picosecond, 1, -12), # 1.52 seconds
    (:femtosecond, 1, -15), # 2.07 seconds
    (:attosecond, 1, -18), # 3.11 seconds
    #    (:zeptosecond,        1,    -21), # 5.35 seconds
    #    (:yoctosecond,        1,    -24), # 10.36 seconds
    #    (:rontosecond,        1,    -27), # 23.88 seconds
    #    (:quectosecond,       1,    -30), # 60.37 seconds
)

const TIME_NAMES = (:year, :month, first.(TIME_DIVISION)...)
