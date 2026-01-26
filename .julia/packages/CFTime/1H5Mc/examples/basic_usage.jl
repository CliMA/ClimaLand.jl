# # Basic Usage
#
# This basic example shows how to create a DateTime structure following the
# calendar from the CF conventions.


using CFTime
using Dates

# A DateTime structure following the standard calendar from the CF conventions
# can be instantiated by providing the year, month, day, hour, minute, and second:

dt0 = DateTimeStandard(1999, 12, 31, 23, 59, 59)

# A DateTimeStandard structure can also be created from a string by providing
# the format. See `Dates.DateFormat` for details about format specifiers.

dt0 = DateTimeStandard("1999-12-31T23:59:59", dateformat"yyyy-mm-ddTHH:MM:SS")

# A DateTime structure can be converted to a string explicitly using the
# `string` function, which follows the ISO 8601 format:

string(dt0)

# The conversion can also be done implicitly using string interpolation:

"$dt0"

# Different string formats are supported using the `Dates.format` function:

Dates.format(dt0, "yyyymmdd-HHMMSS")

# Year, month, day, etc. can be extracted from the DateTime structure using
# the corresponding functions:

year = Dates.year(dt0);
month = Dates.month(dt0);
day = Dates.day(dt0);
hour = Dates.hour(dt0);
minute = Dates.minute(dt0);
second = Dates.second(dt0);


# Basic arithmetic operations are supported,
# such as adding a duration to a DateTime:

dt1 = DateTimeStandard(2000, 1, 1)
dt0 + Dates.Second(1) == dt1

# Or subtracting DateTimes to obtain a duration:

dt1 - dt0

# These durations can be converted to periods of the Dates module:

Dates.Second(dt1 - dt0)

# Higher-resolution time units are also supported (the smallest resolution
# is attoseconds). To avoid overflows, one can use `Int128` or `BigInt`
# as storage types.

y, m, d = (2000, 1, 1)
hour, minute, sec = (0, 0, 0)
msec, µsec, nsec = (0, 0, 1)

dt = DateTimeStandard(
    Int128, y, m, d, hour, minute, sec, msec, µsec, nsec;
    units = :nanosecond
)

# In the example above, `units = :nanosecond` can also be omitted as a nanosecond
# argument (`nsec`) is provided.
