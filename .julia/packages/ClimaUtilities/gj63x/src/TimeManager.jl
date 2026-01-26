"""
    TimeManager

This module facilitates calendar functions and temporal interpolations
of data.
"""
module TimeManager

import Dates

export to_datetime,
    strdate_to_datetime,
    datetime_to_strdate,
    trigger_callback,
    Monthly,
    EveryTimestep

include("ITime.jl")

"""
    to_datetime(date)

Convert a `DateTime`-like object (e.g. `DateTimeNoLeap`) to a `DateTime`.
We need this since some data files we use contain
`DateTimeNoLeap` objects for dates, which can't be used for math with `DateTime`s.
The `DateTimeNoLeap` type uses the Gregorian calendar without leap years, while
the `DateTime` type uses Gregorian calendar with leap years.

For consistency, all input data files should have dates converted to `DateTime`
before being used in a simulation.

This function is similar to `reinterpret` in CFTime.jl.

# Arguments
- `date`: `DateTime`-like object to be converted to `DateTime`
"""
function to_datetime(date)
    return Dates.DateTime(
        Dates.year(date),
        Dates.month(date),
        Dates.day(date),
        Dates.hour(date),
        Dates.minute(date),
        Dates.second(date),
        Dates.millisecond(date),
    )
end

"""
    strdate_to_datetime(strdate::String)

Convert from String ("YYYYMMDD") to Date format,
required by the official AMIP input files.
"""
strdate_to_datetime(strdate::String) = Dates.DateTime(
    parse(Int, strdate[1:4]),
    parse(Int, strdate[5:6]),
    parse(Int, strdate[7:8]),
)

"""
    datetime_to_strdate(datetime::Dates.DateTime)

Convert from DateTime to String ("YYYYMMDD") format.
"""
datetime_to_strdate(datetime::Dates.DateTime) =
    string(lpad(Dates.year(datetime), 4, "0")) *
    string(string(lpad(Dates.month(datetime), 2, "0"))) *
    string(lpad(Dates.day(datetime), 2, "0"))

abstract type AbstractFrequency end

"Struct used to dispatch callback that is triggered monthly"
struct Monthly <: AbstractFrequency end

"Struct used to dispatch callback that is triggered every timestep"
struct EveryTimestep <: AbstractFrequency end

"""
    trigger_callback(date_nextcall::Dates.DateTime,
        date_current::Dates.DateTime,
        ::Monthly,
        func::Function,)

If the current date is equal to or later than the "next call" date at time
00:00:00, call the callback function and increment the next call date by one
month. Otherwise, do nothing and leave the next call date unchanged.

The tuple of arguments `func_args` must match the types, number, and order
of arguments expected by `func`.

# Arguments
- `date_nextcall::DateTime` the next date to call the callback function at or after
- `date_current::DateTime` the current date of the simulation
- `save_freq::AbstractFrequency` frequency with which to trigger callback
- `func::Function` function to be triggered if date is at or past the next call date
- `func_args::Tuple` a tuple of arguments to be passed into the callback function
"""
function trigger_callback(
    date_nextcall::Dates.DateTime,
    date_current::Dates.DateTime,
    ::Monthly,
    func::Function,
    func_args::Tuple,
)
    if date_current >= date_nextcall
        func(func_args...)
        return date_nextcall + Dates.Month(1)
    else
        return date_nextcall
    end
end

end # module TimeManager
