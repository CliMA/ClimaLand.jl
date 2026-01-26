# query functions

# https://web.archive.org/web/20240304171516/https://unidata.github.io/cftime/api.html
# https://github.com/cf-convention/cf-conventions/issues/298
@inline _hasyear0(::Type{T}) where {T <: Union{DateTimeJulian, DateTimeStandard}} = false
@inline _hasyear0(::Type{T}) where {T <: AbstractCFDateTime} = true


"""
    monthlength = daysinmonth(::Type{DT},y,m)

Returns the number of days in a month for the year `y` and the month `m`
according to the calendar given by the type `DT` (any subtype of [`AbstractCFDateTime`](@ref)).

Example
```julia-repl
julia> daysinmonth(DateTimeAllLeap,2001,2)
29
```

"""
function daysinmonth(::Type{DT}, y, m) where {DT <: AbstractCFDateTime}
    t = DT(y, m, 1)
    return Dates.value((t + Dates.Month(1)) - t) ÷ (24 * 60 * 60 * 1000)
end

"""
    monthlength = daysinmonth(t)

Returns the number of days in a month containing the date `t`

Example
```julia-repl
julia> daysinmonth(DateTimeAllLeap(2001,2,1))
29
```
"""
function daysinmonth(t::DT) where {DT <: AbstractCFDateTime}
    return daysinmonth(DT, Dates.year(t), Dates.month(t))
end

"""
    yearlength = daysinyear(::Type{DT},y)

Returns the number of days in a year for the year `y`
according to the calendar given by the type `DT` (any subtype of [`AbstractCFDateTime`](@ref)).

Example
```julia-repl
julia> daysinyear(DateTimeAllLeap,2001,2)
366
```

"""
function daysinyear(::Type{DT}, y) where {DT <: AbstractCFDateTime}
    t = DT(y, 1, 1)
    return Dates.value((t + Dates.Year(1)) - t) ÷ (24 * 60 * 60 * 1000)
end

"""
    yearlength = daysinyear(t)

Returns the number of days in a year containing the date `t`

Example
```julia-repl
julia> daysinyear(DateTimeAllLeap(2001,2,1))
366
```
"""
function daysinyear(t::DT) where {DT <: AbstractCFDateTime}
    return daysinyear(DT, Dates.year(t))
end

"""
    yearmonthday(dt::AbstractCFDateTime)

Simultaneously return the year, month and day parts of `dt`.
"""
yearmonthday(dt::AbstractCFDateTime) = _datetuple(dt)[1:3]

"""
    yearmonth(dt::AbstractCFDateTime)

Simultaneously return the year and month parts of `dt`.
"""
yearmonth(dt::AbstractCFDateTime) = _datetuple(dt)[1:2]

"""
    monthday(dt::AbstractCFDateTime)

Simultaneously return the month and day parts of `dt`.
"""
monthday(dt::AbstractCFDateTime) = _datetuple(dt)[2:3]


"""
    firstdayofyear(dt::AbstractCFDateTime)

Return the first day of the year including the date `dt`
"""
@inline firstdayofyear(dt::T) where {T <: AbstractCFDateTime} = T(Dates.year(dt), 1, 1, 0, 0, 0)


"""
    dayofyear(dt::AbstractCFDateTime)

Return the day of the year for dt with January 1st being day 1.
"""
function dayofyear(dt::AbstractCFDateTime)
    t0 = firstdayofyear(dt)
    return Dates.value(floor(dt - t0, Dates.Day)) + 1
end
