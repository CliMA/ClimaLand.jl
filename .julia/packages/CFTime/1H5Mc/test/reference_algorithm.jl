module Reference

import CFTime
import CFTime:
    DateTimeJulian,
    DateTimeProlepticGregorian,
    DateTimeStandard,
    DateTimeNoLeap,
    DateTimeAllLeap,
    DateTime360Day,
    _hasyear0

@inline isleap(::Type{DateTimeAllLeap}, year, has_year_zero) = true
@inline isleap(::Type{DateTimeNoLeap}, year, has_year_zero) = false
@inline isleap(::Type{DateTime360Day}, year, has_year_zero) = false

@inline function isleap(::Type{DateTimeProlepticGregorian}, year, has_year_zero)
    if (year < 0) && !has_year_zero
        year = year + 1
    end
    return (year % 400 == 0) || ((year % 4 == 0) && (year % 100 !== 0))
end

@inline function isleap(::Type{DateTimeJulian}, year, has_year_zero)
    if (year < 0) && !has_year_zero
        year = year + 1
    end
    return year % 4 == 0
end

@inline function isleap(::Type{DateTimeStandard}, year, has_year_zero)
    return if year < 1582
        isleap(DateTimeJulian, year, has_year_zero)
    else
        isleap(DateTimeProlepticGregorian, year, has_year_zero)
    end
end

@inline function month_lengths(::Type{T}, year::Integer, has_year_zero) where {T}
    if isleap(T, year, has_year_zero)
        return (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    else
        return (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    end
end

@inline function month_lengths(::Type{DateTime360Day}, year::Integer, has_year_zero)
    return ntuple(i -> 30, Val(12))
end

# Adapted
# from https://github.com/Unidata/cftime/blob/dc75368cd02bbcd1352dbecfef10404a58683f94/src/cftime/_cftime.pyx
# Licence MIT
# by Jeff Whitaker (https://github.com/jswhit)

@inline function add_timedelta_ymd(::Type{T}, (year, month, day), delta_days, julian_gregorian_mixed, has_year_zero) where {T}

    month_length = month_lengths(T, year, has_year_zero)

    n_invalid_dates =
    if julian_gregorian_mixed
        10
    else
        0
    end

    # The Julian calendar day Thursday, 4 October 1582 was
    # followed by the first day of the Gregorian calendar,
    # Friday, 15 October 1582

    # counting down
    @inbounds while delta_days < 0
        if (year == 1582) && (month == 10) && (day > 14) && (day + delta_days < 15)
            delta_days -= n_invalid_dates    # skip over invalid dates
        end

        if day + delta_days < 1
            delta_days += day
            # decrement month
            month -= 1
            if month < 1
                month = 12
                year -= 1
                if (year == 0) && !has_year_zero
                    year = -1
                end
                month_length = month_lengths(T, year, has_year_zero)
            end

            day = month_length[month]
        else
            day += delta_days
            delta_days = 0
        end
    end

    @inbounds while delta_days > 0
        if (year == 1582) && (month == 10) && (day < 5) && (day + delta_days > 4)
            delta_days += n_invalid_dates    # skip over invalid dates
        end

        if day + delta_days > month_length[month]
            delta_days -= month_length[month] - (day - 1)
            # increment month
            month += 1
            if month > 12
                month = 1
                year += 1
                if (year == 0) && !has_year_zero
                    year = 1
                end
                month_length = month_lengths(T, year, has_year_zero)
            end
            day = 1
        else
            day += delta_days
            delta_days = 0
        end
    end

    return year, month, day
end

function datetuple_ymd(::Type{T}, Z) where {T}
    has_year_zero = _hasyear0(T)
    julian_gregorian_mixed = T <: DateTimeStandard
    return datetuple_ymd(T, Z, julian_gregorian_mixed, has_year_zero)
end


@inline function datetuple_ymd(::Type{T}, delta_days, julian_gregorian_mixed, has_year_zero) where {T}
    # use the same origin
    year, month, day = CFTime.datetuple_ymd(T, 0)
    return add_timedelta_ymd(T, (year, month, day), delta_days, julian_gregorian_mixed, has_year_zero)
end

end
