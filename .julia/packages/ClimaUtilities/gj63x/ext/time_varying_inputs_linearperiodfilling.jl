# This file is included in TimeVaryingInputsExt.jl

"""
    _period_difference(date_left, date_right, period::Dates.DatePeriod)

Return the difference in periods (in `Dates.Period`) from `date_right` to `date_left`.

For example, if `date_left` is 18/08/1995, `date_right` 17/01/1997, and `period = Year(1)`,
the difference is 2 years.

This difference can also be used move a date in `period_right` to `period_left` (or viceversa).

In the example above, 18/08/1995 + 2 years is in `period_right` (or 17/01/1997 - 2 years is
in `period_left`).
"""
function _period_difference(date_left, date_right, period::Dates.DatePeriod)
    return beginningofperiod(date_right, period) -
           beginningofperiod(date_left, period)
end

"""
    _extract_period(date, period::Dates.DatePeriod)

Return the beginning of the period with duration specified by `period`.

This function is here mostly for clarity (because extracting period clearer in intent than beginning
of period).

Example
========

# TODO This jldoctest would require loading the extension...

```
julia> _extract_period(Date(1993, 8, 18), Year(1))
DateTime(1993, 1, 1)

julia> _extract_period(Date(1993, 8, 18), Month(1))
DateTime(1993, 8, 1)
```
"""
function _extract_period(date, period::Dates.DatePeriod)
    return beginningofperiod(date, period)
end

"""
    _neighboring_periods(date, available_periods, period::Dates.DatePeriod)

Return the two periods within `available_periods` that neighbor `date`.

Example
========

# TODO This jldoctest would require loading the extension...

```
julia> available_periods = [DateTime(1985, 1, 1), DateTime(1995, 1, 1), DateTime(2005, 1, 1)];

julia> _neighboring_periods([DateTime(1993, 8, 18)], available_periods, Year(1))
(DateTime(1985, 1, 1), DateTime(1995, 1, 1))

julia> _extract_period(Date(1993, 8, 18), Month(1))
DateTime(1993, 8, 1)
```
"""
function _neighboring_periods(date, available_periods, period::Dates.DatePeriod)
    index_period_right =
        findfirst(d -> d >= beginningofperiod(date, period), available_periods)
    index_period_left = index_period_right - 1
    return available_periods[index_period_left],
    available_periods[index_period_right]
end

"""
    _move_date_to_period(date, target_period, period::Dates.DatePeriod)

Return a new date so that it is the equivalent date in the `target_period`

Example
========

# TODO This jldoctest would require loading the extension...

```
julia> _move_date_to_period(DateTime(1987, 3, 7), DateTime(1985, 1, 1), Year(1))
DateTime(1985, 3, 7)

julia> _move_date_to_period(DateTime(1987, 3, 7), DateTime(1995, 1, 1), Year(1))
DateTime(1995, 3, 7)

```
"""
function _move_date_to_period(date, target_period, period::Dates.DatePeriod)
    period_offset = _period_difference(date, target_period, period)
    return date + period_offset
end

"""
    _date_in_range(date, left, right)

Check if `date` is between `left` and `right`.

Example
========

# TODO This jldoctest would require loading the extension...

```
julia> _date_in_range(DateTime(1993, 10, 11), DateTime(1993, 1, 11), DateTime(1993, 12, 11))
true

julia> _date_in_range(DateTime(1993, 1, 1), DateTime(1993, 1, 11), DateTime(1993, 12, 11))
false
```
"""
function _date_in_range(date, left, right)
    return left <= date <= right
end

"""
    _interpolable_range(dates_period_left, dates_period_right, period)

Find the two dates that bound the interpolable region.

The dates are returned defined in `period_left`

The interpolable region is defined as the dates when interpolation can be immediately
performed because corresponding dates in both periods are available.

Have a look at the examples below for more information.

Example
========

# TODO This jldoctest would require loading the extension...

```
julia> dates_period_left = [DateTime(1985, 1, 15), DateTime(1985, 2, 17), DateTime(1985, 12, 11)];

julia> dates_period_right = [DateTime(1995, 1, 8), DateTime(199h5, 2, 18), DateTime(1995, 12, 18)];

julia> period = Year(1);

julia> _interpolable_range(dates_period_left, dates_period_right, period)
[DateTime(1985, 1, 15), DateTime(1985, 12, 11)]

julia> _interpolable_range([DateTime(1985, 1, 1), DateTime(1985, 12, 11)], dates_period_left, period)
[DateTime(1985, 1, 8), DateTime(1985, 12, 11)]

julia> _interpolable_range(dates_period_left, [DateTime(1995, 1, 3), DateTime(1985, 12, 21)], period)
[DateTime(1985, 1, 3), DateTime(1985, 12, 21)]
```
"""
function _interpolable_range(
    dates_period_left,
    dates_period_right,
    period::Dates.DatePeriod,
)
    if length(unique(beginningofperiod.(dates_period_left, period))) != 1
        error("dates in dates_period_left belongs to different periods")
    end
    if length(unique(beginningofperiod.(dates_period_right, period))) != 1
        error("dates in dates_period_right belongs to different periods")
    end

    # First, we move all the dates to period_left
    period_offset = _period_difference(
        first(dates_period_left),
        first(dates_period_right),
        period,
    )
    dates_period_right_moved_to_left = dates_period_right .- period_offset

    # Now, to identify the interpolable range, we go through all the dates and find
    # the first that is larger than both the smallest (largest) entries in
    # dates_period_left and dates_period_right_moved_to_left
    #
    # For example, for
    # dates_period_left = [Date(1985, 1, 15), Date(1985, 2, 17), Date(1985, 12, 11)]
    # dates_period_right = [Date(1995, 1, 8), Date(1995, 2, 18), Date(1995, 12, 18)]
    #
    # largest_first_date would be Date(1985, 1, 15) and smallest_last_date Date(1985, 12, 11).
    largest_first_date = maximum([
        minimum(dates_period_left),
        minimum(dates_period_right_moved_to_left),
    ])
    smallest_last_date = minimum([
        maximum(dates_period_left),
        maximum(dates_period_right_moved_to_left),
    ])

    return largest_first_date, smallest_last_date
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time,
    method::LinearPeriodFillingInterpolation,
    args...;
    kwargs...,
)
    # E.g., Year(1)
    period = method.period

    # E.g., [Date(1985, 1, 15), Date(1985, 2, 17), Date(1995, 1, 8), Date(1995, 2, 18),
    # Date(2015, 1, 10)]
    available_dates = DataHandling.available_dates(itp.data_handler)

    # E.g., [Date(1985, 1, 1), Date(1995, 1, 1), Date(2005, 1, 1)]
    available_periods = unique_periods(available_dates, period)

    # E.g., Date(1987, 2, 1)
    target_date = time

    # If the target date falls within an available period, we just interpolate it with
    # LinearInterpolation()
    if _extract_period(target_date, period) in available_periods
        # Here we are just interpolating within some available data, so it is easy
        TimeVaryingInputs.evaluate!(dest, itp, time, LinearInterpolation())
        return nothing
    end

    # We can assume that time is defined in the range of available_times because we handle
    # extrapolation boundary conditions elsewhere. For this reason, we will always have a
    # left and right periods

    # We are not in one of the available_periods, we have two cases: either we are in the
    # interpolable region, or not. The interpolable region is a region where we can directly
    # perform linear interpolation in both the left and right periods, and then interpolate
    # the resulting two values.
    #
    # When we are not in the interpolable region, we have to reduce to that case. So, we
    # have to identify two dates in the interpolable region that bound the desired date and
    # interpolate across them.

    # itp.preallocated_regridded_fields[1:2] is used internally by functions called by this
    # function, so we should not used them here. Also, we will use dest as a working area to
    # avoid extra allocations
    tmp_field1, tmp_field2 = itp.preallocated_regridded_fields[(end - 1):end]

    # E.g, Date(1985, 1, 1), Date(1995, 1, 1)
    period_left, period_right =
        _neighboring_periods(target_date, available_periods, period)

    # E.g., Date(1985, 1, 15), Date(1985, 2, 17)
    dates_period_left = filter(
        d -> period_left <= d <= endofperiod(period_left, period),
        available_dates,
    )

    # E.g., Date(1995, 1, 8), Date(1995, 2, 18)
    dates_period_right = filter(
        d -> period_right <= d <= endofperiod(period_right, period),
        available_dates,
    )

    # Move the target date to the left period, where we can compare it with the interpolable
    # range
    # E.g., Date(1985, 2, 1)
    target_date_in_left_period =
        _move_date_to_period(target_date, period_left, period)

    # E.g., [Date(1985, 1, 8), Date(1985, 2, 17)]
    min_interpolable_range, max_interpolable_range =
        _interpolable_range(dates_period_left, dates_period_right, period)

    in_interpolable_region = _date_in_range(
        target_date_in_left_period,
        min_interpolable_range,
        max_interpolable_range,
    )

    if in_interpolable_region
        date_pre = _move_date_to_period(target_date, period_left, period)
        date_post = _move_date_to_period(target_date, period_right, period)

        # Linear interpolation: y = y0 * (1 - coeff) + coeff * y1
        #
        # coeff here is the period weight, for example, if we are interpolating in 1987
        # from 1985 and 1985, it would be 2/10.

        # E.g., 730 days
        offset_periods_left =
            beginningofperiod(target_date, period) - period_left
        period_offset = period_right - period_left
        coeff = offset_periods_left / period_offset

        TimeVaryingInputs.evaluate!(dest, itp, date_pre, LinearInterpolation())
        dest .*= (1 - coeff)
        TimeVaryingInputs.evaluate!(
            tmp_field1,
            itp,
            date_post,
            LinearInterpolation(),
        )
        tmp_field1 .*= coeff
        dest .+= tmp_field1
        return nothing
    else
        # In this branch, we are not in the interpolable region. This can happen because the
        # target_date is earlier than the left boundary or later than the right boundary of
        # the interpolable region in one/both period_left and period_right.
        #
        # When that happens, we change the interpolation dates to be the boundaries of the
        # interpolable region moved to the relevant period
        #
        # For example, for Date(1987, 1, 1), we would interpolate it with
        # Date(1986, 2, 17) and Date(1987, 1, 15)
        #
        target_date_in_left_period =
            _move_date_to_period(target_date, period_left, period)
        target_date_in_right_period =
            _move_date_to_period(target_date, period_right, period)

        target_period = _extract_period(target_date, period)

        if target_date_in_left_period >= maximum(dates_period_left) ||
           target_date_in_right_period >= maximum(dates_period_right)
            # First case, the date is later than the interpolable region
            #
            # E.g., Date(1987, 12, 31), in this case we interpolate with
            # Date(1987, 2, 17) and Date(1988, 1, 15)
            #
            # For date_pre, we take the max of the interpolable range and bring it to the
            # current period. For date_post, we take the min, and take it to the next period
            date_pre = _move_date_to_period(
                max_interpolable_range,
                target_period,
                period,
            )
            date_post = _move_date_to_period(
                min_interpolable_range,
                target_period + period,
                period,
            )
        elseif minimum(dates_period_left) >= target_date_in_left_period ||
               minimum(dates_period_right) >= target_date_in_right_period
            # Second case, the date is  than the interpolable region
            #
            # E.g., Date(1987, 1, 1), in this case we interpolate with
            # Date(1986, 2, 17) and Date(1987, 1, 15)
            #
            # For date_pre, we take the max of the interpolable range and bring it to the
            # previous period. For date_post, we take the min, and take it to the current period
            date_pre = _move_date_to_period(
                max_interpolable_range,
                target_period + period,
                period,
            )
            date_post = _move_date_to_period(
                min_interpolable_range,
                target_period,
                period,
            )
        else
            error("We should not be here!")
        end

        # y = y0 * (1 - coeff) + coff * y1
        TimeVaryingInputs.evaluate!(dest, itp, date_pre, method)
        coeff = (time - date_pre) / (date_post - date_pre)
        dest .*= (1 - coeff)
        TimeVaryingInputs.evaluate!(tmp_field2, itp, date_post, method)
        tmp_field2 .*= coeff
        dest .+= tmp_field2
        return nothing
    end
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time::Number,
    method::LinearPeriodFillingInterpolation,
    args...;
    kwargs...,
)
    TimeVaryingInputs.evaluate!(
        dest,
        itp,
        Dates.Millisecond(round(1_000 * time)) + itp.data_handler.start_date,
        args...,
        kwargs...,
    )
    return nothing
end

function TimeVaryingInputs.evaluate!(
    dest,
    itp::InterpolatingTimeVaryingInput23D,
    time::ITime,
    method::LinearPeriodFillingInterpolation,
    args...;
    kwargs...,
)
    TimeVaryingInputs.evaluate!(dest, itp, date(time), args..., kwargs...)
    return nothing
end
