import Dates: Hour, Day, Minute, Second

"""
    monthly_maxs(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the monthly max for the given variables.
"""
monthly_maxs(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Month(1), max, output_writer, short_names...)

"""
    monthly_mins(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the monthly min for the given variables.
"""
monthly_mins(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Month(1), min, output_writer, short_names...)

"""
    monthly_averages(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the monthly average for the given variables.
"""
# An average is just a sum with a normalization before output
monthly_averages(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        Month(1),
        (+),
        output_writer,
        start_date,
        short_names...;
        pre_output_hook! = average_pre_output_hook!,
    )

"""
    tendaily_maxs(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the max over ten days for the given variables.
"""
tendaily_maxs(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        10 * Day(1),
        max,
        output_writer,
        start_date,
        short_names...,
    )

"""
    tendaily_mins(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the min over ten days for the given variables.
"""
tendaily_mins(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        10 * Day(1),
        min,
        output_writer,
        start_date,
        short_names...,
    )

"""
    tendaily_averages(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the average over ten days for the given variables.
"""
# An average is just a sum with a normalization before output
tendaily_averages(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        10 * Day(1),
        (+),
        output_writer,
        start_date,
        short_names...;
        pre_output_hook! = average_pre_output_hook!,
    )

"""
    daily_maxs(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the daily max for the given variables.
"""
daily_maxs(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Day(1), max, output_writer, start_date, short_names...)

"""
    daily_mins(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the daily min for the given variables.
"""
daily_mins(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Day(1), min, output_writer, start_date, short_names...)

"""
    daily_averages(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the daily average for the given variables.
"""
# An average is just a sum with a normalization before output
daily_averages(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        Day(1),
        (+),
        output_writer,
        start_date,
        short_names...;
        pre_output_hook! = average_pre_output_hook!,
    )

"""
    hourly_maxs(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the hourly max for the given variables.
"""
hourly_maxs(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Hour(1), max, output_writer, start_date, short_names...)

"""
    hourly_mins(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the hourly min for the given variables.
"""
hourly_mins(FT, short_names...; output_writer, start_date) =
    common_diagnostics(Hour(1), min, output_writer, start_date, short_names...)

# An average is just a sum with a normalization before output
"""
    hourly_averages(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the hourly average for the given variables.
"""
hourly_averages(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        Hour(1),
        (+),
        output_writer,
        start_date,
        short_names...;
        pre_output_hook! = average_pre_output_hook!,
    )

"""
    halfhourly_averages(FT, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the 30 minute average for the given variables.
"""
halfhourly_averages(FT, short_names...; output_writer, start_date) =
    common_diagnostics(
        Minute(30),
        (+),
        output_writer,
        start_date,
        short_names...;
        pre_output_hook! = average_pre_output_hook!,
    )

"""
    every_dt_inst(FT, dt, short_names...; output_writer,  start_date)

Return a list of `ScheduledDiagnostics` that compute the instantaneous value
at every time step for the given variables.
"""
every_dt_inst(FT, dt, short_names...; output_writer, start_date) =
    common_diagnostics(
        Second(dt),
        nothing, # reduction
        output_writer,
        start_date,
        short_names...;
    )
