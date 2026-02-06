"""
    var_missing(x; val = -9999)

A function that checks if the value of `x` is equal to
-9999, which is the value that Fluxnet uses for missing
data. Returns true if x == -9999.
"""
var_missing(x; val = -9999) = x == val

"""
    hour_offset_to_period(hour_offset_from_UTC)

Convert a numerical hour offset (which may be fractional) to a
Dates.Period that can be added to DateTime objects.
Supports both integer and fractional time zones (e.g., 5.5 hours).
"""
function hour_offset_to_period(hour_offset_from_UTC)
    hours = floor(Int, hour_offset_from_UTC)
    minutes = round(Int, (hour_offset_from_UTC - hours) * 60)
    return Dates.Hour(hours) + Dates.Minute(minutes)
end

"""
   mask_data(t, v; val = -9999)

Returns t[id], v[id], where id corresponds
to the indices where v is not equal `val`,
which indicates missing data.

Note that `v` can either be a 1d array or a 2d array;
in the latter case, id refers to rows where
none of the columns of `v` are missing. The returned
array is 2d: v[id, :]
"""
function mask_data(t, v; val = -9999)
    if ndims(v) == 1
        not_missing_mask = .~var_missing.(v; val)
        return t[not_missing_mask], v[not_missing_mask]
    elseif ndims(v) == 2
        not_missing = .~var_missing.(v; val)
        not_missing_mask = all(not_missing, dims = 2)[:]
        return t[not_missing_mask], v[not_missing_mask, :]
    else
        @error("Dimension of data is not 1 or 2.")
    end
end

"""
    time_varying_input_from_data(
        data,
        varname::String,
        column_name_map,
        time_in_seconds::Vector;
        preprocess_func = identity,
        val=-9999,
    )

Returns the TimeVaryingInput object corresponding
to `varname` in the data matrix, the `column_name_map`
should map between varname and column id, and the `time_in_seconds`
should be the timestamp in seconds relative to the start_date of the
simulation corresponding to each row of data.

If you need to preprocess the data (e.g., unit conversion), you must pass
a pointwise function preprocess_func(var) as a keyword argument.

Note that this function handles missing data by removing it (assuming it is
marked by missing by a given value equal to `val`), because the
TimeVaryingInput object is an interpolating object (in time).
"""
function time_varying_input_from_data(
    data,
    varname::String,
    column_name_map,
    time_in_seconds::Vector;
    preprocess_func = identity,
    val = -9999,
)
    var_data = data[:, column_name_map[varname]]
    # The time varying input object interpolates over gaps
    # as needed, so we remove data that is marked as missing here
    t, v = mask_data(time_in_seconds, var_data; val)
    return TimeVaryingInput(t, preprocess_func.(v))
end

"""
    time_varying_input_from_data(
        data,
        varnames::Vector{String},
        column_name_map,
        time_in_seconds::Vector,
        preprocess_func,
        val=-9999,
    )

Returns a TimeVaryingInput object which is computing using
`preprocess_func` as a pointwise function of -in order - the columns in
`data` specified by `varnames`.

For example, if you wish to compute specific humidity from temperature,
pressure, and vpd, you would do:
varnames = ["TA_F", "PA_F", "VPD_F"] along with a preprocess function of
the form
function preprocess_func(T,P,VPD)
# carries out unit conversion
# computes q from T, P, VPD
# return q
end

The `column_name_map`
should map between varname and column id, and the `time_in_seconds`
should be the timestamp in seconds relative to the start_date of the
simulation corresponding to each row of data.

Note that this function handles missing data by removing it (assuming it is marked by missing by a given value equal to `val`), because the
TimeVaryingInput object is an interpolating object (in time).
"""
function time_varying_input_from_data(
    data,
    varnames::Vector{String},
    column_name_map,
    time_in_seconds::Vector;
    preprocess_func = identity,
    val = -9999,
)
    var_ids = [column_name_map[varname] for varname in varnames]
    var_data = data[:, var_ids]
    # The time varying input object interpolates over gaps
    # as needed, so we remove data that is marked as missing here
    t, v = mask_data(time_in_seconds, var_data; val)
    return TimeVaryingInput(t, preprocess_func.(eachcol(v)...))
end

"""
     get_data_at_start_date(
        v::Vector,
        Δ_date::Vector;
        preprocess_func = identity,
        val = -9999,
    )

Returns the value in the raw data `v` closest to where
|Δ_date| = 0, after preprocessing the data using preprocess_func
(a pointwise function) and removing missing values.

If Δ_date corresponds to a vector of dates relative to the start_date
of the simulation, the returned value can be used as an initial condition.
"""
function get_data_at_start_date(
    v::Vector,
    Δ_date::Vector;
    preprocess_func = identity,
    val = -9999,
)
    Δ_date, v = mask_data(Δ_date, v; val)
    idx_start = argmin(abs.(Δ_date))
    return preprocess_func(v[idx_start])
end

"""
    read_fluxnet_data(site_ID)

Reads Fluxnet CSV data for the provided site, and returns a tuple `(data, columns)`
where `data` is a `n_rows x n_cols` Matrix of the data, and `columns` is a
`1 x n_cols` Matrix of column names.

This function now supports both standard Fluxnet sites and CalMIP sites.
For CalMIP sites (e.g., "DK-Sor"), it uses the callmip_data_path function.
"""
function read_fluxnet_data(site_ID)
    # Check if this is a CalMIP site
    callmip_sites = ("DK-Sor",)  # Add more CalMIP sites here as needed
    
    if site_ID ∈ callmip_sites
        fluxnet_csv_path = ClimaLand.Artifacts.callmip_data_path(site_ID)
    else
        fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    end
    
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    return (data, columns)
end

"""
    get_UTC_datetimes(hour_offset_from_UTC, data, column_name_map; timestamp_name = "TIMESTAMP")

Converts local timestamps from FLUXNET data to UTC timestamps.

Fluxnet as has three timestamps associated with it: TIMESTAMP, TIMESTAMP_START, and TIMESTAMP_END,
with the first referring to an unspecified timestamp of the data, the second referring to the start
of the averaging period, and the third the end.

For forcing data, we want to get the timestamps at the start and end of each averaging period,
then average them to get the timestamp at the midpoints. This can be done by calling this function
twice, once with `timestamp_name = "TIMESTAMP_START"` and once with `"TIMESTAMP_END"`,
then averaging the results.

For comparison data, we want to get the timestamp at the start of the averaging period (described below).
This can be done by calling this function with `timestamp_name = "TIMESTAMP_START"`.

ClimaLand diagnostics are reduced over a time period (e.g. hourly, daily, monthly), and saved with the
first date of period. For example, the hourly average from 11-noon is saved with a
timestamp of 11. To make a true comparison to Fluxnet data, therefore, we must use halfhourly
diagnostics in ClimaLand, and return the UTC time that corresponds to the start of the averaging period in Fluxnet.

# Arguments
- `hour_offset_from_UTC`: The hour offset from UTC for the local timezone
        (may be fractional)
- `data`: The data matrix containing timestamp information
- `column_name_map`: Dictionary mapping variable names to column indices in the data
- `timestamp_name`: The timestamp column name to use
        (options: "TIMESTAMP_START" (default), "TIMESTAMP_END", "TIMESTAMP")

# Returns
- `UTC_datetimes`: Vector of UTC DateTime objects

# Example
```julia
UTC_datetimes = get_UTC_datetimes(-5, data, column_name_map)  # EST to UTC
```
"""
function get_UTC_datetimes(
    hour_offset_from_UTC,
    data,
    column_name_map;
    timestamp_name = "TIMESTAMP_START",
)
    @assert timestamp_name in ("TIMESTAMP", "TIMESTAMP_START", "TIMESTAMP_END") "Invalid timestamp column name $timestamp_name"

    col = column_name_map[timestamp_name]
    local_datetime = DateTime.(string.(Int.(data[:, col])), "yyyymmddHHMM")
    UTC_datetimes =
        local_datetime .- hour_offset_to_period(hour_offset_from_UTC)
    return UTC_datetimes
end

"""
    get_column_name_map(varnames, columns; error_on_missing = true)

Returns a dictionary mapping variable names to column indices in the data.

The option `error_on_missing` is used to control whether to error if any variables are missing.
If `error_on_missing` is false, missing variable names will be mapped to `nothing`.

# Arguments
- `varnames`: Vector of variable names
- `columns`: Vector of column names
- `error_on_missing`: Whether to error if any variables are missing.
    For forcing data, we need all variables to be present so we want to error if any are missing.
    For comparison data, we don't need all variables to be present so we don't want to error.

# Returns
- `column_name_map`: Dictionary mapping variable names to column indices
"""
function get_column_name_map(varnames, columns; error_on_missing = true)
    column_name_map = Dict(
        varname => findfirst(columns[:] .== varname) for varname in varnames
    )

    # If any variables are missing, error if requested
    if error_on_missing
        missing_vars = varnames[findall(values(column_name_map) .== nothing)]
        @assert isempty(missing_vars) "Data is missing required variables: $missing_vars"
    end
    return column_name_map
end

"""
    get_comparison_data(
        data::Matrix,
        varname::String,
        column_name_map::Dict,
        climaland_shortname::String;
        preprocess_func = identity,
        val = -9999,
    )

Gets and returns a NamedTuple with the data identified
by `varname` in the `data` matrix by looking up the
column index of varname
using the column_name_map, replacing missing data with the mean
of the non-missing data, and preprocessing the data using the
`preprocess_func`, which should be a pointwise function. The
key of the NamedTuple should be the shortname of the corresponding
variable using the shortname of ClimaDiagnostics: `climaland_shortname`.

If the data column
is missing completely, or if the column is present
but the data is missing in each row, an empty NamedTuple
is returned.

In the future, we can explore dropping the missing values to be
consistent with what we do above, but this is consistent with the
current Fluxnet runs.
"""
function FluxnetSimulations.get_comparison_data(
    data::Matrix,
    varname::String,
    column_name_map::Dict,
    climaland_shortname::String;
    preprocess_func = identity,
    val = -9999,
)
    idx = column_name_map[varname]
    if isnothing(idx)
        return (;)
    else
        v = data[:, idx]
        missing_mask = var_missing.(v; val)
        not_missing_mask = .~missing_mask
        n_missing = sum(missing_mask)
        if n_missing < length(v) # some data are present
            v[missing_mask] .= sum(v[not_missing_mask]) / sum(not_missing_mask)
            return (; Symbol(climaland_shortname) => preprocess_func.(v))
        else
            return (;)
        end
    end
end
