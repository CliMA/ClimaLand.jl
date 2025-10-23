"""
    var_missing(x; val = -9999)

A function that checks if the value of `x` is equal to
-9999, which is the value that Fluxnet uses for missing
data. Returns true if x == -9999.
"""
var_missing(x; val = -9999) = x == val

"""
    hour_offset_to_period(hour_offset_from_UTC)

Convert a hour offset (which may be fractional) to a Dates.Period
that can be added to DateTime objects.
Supports both integer and half-hour time zones (e.g., 5.5 hours).
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
    get_required_varnames(construct_prescribed_ground::Bool)

Returns the tuple of required variable names for Fluxnet data processing.
If `construct_prescribed_ground` is true, includes soil moisture and temperature variables.
"""
function get_required_varnames(construct_prescribed_ground::Bool)
    varnames = (
        "TIMESTAMP_START",
        "TIMESTAMP_END",
        "TA_F",
        "VPD_F",
        "PA_F",
        "P_F",
        "WS_F",
        "LW_IN_F",
        "SW_IN_F",
        "CO2_F_MDS",
    )
    construct_prescribed_ground &&
        (varnames = (varnames..., "SWC_F_MDS_1", "TS_F_MDS_1"))
    return varnames
end

"""
    create_column_name_map(columns, varnames; is_comparison_data)

Creates a dictionary mapping variable names to their column indices in the data.
The input `columns` is expected to be a vector of column names from the data,
and `varnames` is a tuple of variable names to look for.

If `is_comparison_data` is false (i.e. we're processing forcing data) and any
input varname is missing in the columns, an error is logged.

Returns a Dict{String, Union{Int, Nothing}} where each key is a varname
and the value is the corresponding column index.
"""
function create_column_name_map(columns, varnames; is_comparison_data)
    column_name_map = Dict(
        varname => findfirst(columns[:] .== varname) for varname in varnames
    )

    # Check for missing columns, if requested
    if !is_comparison_data
        nothing_id = findall(collect(values(column_name_map)) .== nothing)
        if !isempty(nothing_id)
            missing_vars = [varnames[i] for i in nothing_id]
            @error("$(missing_vars) is missing in the data, but required.")
        end
    end

    return column_name_map
end

"""
    find_rows_with_all_variables_available(data, column_name_map, varnames; val = -9999)

Return a boolean vector marking rows where all required variables are available as `true`.
For each name in `varnames`, the column index is looked up in `column_name_map`,
and a row is `true` only if none of those columns equals `val` in that row.

Arguments
- data: matrix with observations in rows and variables in columns.
- column_name_map: Dict mapping variable names (String) to column indices (Int).
- varnames: tuple of variable names to check.
- val: sentinel used to indicate missing data (default -9999).

Returns
- Vector{Bool} of length size(data, 1), with `true` indicating rows where all
    required variables are available.
"""
function find_rows_with_all_variables_available(
    data,
    column_name_map,
    varnames;
    val = -9999,
)
    # Get the column indices for the input variables
    cols = getindex.(Ref(column_name_map), varnames)

    # Row is valid (true) only if all columns are not equal to val
    valid_rows = trues(size(data, 1))
    @inbounds for col_idx in cols
        valid_rows .&= (data[:, col_idx] .!= val)
    end
    return valid_rows
end

"""
    read_fluxnet_data(
        site_ID;
        construct_prescribed_ground = false,
        varnames = get_required_varnames(construct_prescribed_ground),
        is_comparison_data = false,
        hour_offset_from_UTC = nothing,
    )

Reads Fluxnet CSV data and processes timestamps, converting local time to UTC.
Returns (data, column_name_map, UTC_datetime).

If `hour_offset_from_UTC` is nothing, the UTC datetime is not returned.

Fluxnet as has two timestamps associated with it: TIMESTAMP_START and TIMESTAMP_END,
with the first referring to the start of the averaging period, and the second the end.

The flag `is_comparison_data` controls which timestamp is used to construct the local and UTC
datetimes, and whether all variables are required to be present in the data.

For forcing data, we want to use a timestamp corresponding to the middle of the period
`[TIMESTAMP_START, TIMESTAMP_END]`, since the data is provided as an average over that period.
For comparison data, we want to use the `TIMESTAMP_END` value to compare to ClimaLand
diagnostics, which are saved at the end of the averaging period.
For example, the hourly average from 10:00-11:00 is saved with a timestamp of 11:00.
To make a true comparison to Fluxnet data, therefore, we must use halfhourly diagnostics in ClimaLand,
and return the UTC time that corresponds to the end of the averaging period in Fluxnet.
"""
function read_fluxnet_data(
    site_ID;
    hour_offset_from_UTC,
    construct_prescribed_ground = false,
    varnames = get_required_varnames(construct_prescribed_ground),
    is_comparison_data = false,
)
    # Read the data from the CSV file for this site
    fluxnet_csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
    (data, columns) = readdlm(fluxnet_csv_path, ','; header = true)

    # Get map from variable name to column index
    column_name_map =
        create_column_name_map(columns, varnames; is_comparison_data)

    # Find rows with all required variables available
    valid_rows =
        find_rows_with_all_variables_available(data, column_name_map, varnames)

    # Get local datetime at the end or middle of the averaging period
    local_datetime_end =
        DateTime.(
            string.(Int.(data[valid_rows, column_name_map["TIMESTAMP_END"]])),
            "yyyymmddHHMM",
        )
    if is_comparison_data
        local_datetime = local_datetime_end
    else
        local_datetime_start =
            DateTime.(
                string.(
                    Int.(data[valid_rows, column_name_map["TIMESTAMP_START"]])
                ),
                "yyyymmddHHMM",
            )
        local_datetime =
            local_datetime_start .+
            (local_datetime_end .- local_datetime_start) ./ 2
    end

    # Convert to UTC if hour offset is given
    if isnothing(hour_offset_from_UTC)
        UTC_datetime = nothing
    else
        UTC_datetime =
            local_datetime .+ hour_offset_to_period(hour_offset_from_UTC)
    end

    return (data, column_name_map, UTC_datetime)
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
