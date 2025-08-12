"""
    var_missing(x; val = -9999)

A function that checks if the value of `x` is equal to 
-9999, which is the value that Fluxnet uses for missing 
data. Returns true if x == -9999.
"""
var_missing(x; val = -9999) = x == val

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

Note that this function handles missing data by removing it (assuming it is marked by missing by a given value equal to `val`), because the 
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
