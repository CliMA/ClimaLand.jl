#=
The following functions are used for scraping data:
    - `sitedata_daily`
    - `sitedata_hourly`
    - `snotel_metadata`
    - `scrape_site_paper` (follow same process as paper)

The following functions are used for filtering/processing scraped data:
    - `apply_bounds`
    - `hourly2daily`
    - `rectify_daily_hourly`
    - `scale_cols`
    - `makediffs`
    - `rolldata`

The following functions are used for quality control:
    - `serreze_qc`
    - `yan_qc`
    - `fix_temp`
    - `bcqc_daily` (wrapper around all steps for daily data)
    - `bcqc_hourly` (wrapper around all steps for hourly data)
    - `d_impute`
    - `qc_filter`

The following functions are used for data preprocessing for a model:
    - `snowsplit`
    - `train_filter`
    - `prep_data`
    - `make_data`
=#
module DataTools
using DataFrames
using Dates, Statistics
export sitedata_daily,
    sitedata_hourly,
    apply_bounds,
    hourly2daily,
    rectify_daily_hourly,
    scale_cols,
    makediffs,
    rolldata,
    snowsplit,
    train_filter,
    prep_data,
    make_data,
    snotel_metadata,
    serreze_qc,
    yan_qc,
    fix_temp,
    bcqc_daily,
    bcqc_hourly,
    d_impute,
    qc_filter,
    scrape_site_paper

include("./WebTools.jl")

"""
    missingfirst(input)

Helper function for getting function value among missing values. Equivalent to:
    `isempty(collect(skipmissing(input))) ? missing : first(skipmissing(input))`
"""
missingfirst(input) =
    isempty(collect(skipmissing(input))) ? missing : first(skipmissing(input))

"""
    missingmean(input)

Helper function for getting function value among missing values. Equivalent to:
    `isempty(collect(skipmissing(input))) ? missing : mean(skipmissing(input))`
"""
missingmean(input) =
    isempty(collect(skipmissing(input))) ? missing : mean(skipmissing(input))

"""
    missingmax(input)

Helper function for getting function value among missing values. Equivalent to:
    `isempty(collect(skipmissing(input))) ? missing : maximum(skipmissing(input))`
"""
missingmax(input) =
    isempty(collect(skipmissing(input))) ? missing : maximum(skipmissing(input))

"""
    data_url(id, state, fields; metric, hourly, start_date, end_date)

Return the url to gather SNOTEL data from a particular station.

# Arguments
- `id::AbstractString`: the SNOTEL ID code of the station to access.
- `state::AbstractString`: the state abbreviation of the station to access, e.g., `"MT"` for Montana.
- `fields::Vector{<:AbstractString}`: a list of strings specifying the SNOTEL parameters to collect.
- `metric::Bool`: a boolean specifying whether to gather data in metric (true) or imperial (false) units. Default is `false`.
- `hourly::Bool`: a boolean specifying whether to scrape hourly (true) or daily (false) data. Default is `false`.
- `start_date::AbstractString`: A string specifying a start date for scraping data, or `"start"` for "first available data date"
- `end_date::AbstractString`: A string specifying an end date for scraping data, or `"end"` for `"2024-02-01"`
"""
function data_url(
    id::Int,
    state::AbstractString,
    fields::Vector{<:AbstractString};
    metric::Bool = false,
    hourly::Bool = false,
    start_date::AbstractString = "start",
    end_date::AbstractString = "end",
)
    start = (start_date == "start") ? "1850-01-01" : start_date
    fin = (end_date == "end") ? "2024-02-01" : end_date
    link = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customSingleStationReport"
    link =
        link *
        (metric ? ",metric/" : "/") *
        (hourly ? "hourly" : "daily") *
        "/start_of_period/"
    link = link * string(id) * ":" * state * ":SNTL/" * start * "," * fin * "/"
    for field in fields
        link = link * field * "::value,"
    end
    return link
end


"""
    get_data(id, state, fields; metric, hourly, start, finish)

Return the SNOTEL data from a station as a `DataFrame`.

# Arguments
- `id::AbstractString`: the SNOTEL ID code of the station to access.
- `state::AbstractString`: the state abbreviation of the station to access.
- `fields::Vector{<:AbstractString}`: a list of strings specifying the SNOTEL parameters to collect.
- `metric::Bool = false`: a boolean specifying whether to gather data in metric (true) or imperial (false) units.
- `hourly::Bool = false`: a boolean specifying whether to scrape hourly (true) or daily (false) data.
- `start::AbstractString`: Optional string to specify starting date of data collection. Default is `"start"`.
- `finish::AbstractString`: Optional string to specify ending date of data collection. Default is `"end"`.
"""
function get_data(
    id::Int,
    state::AbstractString,
    fields::Vector{<:AbstractString};
    metric::Bool = false,
    hourly::Bool = false,
    start::AbstractString = "start",
    finish::AbstractString = "end",
)
    url = data_url(
        id,
        state,
        fields,
        metric = metric,
        hourly = hourly,
        start_date = start,
        end_date = finish,
    )
    #data = CSV.read(HTTP.get(url).body, DataFrame, comment = "#", delim = ",")
    data = df_from_url(url)
    return data
end


"""
    sitedata_daily(id, state; imp_fields, metric_fields, colnames, start, finish)

Return the daily SNOTEL data from a station as a `DataFrame`.

# Arguments
- `id::AbstractString`: the id code of the station to access.
- `state::AbstractString`: the state abbreviation of the station to access.
- `imp_fields::Vector{<:AbstractString}`: Parameters to return in imperial units. Default is `["WTEQ", "SNWD", "PREC"]`.
-  `metric_fields::Vector{<:AbstractString}`: Parameters to return in metric units. Default is `["RHUMV", "SRADV", "WSPDV", "TAVG"]`.
- `colnames::Vector{<:AbstractString}`: Optional column names to change header after scraping data. Default follows that of the paper,
which is `["date", "SWE", "z", "precip", "rel_hum_avg", "sol_rad_avg", "wind_speed_avg", "air_temp_avg"]`.
- `start::AbstractString`: Optional string to specify starting date of data collection. Default is `"start"`.
- `finish::AbstractString`: Optional string to specify ending date of data collection. Default is `"end"`.
"""
function sitedata_daily(
    id::Int,
    state::AbstractString;
    imp_fields::Vector{<:AbstractString} = ["WTEQ", "SNWD", "PREC"],
    metric_fields::Vector{<:AbstractString} = [
        "RHUMV",
        "SRADV",
        "WSPDV",
        "TAVG",
    ],
    colnames::Vector{<:AbstractString} = [
        "date",
        "SWE",
        "z",
        "precip",
        "rel_hum_avg",
        "sol_rad_avg",
        "wind_speed_avg",
        "air_temp_avg",
    ],
    start::AbstractString = "start",
    finish::AbstractString = "end",
)
    imp_data = get_data(id, state, imp_fields, start = start, finish = finish)
    metric_data = get_data(
        id,
        state,
        metric_fields,
        metric = true,
        start = start,
        finish = finish,
    )
    data = DataFrames.innerjoin(imp_data, metric_data, on = :Date)
    DataFrames.rename!(data, Symbol.(colnames))
    return data
end


"""
    sitedata_hourly(id, state; imp_fields, metric_fields, colnames, start, finish)

Return the hourly SNOTEL data from a station as a `DataFrame`.

# Arguments
- `id::AbstractString`: the id code of the station to access.
- `state::AbstractString`: the state abbreviation of the station to access.
- `imp_fields::Vector{<:AbstractString}`: Parameters to return in imperial units. Default is `["WTEQ", "SNWD", "PREC"]`.
-  `metric_fields::Vector{<:AbstractString}`: Parameters to return in metric units. Default is `["RHUMV", "SRADV", "WSPDV", "TAVG"]`.
- `colnames::Vector{<:AbstractString}`: Optional column names to change header after scraping data. Default follows that of the paper,
which is `["date", "SWE", "z", "precip", "rel_hum_avg", "sol_rad_avg", "wind_speed_avg", "air_temp_avg"]`.
- `start::AbstractString`: Optional string to specify starting date of data collection. Default is `"start"`.
- `finish::AbstractString`: Optional string to specify ending date of data collection. Default is `"end"`.
"""
function sitedata_hourly(
    id::Int,
    state::AbstractString;
    imp_fields::Vector{<:AbstractString} = ["WTEQ", "SNWD", "PREC"],
    metric_fields::Vector{<:AbstractString} = [
        "RHUMV",
        "SRADV",
        "WSPDV",
        "TOBS",
    ],
    colnames::Vector{<:AbstractString} = [
        "date",
        "SWE",
        "z",
        "precip",
        "rel_hum_avg",
        "sol_rad_avg",
        "wind_speed_avg",
        "air_temp_avg",
    ],
    start::AbstractString = "start",
    finish::AbstractString = "end",
)
    imp_data = get_data(
        id,
        state,
        imp_fields,
        hourly = true,
        start = start,
        finish = finish,
    )
    metric_data = get_data(
        id,
        state,
        metric_fields,
        metric = true,
        hourly = true,
        start = start,
        finish = finish,
    )
    data = DataFrames.innerjoin(imp_data, metric_data, on = :Date)
    DataFrames.rename!(data, Symbol.(colnames))
    data[!, :date] .= DateTime.(data[!, :date], "yyyy-mm-dd HH:MM")
    return data
end


"""
    apply_bounds(data, bounds)

Threshold data in a min/max form from a database using a dictionary of specified bounds.
Values outside bounds are converted to 'missing'.

# Arguments
- `data::DataFrame`: the data over which to apply bounds.
- `bounds::Dict{Symbol, Tuple{Real, Real}}`: The dictionary specifying which columns to threshold, as
well as the thresholds to apply via (min, max).
"""
function apply_bounds(data::DataFrame, bounds::Dict{Symbol, Tuple{Real, Real}})
    bounded = deepcopy(data)
    for colname in names(bounded)
        coltype = eltype(bounded[!, Symbol(colname)])
        bounded[!, Symbol(colname)] = convert(
            Vector{Union{Missing, coltype}},
            bounded[!, Symbol(colname)],
        )
    end
    for key in keys(bounds)
        below_min = findall(.<(bounds[key][1]), skipmissing(bounded[!, key]))
        above_max = findall(.>(bounds[key][2]), skipmissing(bounded[!, key]))
        if !isempty(above_max)
            bounded[above_max, key] .= missing
        end
        if !isempty(below_min)
            bounded[below_min, key] .= missing
        end
    end
    return bounded
end


"""
    hourly2daily(hourlydata; meanvars, firstvars, blockhours = true, blocknum = 3)

Convert an hourly SNOTEL dataset to a daily SNOTEL dataset, averaging humidity, radiation, wind, and temperature
data, but maintaining start-of-window SWE, z, and precipitation data.
Makes use of custom `missingmean` and `missingfirst` functions to handle missing values when finding the mean and first non-missing
values, respectively, when working with hourly datasets with missing values.

**Note: requires time field to be named "date". Does not extend to data beyond that for the paper without usage of optional args.

# Arguments
- `hourlydata::DataFrame`: the hourly data over which to average by-day.
- `meanvars::Vector{Symbol}`: the variables to average over
- `firstvars::Vector{Symbol}`: the variables to take the first available value
- `blockhours::Bool`: indicating whether to take values in hourly chunks (to avoid biases depending on which hours are available). Default is `true`.
- `blocknum::Int`: indicates the number of chunks to use when chunking by hour (Default is 3)
"""
function hourly2daily(
    hourlydata::DataFrame;
    meanvars = [:rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg],
    firstvars = [:date, :SWE, :z, :precip],
    blockhours = true,
    blocknum = 3,
)
    data = deepcopy(hourlydata)
    divisor = div(24, blocknum)
    block =
        Dates.date2epochdays.(Date.(data[!, :date])) .+
        ((div.(Dates.hour.(data[!, :date]), divisor) .+ 1) ./ (blocknum + 1))
    data[!, :block] .= block
    data[!, :date] .= Date.(data[!, :date])
    groupvar = blockhours ? :block : :date
    id = data[1, :id]
    dailydata = combine(
        groupby(data, groupvar),
        firstvars .=> missingfirst,
        meanvars .=> missingmean,
        renamecols = false,
    )
    if blockhours
        dailydata = combine(
            groupby(dailydata, :date),
            firstvars .=> missingfirst,
            meanvars .=> mean,
            renamecols = false,
        )
    end
    dailydata[!, :id] .= id
    return dailydata
end

"""
    rectify_daily_hourly(daily_data, hourly_data; prioritize_hour, daily_only, hourly_only)

Use one SNOTEL data stream (hourly data) to fill holes in another SNOTEL data stream (daily data).

**Note: requires time column to be named "date".

# Arguments
- `daily_data::DataFrame`: the main (daily) data over which to fill missing holes.
- `hourly_data::DataFrame`: the (converted to daily-format, see hourly2daily) hourly data over which to fill holes in the daily data.
- `prioritize_hour::Vector{Symbol}`: an optional list of fields for which to prioritize the hourly data over the daily when coalescing.
- `daily_only::Vector{Symbol}`: an optional list of fields for which to only use daily data when coalescing.
    Default is `[:SWE, :precip]`
- `hourly_only::Vector{Symbol}`: an optional list of fields for which to only use hourly data when coalescing.
    Default is `[:sol_rad_avg, :wind_speed_avg, :rel_hum_avg]`
"""
function rectify_daily_hourly(
    daily_data::DataFrame,
    hourly_data::DataFrame;
    prioritize_hour::Vector{Symbol} = Symbol[],
    daily_only::Vector{Symbol} = [:SWE, :precip],
    hourly_only::Vector{Symbol} = [:sol_rad_avg, :wind_speed_avg, :rel_hum_avg],
)
    vars = names(daily_data)
    combined = outerjoin(daily_data, hourly_data, on = :date, makeunique = true)
    for var in vars[2:end]
        if Symbol(var) in daily_only
            continue
        elseif Symbol(var) in hourly_only
            combined[!, Symbol(var)] .= combined[!, Symbol(var * "_1")]
        elseif Symbol(var) in prioritize_hour
            combined[!, Symbol(var)] .=
                coalesce.(
                    combined[!, Symbol(var * "_1")],
                    combined[!, Symbol(var)],
                )
        else
            combined[!, Symbol(var)] .=
                coalesce.(
                    combined[!, Symbol(var)],
                    combined[!, Symbol(var * "_1")],
                )
        end
    end
    return sort(select(combined, vars), :date)
end

"""
    scale_cols(data, scales)

Apply multiplicative scaling to select columns of a data frame.

# Arguments
- `data::DataFrame`: the data frame over which to apply the scaling.
- `scales::Dict{Symbol, Real}`: The dictionary specifying which columns to scale, as
well as the multitiplicative constant to apply to that column.
"""
function scale_cols(data::DataFrame, scales::Dict{Symbol, Real})
    scaled = deepcopy(data)
    for key in keys(scales)
        scaled[!, key] .= scaled[!, key] .* scales[key]
    end
    return scaled
end

"""
    makediffs(data, Δt; diffvars)

Turn accumulated data fields into differential data fields.
Turns data fields which represent accumulated values of variables into a rate of change of that variable, 
`(data[i+1]-data[i])/(time[i+1]-time[i])`, but only in the case where `time[i+1]-time[i] = Δt`.

**Note: apply after scaling and getting rid of missing values. Assumes time column has name "date" and has `Date`/`DateTime` units, in sequential order.

# Arguments
- `data::DataFrame`: the data over which to apply differencing.
- `Δt::Period`: the amount of time representing one unit timestep in the data.
- `diffvars::Vector{Symbol}`: Columns to apply differencing to. Default is `[:SWE, :z, :precip]`.
"""
function makediffs(
    data::DataFrame,
    Δt::Period;
    diffvars::Vector{Symbol} = [:SWE, :z, :precip],
)
    diffdata = deepcopy(data[1:(end - 1), :])
    dt = data[2:end, :date] .- data[1:(end - 1), :date]
    good_t = (dt .== Δt)
    for var in diffvars
        dvar = data[2:end, var] .- data[1:(end - 1), var]
        diffdata[!, Symbol("d" * string(var) * "dt")] =
            dvar ./ Dates.value.(Second.(dt))
    end
    diffdata = diffdata[good_t, :]
    return diffdata
end

"""
    stack2DF(stack, colnames)

Convert a vector of vectors into a `DataFrame`, with specified column names.

# Arguments
- `stack::Vector{Vector{Any}}`: the data stack to convert.
- `colnames::Vector{<:AbstractString}`: The names to give the columns of the `DataFrame`.
"""
function stack2DF(stack, colnames)
    return DataFrame(mapreduce(permutedims, vcat, stack), colnames)
end

"""
    rolldata(data, Δt, N; takefirst)

Apply a moving average of N timesteps to all data, except for the variables
specified in `takefirst`, for which the leading value is maintained.

**Note: assumes the time column is named "date" and has `Date`/`DateTime` units, in sequential order.

# Arguments
- `data::DataFrame`: the data over which to apply averaging.
- `Δt::Period`: the amount of time representing one unit timestep in the data.
- `N::Int`: the number of intervals (timesteps) to include in the average
- `takefirst::Vector{Symbol}`: Columns to apply differencing to. Default is `[:date, :SWE, :z, :precip]`.
"""
function rolldata(
    data,
    Δt::Period,
    N::Int;
    takefirst::Vector{Symbol} = [:date, :SWE, :z, :precip, :id],
)
    dt = (data[N:end, :date] .- data[1:(end - N + 1), :date]) ./ Δt
    rolleddata = Any[]
    for i in 1:length(dt)
        if dt[i] <= N - 1
            row = Vector{Any}()
            for var in names(data)
                if Symbol(var) in takefirst
                    push!(row, data[i, Symbol(var)])
                else
                    push!(row, mean(data[i:(i + N - 1), Symbol(var)]))
                end
            end
            push!(rolleddata, row)
        end
    end
    df = stack2DF(rolleddata, names(data))
    for j in 1:size(df)[2]
        df[!, j] .= (eltype(data[!, j])).(df[!, j])
    end
    return df
end

"""
    snowsplit(air_temp, hum, precip)

Engineer total water content of precipitation into snow and rain portions,
accoridng to the paper outlined in https://www.nature.com/articles/s41467-018-03629-7.

**Note: assumes relative humidity is in range 0-1.

# Arguments
- `air_temp::Vector{<:Real}`: the air temperature data.
- `hum::Vector{<:Real}`: the relative humidity data.
- `precip::Vector{<:Real}`: the precipitation data.
"""
function snowsplit(
    air_temp::Vector{<:Real},
    hum::Vector{<:Real},
    precip::Vector{<:Real},
)
    α = -10.04
    β = 1.41
    γ = 0.09 * 100
    snow_frac = (1.0 ./ (1.0 .+ exp.(α .+ β .* air_temp .+ γ .* hum)))
    dprecipdt_snow = snow_frac .* precip
    dprecipdt_rain = (1.0 .- snow_frac) .* precip
    return (dprecipdt_snow, dprecipdt_rain)
end

"""
    train_filter(data; eps)

Apply filtering to training timeseries to generate the regression data for training.
Namely, remove days "without a snowpack" or small snowpack (z, SWE < eps, dprecipdt * (1 day) < eps),
and remove days when z != 0 but SWE = 0. Defauly eps is 0.005 m (5 cm, since z is discretized to the inch).

**Note: requires column names to match those of the paper for usage.

# Arguments
- `data::DataFrame`: the cleaned and processed data.
- `eps::Real`: a filtering threshold (in meters), set to 0.005.
"""
function train_filter(data::DataFrame; eps::Real = 0.005)
    zero_condition1 =
        (data[!, :SWE] .< eps) .&
        (data[!, :z] .< eps) .&
        (data[!, :dprecipdt] .< eps / 86400.0)
    data = data[Not(zero_condition1), :]
    zero_condition2 = (data[!, :z] .!= 0.0f0) .& (data[!, :SWE] .== 0.0f0)
    data = data[Not(zero_condition2), :]
    return data
end

"""
    prep_data(data; extract_vars, make_snow_split, physical_filter, eps)

Prepare a cleaned (scaled & gap-filled, potentially rolled) data stream or non-scraped data stream for model usage.

**Note: Requires column names to match those of the paper for usage.

# Arguments
- `data::DataFrame`: the cleaned and processed data.
- `extract_vars::Vector{Symbol}`: The list of columns to be used in the model. Default is the variables used in the paper.
- `make_snow_split::Bool`: Boolean indicating whether to split precipitation into water and snow. Default is `true`.
- `train_filt::Bool`: Boolean indicating whether to filter data using train_filter. Default is `true`.
- `eps::Real`: a filtering threshold, set to 0.005 m (see `train_filter()`).
"""
function prep_data(
    data::DataFrame;
    extract_vars::Vector{Symbol} = [
        :date,
        :id,
        :z,
        :SWE,
        :rel_hum_avg,
        :sol_rad_avg,
        :wind_speed_avg,
        :air_temp_avg,
        :dprecipdt_snow,
        :dprecipdt_rain,
        :dzdt,
        :dSWEdt,
    ],
    make_snow_split::Bool = true,
    train_filt::Bool = true,
    eps::Real = 0.005,
)
    if train_filt
        data = train_filter(data, eps = eps)
    end
    if make_snow_split
        (dprecipdt_snow, dprecipdt_rain) = snowsplit(
            data[!, :air_temp_avg],
            data[!, :rel_hum_avg],
            data[!, :dprecipdt],
        )
        data[!, :dprecipdt_snow] = dprecipdt_snow
        data[!, :dprecipdt_rain] = dprecipdt_rain
    end
    return select(data, extract_vars)
end

"""
    make_data(data, input_vars, target, target_scale; testidx, dtype)

Create training and testing matrices for the models from prepped input data, to be passed to a Flux `DataLoader` call for batching.

# Arguments
- `data::DataFrame`: the data to be used for training and validation data.
- `input_vars::Vector{Symbol}`: The set of features to be extracted as input data
- `target::Symbol`: The feature to be extracted as the target variable
- `target_scale::Real`: The scaling to apply to the output data
- `testidx': a boolean array to demarcate indices of testing/validation data, or `nothing`. Default is `nothing`
- `dtype::Type`: The data type consistent with the model. Default is `Float32`.
"""
function make_data(
    data::DataFrame,
    input_vars::Vector{Symbol},
    target::Symbol,
    target_scale::Real;
    testidx = nothing,
    dtype::Type = Float32,
)
    vardata = select(data, cat(dims = 1, input_vars, target))
    traindata = (isnothing(testidx)) ? vardata : vardata[(!).(testidx), :]
    x_train = Matrix{dtype}(traindata[!, Not(target)])'
    y_train = Vector{dtype}(traindata[!, target])' ./ dtype(target_scale)
    if isnothing(testidx)
        return Matrix{dtype}(x_train), Matrix{dtype}(y_train)
    end
    testdata = vardata[testidx, :]
    x_test = Matrix{dtype}(testdata[!, Not(target)])'
    y_test = Vector{dtype}(testdata[!, target])' ./ dtype(target_scale)
    return Matrix{dtype}(x_train),
    Matrix{dtype}(y_train),
    Matrix{dtype}(x_test),
    Matrix{dtype}(y_test)
end

"""
    snotel_metadata(; fields)

Return a database of snotel station metadata for usage in dataset creation.

# Arguments
- `fields::Vector{<:AbstractString}`: optional list of specific metadata fields to extract. Default is
`["stationID", "name", "state.code", "elevation", "latitude", "longitude"]`.
"""
function snotel_metadata(;
    fields::Vector{<:AbstractString} = [
        "stationId",
        "name",
        "state.code",
        "elevation",
        "latitude,longitude",
    ],
)
    link = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultipleStationReport,metric/daily/start_of_period/network=%22SNTL%22%7Cname/0,0/"
    for field in fields
        link = link * field * ","
    end
    data = df_from_url(link)
    data[!, :Elevation] .= data[!, :Elevation] * 12 * 2.54 / 100
    return data
end


"""
    fix_temp(temps; alaska)

Apply temperature bias corrections to raw SNOTEL temperature data, using temperature correction
formula developed by Atwood et. al. (2023)

# Arguments
- `temps::Vector{<:Real}`: the temperature data to fix.
- `alaska::Bool`: Specifies whether to apply the Alaska fix (as opposed to the CONUS fix). Default
 is `false`.
"""
function fix_temp(temps::Vector{<:Union{Real, Missing}}; alaska::Bool = false)
    tc = alaska ? (temps .+ 67.296) ./ 197.278 : (temps .+ 65.929) ./ 194.45
    c1 = 610558.226380138
    c2 = -2056177.65461394
    c3 = 2937046.42906361
    c4 = -2319657.12916417
    c5 = 1111854.33825836
    c6 = -337069.883250001
    c7 = 66105.7015922199
    c8 = -8386.78320604513
    c9 = 824.818021779729
    c0 = -86.7321006757439

    new_t =
        c1 .* (tc .^ 9) .+ c2 .* (tc .^ 8) .+ c3 .* (tc .^ 7) .+
        c4 .* (tc .^ 6) .+ c5 .* (tc .^ 5) .+ c6 .* (tc .^ 4) .+
        c7 .* (tc .^ 3) .+ c8 .* (tc .^ 2) .+ c9 .* tc .+ c0

    return new_t
end

"""
    serreze_qc(input, id, state)

Apply basic quality control measures to the daily data, using procedures
as developed by Serreze et. al. (1999), as well as an additional step to the
second filter, as applied by Yan et. al. (2018).
This function can only be applied over a valid SNOTEL site.

**Note: requires column names of at least "date", "SWE", "precip", "air_temp_avg" in the data, with data in sequential order.

# Arguments
- `input::DataFrame`: the data over which to apply quality control.
- `id::Int`: The site id of the SNOTEL site
- `state::AbstractString`: The state code (abbreviation) of the chosen SNOTEL site.
"""
function serreze_qc(input::DataFrame, id::Int, state::AbstractString)
    data = deepcopy(input)
    allowmissing!(data)
    years = sort(unique(Dates.year.(data[!, :date])))
    if Dates.month(data[end, :date]) >= 10
        push!(years, years[end] + 1)
    end
    for year in years
        oct1 = Date(string(year - 1) * "-10-01")
        oct15 = Date(string(year - 1) * "-10-15")
        allyear = oct1 .<= data[!, :date] .< (oct1 + Dates.Year(1))
        idx = oct1 .<= data[!, :date] .<= oct15
        testval = data[data[!, :date] .== oct1, :precip]
        case1 = sum(idx) != 15
        case2 = sum(ismissing.(data[idx, :SWE])) == 15
        case3 = sum(ismissing.(data[idx, :precip])) == 15
        case4 =
            (isempty(testval) | ismissing(testval)) ? true : (testval[1] > 5)
        clear_year = case1 | case2 | case3 | case4
        if clear_year
            data[allyear, [:SWE, :precip]] .= missing
            continue
        end
        precips = data[allyear, :precip]
        for i in 2:(length(precips) - 1)
            if sum(ismissing.(precips[(i - 1):(i + 1)])) == 0
                if 0 <= precips[i + 1] - precips[i - 1] <= 0.5
                    if (precips[i] < precips[i - 1]) |
                       (precips[i] > precips[i + 1])
                        precips[i] = (precips[i - 1] + precips[i + 1]) / 2.0
                    end
                end
            end
        end
        data[allyear, :precip] = precips
    end

    maxmin = get_data(id, state, ["TMIN", "TMAX"], metric = true)
    DataFrames.rename!(maxmin, [:date, :tmin, :tmax])
    overlap = intersect(maxmin[!, :date], data[!, :date])
    maxmin = maxmin[in.(maxmin[!, :date], [overlap]), :]
    data = data[in.(data[!, :date], [overlap]), :]
    data[!, :tmin] = maxmin[!, :tmin]
    data[!, :tmax] = maxmin[!, :tmax]
    flags =
        ismissing.(
            data[1:(end - 1), [:SWE, :precip, :air_temp_avg, :tmin, :tmax]],
        )
    flags[!, :date] = data[1:(end - 1), :date]

    dswes = data[2:end, :SWE] - data[1:(end - 1), :SWE]
    dprecips = data[2:end, :precip] - data[1:(end - 1), :precip]

    case1 = Vector{Union{Missing, Bool}}(abs.(dswes) .> 10)
    case1[ismissing.(case1)] .= true
    case2 = Vector{Union{Missing, Bool}}(
        (dswes[1:(end - 1)] .> 2.5) .& (dswes[2:end] .< -2.5),
    )
    case2 = case2 .| ((dswes[1:(end - 1)] .< -2.5) .& (dswes[2:end] .> 2.5))
    case2[ismissing.(case2)] .= true
    push!(case2, ismissing(dswes[end]))
    flags[(case1 .| case2), :SWE] .= true
    case3 = Vector{Union{Missing, Bool}}(-0.5 .< dprecips .< 10)
    case3[ismissing.(case3)] .= false
    flags[(!).(case3), :precip] .= true

    for var in [:air_temp_avg, :tmin, :tmax]
        temps = data[!, var]
        case4 = Vector{Union{Missing, Bool}}(abs.(temps) .> 40)
        case4[ismissing.(case4)] .= true
        case5 =
            Vector{Union{Missing, Bool}}(temps[2:end] .== temps[1:(end - 1)])
        case5[ismissing.(case5)] .= true
        flags[(case4[1:(end - 1)] .| case5), var] .= true
    end
    data = data[1:(end - 1), :]

    pos_dswe = (dswes .> 0)
    pos_dswe[ismissing.(pos_dswe)] .= false
    pos_dswe = Vector{Bool}(pos_dswe)
    zscores_pos, flags_pos_swe = zscore_filter(
        sqrt.(dswes[pos_dswe]),
        data[pos_dswe, :date],
        flags[pos_dswe, :SWE],
    )
    flags[pos_dswe, :SWE] .= (flags[pos_dswe, :SWE] .| flags_pos_swe)
    zscores_pos_swe = zeros(length(dswes))
    zscores_pos_swe[pos_dswe] .= zscores_pos

    ssqrt(x) = sqrt.(abs.(x)) .* sign.(x)
    neg_dswe = (dswes .< 0)
    neg_dswe[ismissing.(neg_dswe)] .= false
    neg_dswe = Vector{Bool}(neg_dswe)
    zscores_neg, flags_neg_swe = zscore_filter(
        ssqrt.(dswes[neg_dswe]),
        data[neg_dswe, :date],
        flags[neg_dswe, :SWE],
    )
    flags[neg_dswe, :SWE] .= (flags[neg_dswe, :SWE] .| flags_neg_swe)
    zscores_neg_swe = (zeros(length(dswes)))
    zscores_neg_swe[neg_dswe] .= zscores_neg

    nonzero_prec = (dprecips .!= 0)
    nonzero_prec[ismissing.(nonzero_prec)] .= false
    nonzero_prec = Vector{Bool}(nonzero_prec)
    zscores_nonzero_prec, flags_nonzero_prec = zscore_filter(
        ssqrt.(dprecips[nonzero_prec]),
        data[nonzero_prec, :date],
        flags[nonzero_prec, :precip],
    )
    flags[nonzero_prec, :precip] .=
        (flags[nonzero_prec, :precip] .| flags_nonzero_prec)
    zscores_prec = zeros(length(dprecips))
    zscores_prec[nonzero_prec] .= zscores_nonzero_prec

    _, flags_maxt =
        zscore_filter(data[!, :tmax], data[!, :date], flags[!, :tmax], z = 3)
    flags[flags_maxt, :tmax] .= true
    _, flags_mint =
        zscore_filter(data[!, :tmin], data[!, :date], flags[!, :tmin], z = 3)
    flags[flags_mint, :tmin] .= true
    zscores_avgt, flags_avgt = zscore_filter(
        data[!, :air_temp_avg],
        data[!, :date],
        flags[!, :air_temp_avg],
        z = 3,
    )
    flags[flags_avgt, :air_temp_avg] .= true

    unflag1 = (zscores_pos_swe .> 5) .& (zscores_prec .> 5)
    unflag2 = (zscores_pos_swe .> 5) .& (zscores_prec .> 3)
    unflag3 = (zscores_pos_swe .> 3) .& (zscores_prec .> 5)
    unflag4 = (zscores_neg_swe .< -5) .& (zscores_avgt .> 3)
    flags[unflag1 .| unflag2 .| unflag3, :SWE] .= false
    flags[unflag1 .| unflag2 .| unflag3, :precip] .= false
    flags[unflag4, :SWE] .= false
    flags[unflag4, :air_temp_avg] .= false

    for year in years
        oct1 = Date(string(year - 1) * "-10-01")
        yearend = oct1 + Dates.Year(1)
        yearflags = flags[oct1 .<= flags[!, :date] .< yearend, :]
        sweflag1 = findfirst(yearflags[!, :SWE])
        sweflag1 = isnothing(sweflag1) ? Inf : sweflag1
        preflag1 = findfirst(yearflags[!, :precip])
        preflag1 = isnothing(preflag1) ? Inf : preflag1
        idx = minimum([sweflag1, preflag1])
        if idx == Inf
            continue
        end
        idx = Int(idx)
        yeardate = yearflags[idx, :date]
        flags[yeardate .<= flags[!, :date] .< yearend, [:SWE, :precip]] .= true
    end

    data[flags[!, :SWE], :SWE] .= missing
    data[flags[!, :precip], :precip] .= missing
    temp = flags[!, :air_temp_avg] .| flags[!, :tmax] .| flags[!, :tmin]
    data[temp, :air_temp_avg] .= missing

    return data[!, Not([:tmin, :tmax])]
end

"""
    zscore_filter(data, dates, flags; z)

Apply filtering by zscore according to paramter z. Used by `serreze_qc()`.

# Arguments
- `data::Vector{<:Real}`: the numbers data to flag.
- `dates::Vector{Date}`: the dates of the data in order to block by month
- `z::Int`: Specifies the threshold of zscore to use in absolute value. Default is 5.
"""
function zscore_filter(data::Vector, dates::Vector, flags::Vector; z::Int = 5)
    monthdays = Dates.monthday.(dates)
    months = Dates.month.(dates)
    zscores = Vector{Union{Missing, Real}}(zeros(length(data)))
    zflags = Vector{Union{Missing, Bool}}(falses(length(data)))
    for month in 1:12
        gooddays = [(month - 1, 15)] .<= monthdays .<= [(month + 1, 15)]
        if month == 1
            gooddays = (monthdays .>= [(12, 15)]) .| (monthdays .>= [(2, 15)])
        elseif month == 12
            gooddays = (monthdays .>= [(11, 15)]) .| (monthdays .<= [(1, 15)])
        end
        usedays = gooddays .& ((!).(flags))
        if sum(usedays) < 30
            continue
        end
        μ = mean(skipmissing(data[usedays]))
        σ = std(skipmissing(data[usedays]))
        zscores[months .== month] = (data[months .== month] .- μ) / σ
        zflags[months .== month] = abs.(zscores[months .== month]) .> z
    end
    zflags[ismissing.(zflags)] .= true
    return zscores, Vector{Bool}(zflags)
end

"""
    yan_qc(input)

Apply quality control for inconsistent water years as proposed by Yan et. al. (2018).
This filter was designed for SNOTEl data but could be applied over any snow series data.

**Note: requires column names of at least "date", "SWE", "precip" in the data

# Arguments
- `input::DataFrame`: the site data over which to apply quality control.
"""
function yan_qc(input::DataFrame)
    data = deepcopy(input)
    allowmissing!(data)
    years = sort(unique(Dates.year.(data[!, :date])))
    if Dates.month(data[end, :date]) >= 10
        push!(years, years[end] + 1)
    end
    totdiff = 0.0
    n_inconsistent = 0.0
    numlost = 0
    for year in years
        oct1 = Date(string(year - 1) * "-10-01")
        yearend = oct1 + Dates.Year(1)
        allyear = oct1 .<= data[!, :date] .< yearend
        if sum(allyear) == 0
            continue
        elseif sum(ismissing, data[allyear, :SWE]) == sum(allyear)
            continue
        end
        (maxswe, maxidx) = findmax(skipmissing(data[allyear, :SWE]))
        sweprec = data[allyear, :precip][maxidx]
        if ismissing(sweprec)
            data[allyear, [:SWE, :precip]] .= missing
            continue
        end
        if maxswe > 1.05 * sweprec
            data[allyear, [:SWE, :precip]] .= missing
            totdiff += maxswe / sweprec
            n_inconsistent += 1
            numlost += sum(allyear)
        end
    end
    avgdiff = totdiff / n_inconsistent
    return data, avgdiff
end

"""
    livneh_bc(input)

Apply undercatch bias control for water years in the manner outlined by Livneh (2014),
as is done by Yan et. al. (2019).
This filter was originally applied to SNOTEL data but can be applied over any site, though the paper
recommends to avoid applying to maritime climate with greater occurrence of storms with mixed precipitation.

**Note: requires column names of at least "date", "SWE", "precip" in the data, ordered sequentially.

# Arguments
- `input::DataFrame`: the site data over which to apply quality control.
"""
function livneh_bc(input::DataFrame)
    data = deepcopy(input)
    years = sort(unique(Dates.year.(data[!, :date])))
    if Dates.month(data[end, :date]) >= 10
        push!(years, years[end] + 1)
    end

    for year in years
        oct1 = Date(string(year - 1) * "-10-01")
        yearend = Date(string(year) * "-10-01")
        precipdates = oct1 .<= data[!, :date] .< yearend
        yeardata = data[precipdates, [:precip, :SWE]]
        valid_idx = (!).(ismissing.(yeardata[!, :precip]))
        valid_idx = valid_idx .& (!).(ismissing.(yeardata[!, :SWE]))
        if sum(valid_idx) <= 2
            continue
        end
        accum_precip = yeardata[valid_idx, :precip]
        swe = yeardata[valid_idx, :SWE]
        precips = accum_precip[2:end] .- accum_precip[1:(end - 1)]
        Δswes = swe[2:end] .- swe[1:(end - 1)]

        idx_start = -6
        idx_end = 1
        while true
            idx_start += 7
            idx_end += 7
            temp = minimum([idx_end, length(accum_precip)])
            #=
            The following variables are from those in the paper:
            . precip_tot_7day = P7n in paper
            . swe_tot_7day = SWE7n in paper
            . precips_7day = P in paper
            . dswes_7day = dSWE in paper
            =#
            precip_tot_7day = accum_precip[temp] - accum_precip[idx_start]
            swe_tot_7day = swe[temp] - swe[idx_start]
            precips_7day = precips[idx_start:(temp - 1)]
            dswes_7day = Δswes[idx_start:(temp - 1)]
            if ismissing(precip_tot_7day) | ismissing(swe_tot_7day)
                (idx_end <= length(accum_precip)) || break
                continue
            end
            if swe_tot_7day > precip_tot_7day
                precips_7day[dswes_7day .>= 0] .= dswes_7day[dswes_7day .>= 0]
                precips_7day[dswes_7day .< 0] .=
                    [maximum([x, 0.0]) for x in precips_7day[dswes_7day .< 0]]
            else
                temp2 = (precips_7day .== 0)
                precips_7day =
                    maximum(cat(precips_7day, dswes_7day, dims = 2), dims = 2)
                precips_7day[temp2] .= 0.0
            end
            precips[idx_start:(temp - 1)] = precips_7day
            (idx_end <= length(accum_precip)) || break
        end
        pushfirst!(precips, accum_precip[1])
        new_accum_precip = cumsum(precips)
        yeardata[valid_idx, :precip] .= new_accum_precip
        data[precipdates, :precip] .= yeardata[!, :precip]
    end

    return data
end

"""
    z_filter(data; gap_days, thresh, rut_thresh, useswe)

Apply custom quality control measures to hourly depth data from a SNOTEL site.
Used by `z_qc()`.

**Note: requires column names of at least "date", "z", in the data. Also "SWE" if useswe is true. Dates must be in sequential order.

# Arguments
- `data::DataFrame`: the site data over which to apply quality control.
- `gap_days::Int`: the number of missing DAYS before all data is zeroed for that water year. Default is 30.
- `thresh::Real`: The threshold for allowable depth change in one sequential step, in inches. Default is 20 in.
- `rut_thresh::Real`: The amount of sequential steps in a rut before all data is zeroed for that water year. Default is 20.
- `useswe::Bool`: Use SWE data to apply an additional quality check based on possible density values. Default is true.
"""
function z_filter(
    data::DataFrame;
    gap_days::Int = 30,
    thresh::Real = 20,
    rut_thresh::Int = 20,
    useswe::Bool = true,
)
    flags = ismissing.(data[!, :z])

    nomiss = (!).(ismissing.(data[!, :z]))
    zs = data[nomiss, :z]
    ts = data[nomiss, :date]
    zflags = falses(size(zs))

    in_rut = false
    new_year = true
    wait_for_zero = false
    n_rut = 0
    spring_thaw = (Dates.month(ts[1]) in 4:8)
    for i in 2:(size(zs)[1] - 1)
        z = zs[i]
        dzp = zs[i + 1] - zs[i]
        dzm = zs[i] - zs[i - 1]
        dtp = ts[i + 1] - ts[i]
        dtp = Dates.value(Dates.Hour(dtp))

        if in_rut
            n_rut += 1
            if n_rut > rut_thresh
                wait_for_zero = true
            end
        end

        if wait_for_zero
            if z > 1
                zflags[i] = true
            else
                wait_for_zero = false
                in_rut = false
                n_rut = 0
                new_year = true
            end
            continue
        end

        if dtp > 24 * gap_days
            if in_rut
                zflags[i] = true
                in_rut = false
                n_rut = 0
                wait_for_zero = true
            end
            continue
        end

        if (Dates.month(ts[i]) == 10) & (!).(new_year)
            in_rut == false
            n_rut = 0
            new_year = true
            spring_thaw = false
        elseif (Dates.month(ts[i]) != 10) & new_year
            new_year = false
        end

        #handling summer issues
        if Dates.month(ts[i]) in 4:8
            if (!).(spring_thaw) & (z <= 1)
                spring_thaw = true
            elseif spring_thaw & (z != 0)
                zflags[i] = true
            end
        else
            if spring_thaw & (dzp > 0)
                spring_thaw = false
            end
        end

        if (dzm >= thresh) & (dzp <= -thresh)
            zflags[i] = true
        elseif (dzm <= -thresh) & (dzp >= thresh)
            zflags[i] = true
        end

        if (dzm >= thresh) & (abs(dzp) < thresh)
            in_rut = true
            zflags[i] = true
        elseif in_rut & (abs(dzp) < thresh)
            zflags[i] = true
        elseif in_rut & (abs(dzp) >= thresh)
            zflags[i] = true
            in_rut = false
            n_rut = 0
        end
    end
    zflags[end] = in_rut
    flags[nomiss] .= zflags

    case0 = data[!, :z] .> 175
    case0[ismissing.(case0)] .= true
    case0 = Vector{Bool}(case0)
    flags = flags .| case0

    if useswe
        dens = data[!, :z] ./ data[!, :SWE]
        case1 = dens .< 1.0
        case1 = case1 .| (dens .>= 50)
        case1[ismissing.(dens)] .= true
        flags = flags .| case1
    end

    return flags
end

"""
    z_qc(input; thresh, gd, zthresh, rthresh, useswe)

Quality control measures for hourly depth data in the paper. Iteratively applies `z_filter()` until
the number of flagged values increases by strictly less than `thresh`.

**Note: requires column names of at least "date", "z", in the data. Also "SWE" if `useswe` is true.

# Arguments
- `input::DataFrame`: the site data over which to apply quality control.
- `thresh::Int`: The permissible change in flagged values to stop the iterative application. Default is 1.
- `gd::Int`: the number of missing DAYS before all data is zeroed for that water year. Default is 30.
- `zthresh::Real`: The threshold for allowable depth change in one sequential step, in inches. Default is 20 in.
- `rthresh::Real`: The amount of sequential steps in a rut before all data is zeroed for that water year. Default is 20.
- `useswe::Bool`: Use SWE data to apply an additional quality check based on possible density values. Default is true.
"""
function z_qc(
    input;
    thresh::Int = 1,
    gd::Int = 30,
    zthresh::Real = 20,
    rthresh::Int = 20,
    useswe::Bool = true,
)
    data = deepcopy(input)
    flags = ismissing.(data[!, :z])
    numflags = count(ismissing, data[!, :z])
    n = 0
    while true
        n += 1
        zflags = z_filter(
            data,
            gap_days = gd,
            thresh = zthresh,
            rut_thresh = rthresh,
            useswe = useswe,
        )
        data[zflags, :z] .= missing
        nflags = sum(zflags)
        if nflags - numflags < thresh
            flags[:] .= zflags
            break
        end
        numflags += nflags - numflags
    end
    return flags
end

"""
    manual_filter(data)

The list of manually determined suspect data periods in the raw SNOTEL data, by site. Utilized in the paper.
"`s`" demarcation within the code refers to solar data, "`h`" for relative humidity data, "`w`" for wind data.
"`f`" demarcation of a period refers to periods that are suspect in multiple data streams.

**Note: requires column names of at least "date". 

# Arguments
- `data::DataFrame`: the site data over which to apply the filter.
"""
function manual_filter(data::DataFrame)
    sflags = falses(size(data)[1])
    hflags = falses(size(data)[1])
    wflags = falses(size(data)[1])
    if data[1, :id] .== 306
        s1 = Date("2010-05-01") .<= data[!, :date] .< Date("2012-09-26")
        s2 = data[!, :date] .>= Date("2019-08-28")
        h1 = Date("2007-10-01") .<= data[!, :date] .< Date("2008-08-01")
        h2 = data[!, :date] .>= Date("2020-01-01")
        sflags[s1 .| s2] .= true
        hflags[h1 .| h2] .= true

    elseif data[1, :id] .== 367
        s1 = data[!, :date] .< Date("2022-09-01")
        s2 = data[!, :date] .>= Date("2023-04-20")
        sflags[s1 .| s2] .= true

    elseif data[1, :id] .== 457
        s1 = Date("2005-08-01") .<= data[!, :date] .< Date("2022-07-01")
        sflags[(!).(s1)] .= true

    elseif data[1, :id] .== 482
        s1 = data[!, :date] .>= Date("2013-06-01")
        sflags[s1] .= true

    elseif data[1, :id] in [491, 1168]
        s1 = data[!, :date] .>= Date("2023-06-01")
        sflags[s1] .= true

    elseif data[1, :id] in [515, 532]
        s1 = data[!, :date] .>= Date("2022-06-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 571
        s1 = Date("2004-08-01") .<= data[!, :date] .< Date("2004-11-05")
        s2 = data[!, :date] .>= Date("2008-01-01")
        h1 = data[!, :date] .< Date("2004-11-01")
        h2 = data[!, :date] .>= Date("2006-07-27")
        sflags[s1 .| s2] .= true
        hflags[h1 .| h2] .= true

    elseif data[1, :id] .== 599
        s1 = data[!, :date] .>= Date("2023-05-01")
        s2 = data[!, :date] .== DateTime("2022-11-03T19:00:00")
        sflags[s1 .| s2] .= true

    elseif data[1, :id] .== 613
        s1 = Date("2015-06-01") .<= data[!, :date] .< Date("2019-10-01")
        s2 = data[!, :date] .< Date("2006-08-03")
        w1 = Date("2018-06-01") .<= data[!, :date] .< Date("2019-08-01")
        h1 = data[!, :date] .< Date("2007-09-01")
        sflags[s1 .| s2] .= true
        hflags[h1] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 665
        s1 = data[!, :date] .> Date("2022-06-01")
        w1 = data[!, :date] .< Date("2013-06-01")
        sflags[s1] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 708
        s1 = Date("2021-06-01") .<= data[!, :date] .< Date("2023-06-01")
        h1 = Date("2017-06-01") .<= data[!, :date] .< Date("2018-06-01")
        h2 = Date("2021-12-01") .<= data[!, :date] .< Date("2022-08-01")
        w1 = Date("2016-08-27") .<= data[!, :date] .< Date("2018-06-01")
        w2 = Date("2021-12-12") .<= data[!, :date] .< Date("2022-07-01")
        sflags[s1] .= true
        hflags[h1 .| h2] .= true
        wflags[w1 .| w2] .= true

    elseif data[1, :id] .== 715
        s1 = data[!, :date] .< Date("2019-10-01")
        s2 = data[!, :date] .>= Date("2023-05-01")
        s3 = Date("2021-02-01") .<= data[!, :date] .< Date("2021-06-01")
        h1 = data[!, :date] .<= Date("2010-01-01")
        w1 = Date("2015-08-01") .<= data[!, :date] .< Date("2017-10-01")
        sflags[s1 .| s2 .| s3] .= true
        hflags[h1] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 744
        s1 = data[!, :date] .< Date("2013-06-01")
        s2 = data[!, :date] .>= Date("2019-06-01")
        h1 = data[!, :date] .<= Date("2012-06-01")
        sflags[s1 .| s2] .= true
        hflags[h1] .= true

    elseif data[1, :id] .== 832
        w1 = Date("2016-08-17") .<= data[!, :date] .< Date("2018-08-01")
        wflags[w1] .= true

    elseif data[1, :id] .== 845
        s1 = data[!, :date] .>= Date("2021-06-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 854
        f1 = data[!, :date] .<= Date("2018-01-01")
        w1 = data[!, :date] .<= Date("2011-10-01")
        sflags[f1] .= true
        hflags[f1] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 921
        s1 = data[!, :date] .< Date("2018-06-01")
        s2 = data[!, :date] .> Date("2023-06-01")
        s3 = Date("2021-01-01") .<= data[!, :date] .< Date("2021-09-01")
        h1 = Date("2014-06-01") .<= data[!, :date] .< Date("2014-10-01")
        w1 = data[!, :date] .< Date("2012-01-01")
        sflags[s1 .| s2 .| s3] .= true
        hflags[h1] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 922
        s1 = Date("2017-06-01") .<= data[!, :date] .< Date("2018-06-01")
        f1 = data[!, :date] .< Date("2014-11-01")
        h1 = data[!, :date] .> Date("2016-06-01")
        sflags[s1 .| f1] .= true
        hflags[h1 .| f1] .= true

    elseif data[1, :id] .== 927
        s1 = Date("2016-06-01") .<= data[!, :date] .< Date("2020-07-01")
        s2 = data[!, :date] .< Date("2013-06-01")
        h1 = Date("2016-01-01") .<= data[!, :date] .< Date("2018-09-01")
        h2 = data[!, :date] .< Date("2014-07-01")
        w1 = data[!, :date] .>= Date("2022-06-01")
        sflags[s1 .| s2] .= true
        hflags[h1 .| h2] .= true
        wflags[w1] .= true

    elseif data[1, :id] .== 969
        f1 = Date("2008-06-01") .<= data[!, :date] .< Date("2009-06-01")
        s1 = Date("2016-05-01") .<= data[!, :date] .< Date("2016-09-01")
        s2 = data[!, :date] .< Date("2006-10-01")
        s3 = data[!, :date] .>= Date("2022-06-01")
        s4 = Date("2016-04-01") .<= data[!, :date] .< Date("2016-08-01")
        w1 = Date("2016-05-01") .<= data[!, :date] .< Date("2017-10-01")
        w2 = Date("2020-08-01") .<= data[!, :date] .< Date("2021-02-01")
        sflags[s1 .| s2 .| s3 .| s4 .| f1] .= true
        hflags[f1] .= true
        wflags[w1 .| w2] .= true

    elseif data[1, :id] .== 974
        w1 = data[!, :date] .< Date("2010-09-01")
        wflags[w1] .= true

    elseif data[1, :id] .== 978
        s1 = data[!, :date] .>= Date("2018-08-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 1030
        h1 = Date("2016-09-01") .<= data[!, :date] .< Date("2017-09-01")
        w1 = Date("2012-01-01") .<= data[!, :date] .< Date("2013-06-01")
        w2 = Date("2021-06-01") .<= data[!, :date] .< Date("2022-10-01")
        s1 = data[!, :date] .< Date("2005-01-01")
        s2 = data[!, :date] .>= Date("2022-06-01")
        w3 = data[!, :date] .< Date("2006-06-01")
        sflags[s1 .| s2] .= true
        hflags[h1] .= true
        wflags[w1 .| w2 .| w3] .= true

    elseif data[1, :id] .== 1053
        s1 = data[!, :date] .>= Date("2020-06-01")
        s2 = Date("2006-06-01") .<= data[!, :date] .< Date("2007-09-01")
        sflags[s1 .| s2] .= true

    elseif data[1, :id] .== 1083
        s1 = data[!, :date] .>= Date("2023-06-01")
        f1 = Date("2012-10-25") .<= data[!, :date] .< Date("2012-11-24")
        sflags[s1 .| f1] .= true
        hflags[f1] .= true

    elseif data[1, :id] .== 1105
        s1 = data[!, :date] .>= Date("2012-06-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 1170
        f1 = data[!, :date] .>= Date("2023-05-01")
        sflags[f1] .= true
        wflags[f1] .= true

    elseif data[1, :id] .== 1254
        f1 = data[!, :date] .>= Date("2023-05-01")
        sflags[f1] .= true
        hflags[f1] .= true

    elseif data[1, :id] .== 1286
        s1 = Date("2021-06-01") .<= data[!, :date] .< Date("2021-10-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 963
        s1 = data[!, :date] .< Date("2007-06-01")
        sflags[s1] .= true

    elseif data[1, :id] .== 1035
        s1 = Date("2009-06-01") .<= data[!, :date] .< Date("2013-06-01")
        sflags[(!).(s1)] .= true

    elseif data[1, :id] .== 1092
        s1 = data[!, :date] .> Date("2018-06-01")
        sflags[s1] .= true
    end
    return sflags, hflags, wflags
end

"""
    impute_data(data, var; t1, t2, dt)

Apply a generalizable imputation step to the data of variable `var` in `data`.
For a timestep of `dt` in the data, fills gaps of length `t1` steps via linear interpolation,
and for gaps longer than `t1` steps but less than `t2` steps, fills gaps with an appropriately determined profile.
Setting `t1` or `t2` as -1 or `-Inf` will not result in that interpolation type being used.
See `make_profile()` for the profile types of `dt = Day(1)` and `dt = Hour(1)`. No others are encoded.

**Note: requires column names of at least "date", and `var` in the data, with dates in sequential order.

# Arguments
- `data::DataFrame`: the site data over which to apply interpolation.
- `var::Symbol`: the specific variable stream over which to fill gaps. Default is `:sol_rad_avg`.
- `t1::Int`: Gap size (units of `dt`) over which to apply linear interpolation. Use -1 for no linear interpolation, default is 6.
- `t2::Int`: Gap size (units of `dt`) over which to apply profile interpolation. Use -1 for no profile interpolation, default is 24.
- `dt::Period`: the length of timestep in the data. Default is `Hour(1)`.
"""
function impute_data(
    data,
    var::Symbol = :sol_rad_avg;
    t1::Int = 6,
    t2::Int = 24,
    dt::Period = Hour(1),
)
    newvar = Vector{Union{Missing, Float32}}(deepcopy(data[!, var]))
    prof = make_profile(data, var, dt)
    bwk = Int.(ceil.((Dates.week.(data[!, :date])) ./ 2))
    bwk[bwk .== 27] .= 26
    s = findfirst((!).(ismissing.(data[!, var])))
    e = findlast((!).(ismissing.(data[!, var])))
    idx_miss = s
    while (idx_miss < e)
        idx_miss = findnext(ismissing.(newvar), idx_miss)
        if isnothing(idx_miss)
            break
        elseif idx_miss > e
            break
        end
        idx_next = findnext((!).(ismissing.(newvar)), idx_miss)
        if isnothing(idx_next)
            break
        elseif idx_next > e
            break
        end
        Δt = (data[idx_next, :date] - data[idx_miss, :date]) / dt
        if Δt <= t1
            interp = Vector{Float32}(
                collect(
                    range(
                        newvar[idx_miss - 1],
                        newvar[idx_next],
                        idx_next - idx_miss + 2,
                    ),
                ),
            )
            newvar[(idx_miss - 1):idx_next] .= interp
        elseif Δt <= t2
            for k in idx_miss:(idx_next - 1)
                a = bwk[k]
                if dt == Hour(1)
                    b = Dates.hour(data[k, :date]) + 1
                    newvar[k] = prof[a, b]
                elseif dt == Day(1)
                    newvar[k] = prof[a]
                end
            end
        end
        idx_miss = idx_next
    end
    return newvar
end

"""
    make_profile(data, var; dt)

Create the profile for interpolation in `impute_data()` of variable `var`.
For `dt = Hour(1)`, makes a 26x24 "hour per biweek" profile across a year of the variable.
For `dt = Day(1)`, makes a 26-element vector "biweek" profile across a year of the variable.
Value of the profile is determined by the mean in this code, using the custom `missingmean` function to handle missing values.
Profiles for using the mean or other `dt`'s must be created by extending/replacing this function.

**Note: requires column names of at least "date", and `var` in the data.

# Arguments
- `data::DataFrame`: the site data over which to apply interpolation.
- `var::Symbol`: the specific variable stream over which to fill gaps. Default is `:sol_rad_avg`.
- `dt::Period': the length of timestep in the data. Default is `Hour(1)`.
"""
function make_profile(data::DataFrame, var::Symbol, dt::Period = Hour(1))
    bwk = ceil.((Dates.week.(data[!, :date])) ./ 2)
    bwk[bwk .== 27.0] .= 26.0
    if dt == Hour(1)
        profile = Matrix{Union{Missing, Real}}(zeros((26, 24)))
        hr = Dates.hour.(data[!, :date]) .+ 1
        for i in 1:26
            for j in 1:24
                pdata = data[(hr .== i) .& (bwk .== j), var]
                profile[i, j] = missingmean(pdata)
            end
        end
        return profile
    elseif dt == Day(1)
        profile = Vector{Union{Missing, Real}}(zeros(26))
        for i in 1:26
            pdata = data[(bwk .== i), var]
            profile[i] = missingmean(pdata)
        end
        return profile
    end
end

"""
    d_impute(data, t1, t2)

Helper/wrapper function to interpolate over daily data, calling `impute_data()` on all paper variables with `dt = Day(1)`.
# Arguments
- `data::DataFrame`: the site data over which to apply interpolation.
- `t1::Int`: value of t1 (linear interpolation) to use for the variables. Set to 3.
- `t2::Int`: value of t2 (profile interpolation) to use for the variables. Set to -1 for no profile interpolation.
"""
function d_impute(data::DataFrame, t1::Int = 3, t2::Int = -1)
    newd = deepcopy(data)
    st = maximum(findfirst.(eachcol((!).(ismissing.(newd)))))
    ed = minimum(findlast.(eachcol((!).(ismissing.(newd)))))
    newd[!, :z] = convert(Vector{Union{Missing, Real}}, newd[!, :z])
    for var in [:z, :SWE, :precip]
        imp = impute_data(newd[st:ed, :], var, t1 = t1, t2 = -1, dt = Day(1))
        newd[st:ed, var] .= imp
    end
    for var in [:sol_rad_avg, :rel_hum_avg, :wind_speed_avg, :air_temp_avg]
        imp = impute_data(newd[st:ed, :], var, t1 = t1, t2 = t2, dt = Day(1))
        newd[st:ed, var] .= imp
    end
    return newd
end

"""
    qc_filter(data, var, t1, nw, verbose)

Assess datastream for suspect extreme values via median `x̄` and of interquartile range `IQR` of weekly maximums of `var` in `data`.
Flag all individual data satisfying `(x - x̄)/IQR > t1`, and flag entire blocks having over 5% flagged values, by
growing and checking "blocks" at a rate of nw entries per check.
Makes use of custom `missingmax` function to handle missing values.

This is called once inside the hourly quality control checks on wind_speed data, and once on solar data after collecting
hourly data into daily data.

# Arguments
- `data::DataFrame`: the site data over which to apply the check.
- `var::Symbol`: the specific variable stream over which to assess. Default is `:wind_speed_avg`.
- `t1::Int`: the score over which to flag a value (modified z-score, using median and interquartile range). Default is 6.
- `nw::Int`: the rate at which to grow/check blocks for sufficient count of missing values. Default is 72.
- `verbose::Bool`: optional helper arg to help with debugging and calibration. Prints statistics, requires site "id" to be a data column. Default is `false`.
"""
function qc_filter(
    data::DataFrame,
    var::Symbol = :wind_speed_avg;
    t1::Int = 6,
    nw::Int = 72,
    verbose::Bool = false,
)
    gdata = deepcopy(data[!, [:date, var]])
    flags = ismissing.(data[!, var])
    n1 = length(flags) - sum(flags)
    tchunk = Vector{Real}(Dates.year.(data[!, :date]))
    if var == :wind_speed_avg
        tchunk .+= Dates.week.(data[!, :date]) / 100.0
    end
    gdata[!, :tchunk] .= tchunk
    mmaxes =
        combine(groupby(gdata, :tchunk), var => missingmax, renamecols = false)
    mu = median(skipmissing(mmaxes[!, 2]))
    sp =
        quantile(skipmissing(mmaxes[!, 2]), 0.75) -
        quantile(skipmissing(mmaxes[!, 2]), 0.25) #iqr(skipmissing(mmaxes[!, 2]))
    scores = ((gdata[!, 2] .- mu)) ./ sp .> t1
    if verbose
        print(
            "X̄, IQR, N_VIOLATIVE: ",
            mu,
            ", ",
            sp,
            " ",
            sum(skipmissing(scores)),
            "\n",
        )
    end
    if sum(skipmissing(scores)) .< 24
        flags[(!).(ismissing.(scores))] .= skipmissing(scores)
        return flags
    else
        flags[findall(skipmissing(scores))] .= true
        idx = findfirst(skipmissing(scores))
        while idx < length(scores)
            idxe = minimum([idx + nw - 1, length(scores)])
            tsum = 1
            pre = 0
            while sum(skipmissing(scores[idx:idxe])) > tsum
                tsum = sum(skipmissing(scores[idx:idxe]))
                idxe = minimum([idxe + nw, length(scores)])
                pre += 1
            end
            if tsum > 0.05 * (idxe - idx + pre)
                idx = maximum([idx - pre, 1])
                flags[idx:idxe] .= true
            end
            newidx = findfirst(skipmissing(scores[idxe:end]))
            idx = isnothing(newidx) ? length(scores) : idxe - 1 + newidx
        end
    end
    n2 = length(flags) - sum(flags)
    if verbose
        print("N_LOST: ", n1 - n2, "\n")
    end
    return flags
end

"""
    bcqc_hourly(input, id, state)

Apply bias control and quality control to a SNOTEL hourly dataset, according to paper.

*Some extra steps are applied outside this function.

Steps are:

. apply `manual_filter()` to flag suspect data

. apply `qc_filter()` to wind data *(also applied to solar data with `t1 = 2`, after combining to daily scale with `hourly2daily()`)

. apply `z_qc()` to depth data

. apply `impute_data()` to solar, humidity, wind, and temperature data

# Arguments
- `input::DataFrame`: the hourly site data over which to apply quality control.
"""
function bcqc_hourly(input::DataFrame)
    data = deepcopy(input)
    sflags, hflags, wflags = manual_filter(data)
    wflags2 = qc_filter(data)
    z_flags = z_qc(data)
    allowmissing!(data)
    data[sflags, :sol_rad_avg] .= missing
    data[wflags .| wflags2, :wind_speed_avg] .= missing
    data[hflags, :rel_hum_avg] .= missing
    data[z_flags, :z] .= missing
    for var in [:sol_rad_avg, :rel_hum_avg, :wind_speed_avg, :air_temp_avg]
        imp = impute_data(data, var, t1 = 6, t2 = 24)
        data[!, var] .= imp
    end
    return data
end

"""
    bcqc_daily(input, id, state)

Apply bias control and quality control to a SNOTEL daily dataset, according to paper.

*Some extra steps are applied outside this function

Steps are:

. Apply `serreze_qc()` to data

. Apply `yan_qc()` to data

. Apply `fix_temp()` to applicable SNOTEL sites

. Apply `livneh_bc()` to data

. Apply a density check over z, SWE data.

(*`d_impute()` is called after this function)

# Arguments
- `input::DataFrame`: the site data over which to apply quality control.
- `id::Int`: The SNOTEL ID of the station over which the QC is applied.
- `state::Int`: The state code of the station over which the QC is applied.
"""
function bcqc_daily(input::DataFrame, id::Int, state::AbstractString)
    data = serreze_qc(input, id, state)
    is_alaska = (state == "AK")
    data, _ = yan_qc(data)
    if (!).(in(id, [1070, 1091, 1092])) #As of May 2024
        data[!, :air_temp_avg] .=
            fix_temp(data[!, :air_temp_avg], alaska = is_alaska)
    end
    data = livneh_bc(data)

    dens = data[!, :z] ./ data[!, :SWE]
    zflags = ismissing.(dens)
    case1 = dens .< 1.0
    case1 = case1 .| (dens .>= 50)
    case1[ismissing.(dens)] .= true
    zflags = zflags .| case1
    data[zflags, :z] .= missing
    return data
end

"""
    scrape_site_paper(id, state)

Easy wrapper to scrape site data in the exact same way as the paper.
For a given site, this could take up to a minute or two.

This contains/encodes all "exceptions" noted in the paper for ease of utility to a user,
though a user looking to investigate optional/additional features from this code is better
off implementing this code manually. Two special circumstances not explicitly visible in this code are
the `z_qc()` filter having a cutoff at z > 175 inches, and only three sites (1070, 1091, 1092)
not requiring temperature fixes (as of May 2024) within `bcqc_daily()`.

**Note: This function will likely not work or give unexpected results for a site not utilized in the paper.

# Arguments
- `id::Int`: the SNOTEL ID of the site to extract.
- `state::AbstractString`: the state code of the site to extract.
"""
function scrape_site_paper(id::Int, state::AbstractString)
    inch2meter = 0.0254
    kmphr2mps = 5.0 / 18.0
    metricfields = ["RHUMV", "SRADV", "WSPDV", "TOBS"]
    site = id
    start_date = "start"
    end_date = "end"

    #Handle exceptions:
    if (state == "AK") | (id == 2170)
        metricfields[1] = "RHUM" #A different data tag was used at these sites
    end
    tlim = (state == "AK") ? -50.0 : -40.0 #handles a handful of individual spikes in some sites
    if id == 515
        end_date = "2023-06-02" #portal throws an error if trying to scrape until 2024, as of May 2024
    end

    #Set up the filter values based on the SNOTEL data sensors:
    filter_val = Dict{Symbol, Tuple{Real, Real}}(
        :SWE => (0.0, 250.0),
        :z => (0.0, 420.0),
        :precip => (0.0, 250.0),
        :rel_hum_avg => (10.0, 100.0),
        :sol_rad_avg => (0.0, 1500.0),
        :wind_speed_avg => (0.0, 216.0),
        :air_temp_avg => (tlim, 40.0),
    )

    #Set up the scaling dictionary used in scaling data:
    scales = Dict{Symbol, Real}(
        :SWE => inch2meter,
        :z => inch2meter,
        :precip => inch2meter,
        :rel_hum_avg => 0.01,
        :wind_speed_avg => kmphr2mps,
    )

    #Gather and handle the hourly data:
    hourly = apply_bounds(
        sitedata_hourly(
            site,
            state,
            metric_fields = metricfields,
            start = start_date,
            finish = end_date,
        ),
        filter_val,
    )
    hourly[!, :id] .= site
    hourly = bcqc_hourly(hourly)
    hourly_d = hourly2daily(hourly)
    allowmissing!(hourly_d)
    sflags = qc_filter(hourly_d, :sol_rad_avg, t1 = 2)
    hourly_d[sflags, :sol_rad_avg] .= missing

    #Gather daily data:
    daily = apply_bounds(
        sitedata_daily(site, state, start = start_date, finish = end_date),
        filter_val,
    )
    daily[!, :id] .= site

    #Combine daily/hourly data and handle exceptions:
    gap_daily = rectify_daily_hourly(daily, hourly_d)
    if site in [306, 978]
        gap_daily = rectify_daily_hourly(
            daily,
            hourly_d,
            prioritize_hour = [:rel_hum_avg],
            hourly_only = [:sol_rad_avg, :wind_speed_avg],
        )
    end
    if site in [1122]
        gap_daily = rectify_daily_hourly(
            daily,
            hourly_d,
            prioritize_hour = [:air_temp_avg],
        )
    end

    #Finish generating the site dataset:
    daily_qc = bcqc_daily(gap_daily, site, state)
    daily_qc = d_impute(daily_qc)
    daily_scaled = scale_cols(daily_qc, scales)
    daily_clean = daily_scaled[completecases(daily_scaled), :]
    daily_clean = makediffs(daily_clean, Day(1))
    good_vals = daily_clean[!, :dprecipdt] .>= 0.0
    daily_clean[(!).(good_vals), :dprecipdt] .= 0.0
    daily_clean = daily_clean[!, Not(:precip)]

    #The tutorial adds metadata to this dataset outside this function:
    return daily_clean
end

"""
    save_df(df, filepath; delim)

Utility function to save DataFrames to files. Custom delimiters can be used with the "delimiter" argument.

# Arguments
- `df::DataFrame`: the DataFrame to save.
- `filepath::AbstractString`: the desired place to save the DataFrame as a delimited file.
- `delim`: the character with which to separate values. Default is ",".
"""
function save_df(
    df::DataFrame,
    filepath::AbstractString;
    delim::AbstractString = ",",
)
    open(filepath, "w") do io
        println(io, join(names(df), delim))
        for row in eachrow(df)
            values = [ismissing(val) ? "" : string(val) for val in row]
            println(io, join(values, delim))
        end
    end
end

end
