module DataTools
using DataFrames, CSV, HTTP, Dates, StatsBase
export sitedata_daily,
    sitedata_hourly,
    apply_bounds,
    hourly2daily,
    rectify_daily_hourly,
    scale_cols,
    makediffs,
    rolldata,
    snowsplit,
    filter_phys!,
    prep_data,
    make_data,
    snotel_metadata

"""
    data_url(id, state, fields; metric, hourly, start_date, end_date)

Return the url to gather SNOTEL data from a particular station.

# Arguments
- `id::String`: the id code of the station to access.
- `state::String`: the state abbreviation of the station to access.
- `fields::Vector{String}`: a list of strings specifying the SNOTEL parameters to collect.
- `metric::Bool`: a boolean specifying whether to gather data in metric (true) or imperial (false) units. Default is false.
- `hourly::Bool`: a boolean specifying whether to scrape hourly (true) or daily (false) data. default is false.
- `start_date::String`:
- `end_date::Stirng`:
"""
function data_url(
    id::Int,
    state::AbstractString,
    fields::Vector{String};
    metric::Bool = false,
    hourly::Bool = false,
    start_date::String = "start",
    end_date::String = "end",
)
    start = (start_date == "start") ? "1850-01-01" : start_date
    fin = (end_date == "end") ? "2023-06-01" : end_date
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

Return the SNOTEL data from a station as a DataFrame.

# Arguments
- `id::String`: the id code of the station to access.
- `state::String`: the state abbreviation of the station to access.
- `fields::Vector{String}`: a list of strings specifying the SNOTEL parameters to collect.
- `metric::Bool = false`: a boolean specifying whether to gather data in metric (true) or imperial (false) units.
- `hourly::Bool = false`: a boolean specifying whether to scrape hourly (true) or daily (false) data.
- `start::String`: Optional string to specify starting date of data collection. Default is "start".
- `finish::String`: Optional string to specify ending date of data collection. Default is "end".
"""
function get_data(
    id::Int,
    state::AbstractString,
    fields::Vector{String};
    metric::Bool = false,
    hourly::Bool = false,
    start::String = "start",
    finish::String = "end",
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
    data = CSV.read(HTTP.get(url).body, DataFrame, comment = "#", delim = ",")
    return data
end


"""
    sitedata_daily(id, state; imp_fields, metric_fields, colnames, start, finish)

Return the daily SNOTEL data from a station as a DataFrame.

# Arguments
- `id::String`: the id code of the station to access.
- `state::String`: the state abbreviation of the station to access.
- `imp_fields::Vector{String}`: Parameters to return in imperial units. Default is ["WTEQ", "SNWD", "PREC"].
-  `metric_fields::Vector{String}`: Parameters to return in metric units. Default is ["RHUMV", "SRADV", "WSPDV", "TAVG"].
- `colnames::Vector{String}`: Optional column names to change header after scraping data. Default follows that of the paper,
which is ["date", "SWE", "z", "precip", "rel_hum_avg", "sol_rad_avg", "wind_speed_avg", "air_temp_avg"].
- `start::String`: Optional string to specify starting date of data collection. Default is "start".
- `finish::String`: Optional string to specify ending date of data collection. Default is "end".
"""
function sitedata_daily(
    id::Int,
    state::AbstractString;
    imp_fields::Vector{String} = ["WTEQ", "SNWD", "PREC"],
    metric_fields::Vector{String} = ["RHUMV", "SRADV", "WSPDV", "TAVG"],
    colnames::Vector{String} = [
        "date",
        "SWE",
        "z",
        "precip",
        "rel_hum_avg",
        "sol_rad_avg",
        "wind_speed_avg",
        "air_temp_avg",
    ],
    start::String = "start",
    finish::String = "end",
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

Return the hourly SNOTEL data from a station as a DataFrame.

# Arguments
- `id::String`: the id code of the station to access.
- `state::String`: the state abbreviation of the station to access.
- `imp_fields::Vector{String}`: Parameters to return in imperial units. Default is ["WTEQ", "SNWD", "PREC"].
-  `metric_fields::Vector{String}`: Parameters to return in metric units. Default is ["RHUMV", "SRADV", "WSPDV", "TAVG"].
- `colnames::Vector{String}`: Optional column names to change header after scraping data. Default follows that of the paper,
which is ["date", "SWE", "z", "precip", "rel_hum_avg", "sol_rad_avg", "wind_speed_avg", "air_temp_avg"].
- `start::String`: Optional string to specify starting date of data collection. Default is "start".
- `finish::String`: Optional string to specify ending date of data collection. Default is "end".
"""
function sitedata_hourly(
    id::Int,
    state::AbstractString;
    imp_fields::Vector{String} = ["WTEQ", "SNWD", "PREC"],
    metric_fields::Vector{String} = ["RHUMV", "SRADV", "WSPDV", "TOBS"],
    colnames::Vector{String} = [
        "date",
        "SWE",
        "z",
        "precip",
        "rel_hum_avg",
        "sol_rad_avg",
        "wind_speed_avg",
        "air_temp_avg",
    ],
    start::String = "start",
    finish::String = "end",
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
    hourly2daily(hourlydata)

Convert an hourly SNOTEL dataset to a daily SNOTEL dataset, averaging humidity, radiation, wind, and temperature
data, but maintaining start-of-window SWE, z, and precipitation data.
Makes use of custom "missingmean" and "missingfirst" functions to handle missing values when finding the mean and first non-missing
values, respectively, when working with hourly datasets with missing values.
**Note: does not extend to dataframes with fields beyond those extracted for the paper.

# Arguments
- `hourlydata::DataFrame`: the hourly data over which to average by-day.
"""
function hourly2daily(hourlydata::DataFrame)
    data = deepcopy(hourlydata)
    data[!, :date] .= Date.(DateTime.(data[!, :date]))
    missingmean(input) =
        isempty(collect(skipmissing(input))) ? missing :
        mean(skipmissing(input))
    missingfirst(input) =
        isempty(collect(skipmissing(input))) ? missing :
        first(skipmissing(input))
    dailydata = combine(
        groupby(data, :date),
        [:SWE, :z, :precip] .=> missingfirst,
        [:rel_hum_avg, :sol_rad_avg, :wind_speed_avg, :air_temp_avg] .=>
            missingmean,
        renamecols = false,
    )
    return dailydata
end

"""
    rectify_daily_hourly(daily_data, hourly_data)

Use one SNOTEL data stream (hourly data) to fill holes in another SNOTEL data stream (daily data).
**Note: requires time column to be named "date".

# Arguments
- `daily_data::DataFrame`: the main (daily) data over which to fill missing holes.
- `hourly_data::DataFrame`: the (converted to daily-format, see hourly2daily) hourly data over which to fill holes in the daily data.
"""
function rectify_daily_hourly(daily_data::DataFrame, hourly_data::DataFrame)
    vars = names(daily_data)
    combined = outerjoin(daily_data, hourly_data, on = :date, makeunique = true)
    for var in vars[2:end]
        combined[!, Symbol(var)] .=
            coalesce.(combined[!, Symbol(var)], combined[!, Symbol(var * "_1")])
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
(data[i+1]-data[i])/(time[i+1]-time[i]), but only in the case where time[i+1]-time[i] = Δt.
**Note: apply after scaling and getting rid of missing values. Assumes time column has name "date" and has Date/DateTime units.

# Arguments
- `data::DataFrame`: the data over which to apply differencing.
- `Δt::Period`: the amount of time representing one unit timestep in the data.
- `diffvars::Vector{Symbol}`: Columns to apply differencing to. Default is [:SWE, :z, :precip].
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

Convert a vector of vectors into a DataFrame, with specified column names.

# Arguments
- `stack::Vector{Vector{Any}}`: the data stack to convert.
- `colnames::Vector{String}`: The names to give the columns of the DataFrame.
"""
function stack2DF(stack, colnames)
    return DataFrame(mapreduce(permutedims, vcat, stack), colnames)
end

"""
    rolldata(data, Δt, N; takefirst)

Apply a moving average of N timesteps to all data, except for the variables
specified in "takefirst", for which the leading value is maintained.
**Note: assumes the time column is named "date" and has Date/DateTime units

# Arguments
- `data::DataFrame`: the data over which to apply averaging.
- `Δt::Period`: the amount of time representing one unit timestep in the data.
- `N::Int`: the number of intervals (timesteps) to include in the average
- `takefirst::Vector{Symbol}`: Columns to apply differencing to. Default is [:date, :SWE, :z, :precip].
"""
function rolldata(
    data,
    Δt::Period,
    N::Int;
    takefirst::Vector{Symbol} = [:date, :SWE, :z, :precip],
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
    return stack2DF(rolleddata, names(data))
end

"""
    snowsplit(air_temp, hum, precip)

Engineer total water content of precipitation into snow and rain portions,
accoridng to the paper outlined in https://www.nature.com/articles/s41467-018-03629-7.

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
    filter_phys!(data; eps)

Filter unphysical/undesirable data points out of the dataset, following physical choices outlined in the paper, such as:
- Removing rows were SWE, z, and precipitation are all less than some threshold
- Removing rows where SWE is zero and z is nonzero
- Removing rows where dz/dt is positive but precipitation is zero
- Removing rows where dSWE/dt > precipitation
- Removing rows where SWE < z
**Note: requires column names to match those of the paper for usage.

# Arguments
- `data::DataFrame`: the cleaned and processed data.
- `eps::Real`: a filtering threshold, set to 0.5 cm.
"""
function filter_phys!(data::DataFrame; eps::Real = 0.005)
    #these need to be modified according to paper commentary and correct units.
    zero_condition1 =
        (data[!, :SWE] .< eps) .&
        (data[!, :z] .< eps) .&
        (data[!, :dprecipdt] .< eps / 86400.0)
    data = data[Not(zero_condition1), :]
    zero_condition2 = (data[!, :z] .!= 0.0f0) .& (data[!, :SWE] .== 0.0f0)
    data = data[Not(zero_condition2), :]
    zero_condition3 = (data[!, :dzdt] .> 0.0) .& (data[!, :dprecipdt] .== 0.0)
    data = data[Not(zero_condition3), :]
    data = data[data[!, :dSWEdt] .<= data[!, :dprecipdt], :]
    data = data[data[!, :SWE] .<= data[!, :z], :]
end

"""
    prep_data(data; extract_vars, make_snow_split, physical_filter, eps)

Prepare a cleaned (scaled & gap-filled, potentially rolled) data stream or non-scraped data stream for model usage.
**Note: Requires column names to match those of the paper for usage.

# Arguments
- `data::DataFrame`: the cleaned and processed data.
- `extract_vars::Vector{Symbol}`: The list of columns to be used in the model. Default is the variables used in the paper.
- `make_snow_split::Bool`: Boolean indicating whether to split precipitation into water and snow. Default is true.
- `physical_filter::Bool`: Boolean indicating whether to filter data using filter_phys!. Default is true.
- `eps::Real`: a filtering threshold, set to 0.5 cm.
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
        :dzdt,
    ],
    make_snow_split::Bool = true,
    physical_filter::Bool = true,
    eps::Real = 0.005,
)
    if physical_filter
        filter_phys!(data, eps = eps)
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

Create training and testing matrices for the models from prepped input data, to be passed to a Flux DataLoader call for batching.

# Arguments
- `data::DataFrame`: the data to be used for training and validation data.
- `input_vars::Vector{Symbol}`: The set of features to be extracted as input data
- `target::Symbol`: The feature to be extracted as the target variable
- `target_scale::Real`: The scaling to apply to the output data
- `testidx': a boolean array to demarcate indices of testing/validation data, or "nothing". Default is "nothing"
- `dtype::Type`: The data type consistent with the model. Default is Float32.
"""
#use MLUtils instead? reason I didn't: allows the choice of test indices to be conditional instead of non-numerical,
# i.e. i can past the set of testindices to be "all data with this site id" or even "all data with z < 0.3" in order
# to analyze different test cases. MLUtils doesn't have this capability but does have interesting class-balancing
# capabilities and other stuff. All this function does is convert DataFrame data to a Matrix, otherwise, plus
# a wrapper around the arbitrary testindex option, so it's not saving/adding any calls to use MLUtils instead
# nor should it become obsolete or require updating.
#rename to split_data? or something to indicate no longer dataframe
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
        return x_train, y_train
    end
    testdata = vardata[testidx, :]
    x_test = Matrix{dtype}(testdata[!, Not(target)])'
    y_test = Vector{dtype}(testdata[!, target])' ./ dtype(target_scale)
    return x_train, y_train, x_test, y_test
end

"""
    snotel_metadata(; fields)

Return a database of snotel station metadata for usage in dataset creation.

# Arguments
- `fields::Vector{String}`: optional list of specific metadata fields to extract. Default is
[stationID, state.code, "elevation", "latitude", "longitude"].
"""
function snotel_metadata(;
    fields::Vector{String} = [
        "stationId",
        "state.code",
        "elevation",
        "latitude,longitude",
    ],
)
    link = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultipleStationReport,metric/daily/start_of_period/network=%22SNTL%22%7Cname/0,0/"
    for field in fields
        link = link * field * ","
    end
    data = CSV.read(HTTP.get(link).body, DataFrame, comment = "#", delim = ",")
    return data
end
end
