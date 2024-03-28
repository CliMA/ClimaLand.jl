# # Scraping SNOTEL Data

# This tutorial shows you how to make use of the code developed
# for scraping SNOTEL site data in order to generate datasets for use
# in training artificial intelligence models for seasonal snow forecasting.
# The code below contains a commented version of the code used to produce
# `cleandata.csv`, which is used in the base tutorial for snow forecasting,
# as well as the paper. However, exploration of the optional arguments
# or requesting of alternative [SNOTEL data codes](https://www.nrcs.usda.gov/wps/portal/wcc/home/dataAccessHelp/webService/webServiceReference#elementCodes) offers
# additional utility in creating alternative data sets for further investigation.

# We begin by importing the data tools module, as well as other required packages:
using ClimaLand
using DataFrames, CSV, HTTP, Dates, Flux, StatsBase, cuDNN

# The code lives in an extenson that we have to manually load. The extension can
# be loaded only if "CSV", "HTTP", "Flux", "StatsBase", "cuDNN" and "ClimaLand"
# are loaded.
DataTools = Base.get_extension(ClimaLand, :NeuralSnowExt).DataTools
ModelTools = Base.get_extension(ClimaLand, :NeuralSnowExt).ModelTools


# We then define constants that will be used in the cleaning of the SNOTEL data,
# such as conversion constants from imperial to metric units, and the sensor limits
# defined in the [SNOTEL Engineering Handbook](https://directives.sc.egov.usda.gov/OpenNonWebContent.aspx?content=27630.wba). Some SNOTEL sensors measure
# in imperial units, and some measure in metric units, and the data portal will round
# converted values if a sensor stream is requested in units other than its original
# measurement. Therefore, we will scrape data in the originally measured units to limit
# systemic errors.

const inch2meter = 0.0254
const kmphr2mps = 5.0 / 18.0

filter_val = Dict{Symbol, Tuple{Real, Real}}(
    :SWE => (0.0, 250.0),
    :z => (0.0, 420.0),
    :precip => (0.0, 250.0),
    :rel_hum_avg => (10.0, 100.0),
    :sol_rad_avg => (0.0, 1500.0),
    :wind_speed_avg => (0.0, 216.0),
    :air_temp_avg => (-55.0, 60.0),
)

scales = Dict{Symbol, Real}(
    :SWE => inch2meter,
    :z => inch2meter,
    :precip => inch2meter,
    :rel_hum_avg => 0.01,
    :wind_speed_avg => kmphr2mps,
);

# We next proceed to outline which stations will be scraped by defining
# a dictionary of station IDs, paired with the date range to be scraped.
# "`start`" refers to 1850-01-01 or the first available date, while "`end`"
# refers to 2023-01-01 or the last available date. Most of these stations
# are commented out for the sake of speed and readability in generating
# the tutorial, but can be uncommented to yield the full dataset found in
# `cleandata.csv` used in the [base tutorial](../base_tutorial/). Stations were
# selected based upon their availability of the features utilized in
# creating the model used in the paper.

good_stations = Dict{Int, Tuple{String, String}}(
    1030 => ("start", "2016-01-01"),
    1053 => ("2010-01-01", "end"),
    #=1083 => ("2013-01-01", "end"),
    1105 => ("start", "2012-06-01"),
    1122 => ("start", "end"),
    1123 => ("start", "end"),
    1159 => ("start", "end"),
    1168 => ("start", "end"),
    1170 => ("start", "end"),
    1254 => ("2018-01-01", "end"),
    1286 => ("start", "end"),
    306 => ("2013-01-01", "2020-01-01"),
    316 => ("start", "end"),
    344 => ("start", "end"),
    367 => ("2022-09-01", "2023-04-01"),
    395 => ("start", "end"),
    457 => ("2008-01-01", "end"),
    482 => ("start", "2010-01-01"),
    491 => ("2013-01-01", "end"),
    532 => ("2013-01-01", "2022-01-01"),
    551 => ("start", "end"),
    571 => ("start", "end"),
    599 => ("start", "end"),
    608 => ("start", "2018-01-01"),
    613 => ("2007-01-01", "2015-01-01"),
    665 => ("start", "2022-01-01"),
    734 => ("start", "end"),
    737 => ("start", "end"),
    744 => ("2014-01-01", "2016-01-01"),
    832 => ("start", "end"),
    845 => ("2014-01-01", "end"),
    854 => ("2019-01-01", "end"),
    857 => ("start", "end"),
    921 => ("2019-01-01", "end"),
    922 => ("2015-01-01", "2018-01-01"),
    942 => ("2009-01-01", "2018-01-01"),
    969 => ("2011-01-01", "end"),
    974 => ("start", "end"),
    978 => ("2005-01-01", "2018-01-01"),=#
);

# We also extract a `DataFrame` matching station ID to various station metadata,
# in order to automate some of the scraping process and pass some station
# metadata that is used for analysis in the paper. This resulting `DataFrame`
# can also be used to see other available SNOTEL station IDs for scraping,
# in order to create custom datasets.
metadata = DataTools.snotel_metadata();
metacols = ["id", "state", "elev", "lat", "lon"]
DataFrames.rename!(metadata, Symbol.(metacols));

# We then loop through each site to scrape and follow an automated data pipeline, consisting of:
# - Extracting the daily and hourly timeseries from the site
# - Applying the sensor bounds over each data timeseries (i.e. remove sensor error)
# - Converting the hourly dataset into a daily dataset
# - Filling holes/features in the daily series with the hourly series
# - Scaling all data to the appropriate metric units
# - Restricting data to complete cases
# - Making the differential variables ( ``\frac{dz}{dt}``, etc.)
# - Removing negative precipitation cases (i.e. where the water year resets, or sensor error)
# - Attaching appropriate metadata

allsites = Any[];
for site in sort(collect(keys(good_stations)))
    state = metadata[metadata[!, :id] .== site, :state][1]
    start_date = good_stations[site][1]
    end_date = good_stations[site][2]
    daily = DataTools.apply_bounds(
        DataTools.sitedata_daily(
            site,
            state,
            start = start_date,
            finish = end_date,
        ),
        filter_val,
    )
    hourly = DataTools.apply_bounds(
        DataTools.sitedata_hourly(
            site,
            state,
            start = start_date,
            finish = end_date,
        ),
        filter_val,
    )
    hourly_d = DataTools.hourly2daily(hourly)
    gap_daily = DataTools.rectify_daily_hourly(daily, hourly_d)
    daily_scaled = DataTools.scale_cols(gap_daily, scales)
    daily_clean = daily_scaled[DataTools.completecases(daily_scaled), :]
    daily_clean = DataTools.makediffs(daily_clean, Day(1))
    good_vals = daily_clean[!, :dprecipdt] .>= 0.0
    daily_clean = daily_clean[good_vals, Not(:precip)]
    #show(describe(daily_clean), allrows = true, allcols = true)
    #print("\nSIZE: ", nrow(daily_clean), "\n")

    daily_clean[!, :id] .= site
    daily_clean[!, :elev] .= metadata[metadata[!, :id] .== site, :elev][1]
    daily_clean[!, :lat] .= metadata[metadata[!, :id] .== site, :lat][1]
    daily_clean[!, :lon] .= metadata[metadata[!, :id] .== site, :lon][1]

    push!(allsites, daily_clean)
end;

# With the sites complete, we condense all sites into a single `DataFrame`,
totaldata = deepcopy(allsites[1])
for site in allsites[2:end]
    append!(totaldata, site)
end
# and a final `CSV.write("newcleanfiletestdata.csv", totaldata)` call will
# save the file.
