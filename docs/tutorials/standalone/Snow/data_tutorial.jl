# # Scraping SNOTEL Data

# This tutorial shows you how to make use of the code developed
# for scraping SNOTEL site data in order to generate datasets for use
# in training artificial intelligence models for seasonal snow forecasting.
# The code below contains a basic version of the code used to produce
# `training_data.csv`, which is used in the [base tutorial](../base_tutorial/) for snow forecasting,
# as well as the paper. However, exploration of the optional arguments
# or requesting of alternative [SNOTEL data codes](https://www.nrcs.usda.gov/wps/portal/wcc/home/dataAccessHelp/webService/webServiceReference#elementCodes) offers
# additional utility in creating alternative data sets for further investigation.

# We begin by importing all required packages:
using ClimaLand
using DataFrames, CSV, HTTP, Dates, Flux, StatsBase, cuDNN

# The code lives in an extenson that we have to manually load. The extension can
# be loaded only if "CSV", "HTTP", "Flux", "StatsBase", "cuDNN" and "ClimaLand"
# are loaded.
DataTools = Base.get_extension(ClimaLand, :NeuralSnowExt).DataTools;

# We first extract a `DataFrame` matching station ID to various station metadata,
# in order to automate some of the scraping process and pass some station
# metadata that is used for analysis in the paper. This resulting `DataFrame`
# can also be used to see other available SNOTEL station IDs for scraping,
# in order to create custom datasets.

# ```julia
#metadata = DataTools.snotel_metadata();
#metacols = ["id", "name", "state", "elev", "lat", "lon"]
#DataFrames.rename!(metadata, Symbol.(metacols));
# ```

# At the most user-friendly level, the function `scrape_site_paper()` provides
# a wrapper to scrape SNOTEL data in the exact same manner as the paper (it may take
# a minute or two per site). This function handles all special cases and data processing,
# allowing the user to only pass a SNOTEL ID number and associated state code to retrieve
# the same data as that used in the paper. However, this will likely not work or yield
# unexpected results for sites not used in the paper. Here is an example for
# how to use the metadata to streamline the process:

# ```julia
# example_ID = 1030
# example_state = metadata[findfirst(==(example_ID), metadata[!, :id]), :state]
# example_data = DataTools.scrape_site_paper(example_ID, example_state);
# ```

# And that's it! This can be iterated within a loop to gather the data for all
# sites. However, while straightforward, this wrapper obfuscates many of the
# underlying steps, or some of the opportunities for using different arguments
# to generate custom datasets. As such, we can reimplement much of the same code
# in more detail below to enable more advanced usage.

# We first define constants that will be used in the cleaning of the SNOTEL data,
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
    :air_temp_avg => (-40.0, 60.0),
)

scales = Dict{Symbol, Real}(
    :SWE => inch2meter,
    :z => inch2meter,
    :precip => inch2meter,
    :rel_hum_avg => 0.01,
    :wind_speed_avg => kmphr2mps,
);

# We next proceed to outline which stations will be scraped by defining
# a dictionary of station IDs, paired with the date range to be scraped if
# a custom range is desired. `"start"` refers to 1850-01-01 or the first available
# date, while `"end"` refers to the earlier option bewteen 2024-02-01 or the last available date. Most of these stations
# are commented out for the sake of speed and readability in generating
# the tutorial, or due to special handling required, but can be uncommented
# to yield the full dataset (if special cases are handled) found in
# `training_data.csv` used in the [base tutorial](../base_tutorial/). Stations were
# selected based upon their availability of the features utilized in
# creating the model used in the paper:

# - `*` Indicates alternative handling of the `rectify_daily_hourly()` function.

# - `^` Indicates usage of `RHUM` flag instead of `RHUMV` flag for relative humidity.

# - `A` Indicates an Alaskan site, which is in the testing data, not the training data, and uses a lower temperature bound of -50 instead of -40 in `filter_val`.

# - `T` Requires a site that already has had the temperature bias correction at the portal level as of May 2024.

# - `X` Indicates a SNOTEL portal error when trying to scrape into 2024, as of May 2024.

good_stations = Dict{Int, Tuple{String, String}}(
    #306 => ("start", "end"), #*
    316 => ("start", "end"),
    344 => ("start", "end"),
    #=367 => ("start", "end"),
    395 => ("start", "end"),
    457 => ("start", "end"),
    482 => ("start", "end"),
    491 => ("start", "end"),
    515 => ("start", "2023-06-02"), #X
    532 => ("start", "end"),
    551 => ("start", "end"),
    571 => ("start", "end"),
    599 => ("start", "end"),
    608 => ("start", "end"),
    613 => ("start", "end"),
    641 => ("start", "end"), #A^
    665 => ("start", "end"),
    708 => ("start", "end"),
    715 => ("start", "end"),
    734 => ("start", "end"),
    737 => ("start", "end"),
    744 => ("start", "end"),
    832 => ("start", "end"),
    845 => ("start", "end"),
    854 => ("start", "end"),
    857 => ("start", "end"),
    921 => ("start", "end"),
    922 => ("start", "end"),
    927 => ("start", "end"),
    942 => ("start", "end"),
    963 => ("start", "end"), #A^
    969 => ("start", "end"),
    974 => ("start", "end"),
    978 => ("start", "end"), #*
    1030 => ("start", "end"),
    1035 => ("start", "end"), #A^
    1053 => ("start", "end"),
    1070 => ("start", "end"), #A^T
    1083 => ("start", "end"),
    1091 => ("start", "end"), #A^T
    1092 => ("start", "end"), #A^T
    1105 => ("start", "end"),
    1122 => ("start", "end"), #*
    1123 => ("start", "end"),
    1159 => ("start", "end"),
    1168 => ("start", "end"),
    1170 => ("start", "end"),
    1254 => ("start", "end"),
    1286 => ("start", "end"),
    2080 => ("start", "end"), #A^
    2170 => ("start", "end"), #^
    =#
);

# We then loop through each site to scrape and follow an automated data pipeline, consisting of:
# - Extracting the daily and hourly timeseries from the site
# - Applying the sensor bounds over each data timeseries (i.e. remove sensor error)
# - Converting the hourly dataset into a daily dataset
# - Coalescing the converted-hourly and daily data into one dataset
# - Scaling all data to the appropriate metric units
# - Restricting data to complete cases
# - Making the differential variables ( ``\frac{dz}{dt}``, etc.)
# - Resetting negative precipitation cases (i.e. where the water year resets), and using daily precipitation rates `dprecipdt` instead of accumulated precipitation `precip`
# - Attaching appropriate metadata

# A few steps are commented out, which indicate steps implemented in `scrape_site_paper()` 
# like quality-control measures, which could be substituted with other user-defined steps.

# ```julia
# allsites = Any[];
# for site in sort(collect(keys(good_stations)))
#     state = metadata[metadata[!, :id] .== site, :state][1]
#     start_date = good_stations[site][1]
#     end_date = good_stations[site][2]
# 
#     hourly = DataTools.apply_bounds(
#         DataTools.sitedata_hourly(
#              site,
#             state,
#             start = start_date,
#             finish = end_date,
#         ),
#        filter_val,
#     )
#     hourly[!, :id] .= site
#     #hourly = DataTools.bcqc_hourly(hourly)
#     hourly_d = DataTools.hourly2daily(hourly)
#     #DataFrames.allowmissing!(hourly_d) 
#     #sflags = DataTools.qc_filter(hourly_d, :sol_rad_avg, t1 = 2)
#    #hourly_d[sflags, :sol_rad_avg] .= missing
# 
#     daily = DataTools.apply_bounds(
#         DataTools.sitedata_daily(
#             site,
#             state,
#             start = start_date,
#             finish = end_date,
#         ),
#         filter_val,
#     )
#     daily[!, :id] .= site
#     gap_daily = DataTools.rectify_daily_hourly(daily, hourly_d)
#     #gap_daily = DataTools.bcqc_daily(gap_daily, site, state)
#     #gap_daily = DataTools.d_impute(gap_daily)
#     daily_scaled = DataTools.scale_cols(gap_daily, scales)
#     daily_clean = daily_scaled[completecases(daily_scaled), :]
#     daily_clean = DataTools.makediffs(daily_clean, Day(1))
#     good_vals = daily_clean[!, :dprecipdt] .>= 0.0
#     daily_clean[(!).(good_vals), :dprecipdt] .= 0.0
#     daily_clean = daily_clean[!, Not(:precip)]
#     #show(describe(daily_clean), allrows = true, allcols = true)
#     #print("\nSIZE: ", nrow(daily_clean), "\n")
# 
#     daily_clean[!, :id] .= site
#     daily_clean[!, :elev] .= metadata[metadata[!, :id] .== site, :elev][1]
#     daily_clean[!, :lat] .= metadata[metadata[!, :id] .== site, :lat][1]
#     daily_clean[!, :lon] .= metadata[metadata[!, :id] .== site, :lon][1]
# 
#     push!(allsites, daily_clean)
# end;
# ```

# With the sites complete, we condense all sites into a single `DataFrame`,
# ```julia
# totaldata = deepcopy(allsites[1])
# for site in allsites[2:end]
#     append!(totaldata, site)
# end
# ```

# and a final `CSV.write("data.csv", totaldata)` call will
# save the file.

# Many of the functions above contain default or optional arguments which can
# be explored to obtain a richer set of functionality, or implement some of
# the special cases mentioned above. Such options can be explored in the [code documentation](https://github.com/CliMA/ClimaLand.jl/blob/main/ext/neural_snow/DataTools.jl).
