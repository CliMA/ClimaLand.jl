""" This file provides methods for interacting with the MODIS REST API at a site
    in order to downloads specific subsets of MODIS products at a site level."""

using HTTP
using JSON
using Statistics
using Dates

function send_get_subset(
    product::String,
    start_date::String,
    end_date::String,
    site::String;
    network::String = "AMERIFLUX",
    band::String = "",
)
    """Sends a GET request to the MODIS REST API for a specific product, time, 
    and site. Returns an HTTP.messages.response object."""
    url = "https://modis.ornl.gov/rst/api/v1/$product/$network/$site/subset?"
    band = band == "" ? band : "band=$band"
    start_date = "&startDate=$start_date"
    end_date = "&endDate=$end_date"
    url *= "$band$start_date$end_date"
    return HTTP.get(url)
end

function check_response(response::HTTP.Messages.Response)
    """Checks the response from the MODIS REST API for errors. Returns the 
    response if no errors are found, otherwise throws an error."""
    if response.status != 200
        throw(ErrorException("Error in MODIS pull: $(response.status)"))
    end
    return response
end

function parse_response(response::HTTP.Messages.Response)
    """Parses the response from the MODIS REST API into a nested Julia 
    dictionary using the JSON module and returns the dictionary."""
    return JSON.parse(String(response.body))["subset"]
end

function single_col_data_matrix(JSON_data::Vector)
    """Takes in a vector of JSON data dicitonaries of a single column of chunked
    MODIS data and assembles it into a data matrix with a time column and a
    single data column. The data column gives the average of the grid cell 
    arround the site for each day, and the time column gives the Dates.DateTime 
    object for each day. In averaging, MODIS data marked as missing (>100) is 
    ignored."""
    cleansed_dat = [
        [
            chunk["data"][i] < 100 ? chunk["data"][i] : missing for
            i in 1:length(chunk["data"])
        ] for chunk in JSON_data
    ]
    data_col =
        [mean(filter(!ismissing, day_dat)) / 10 for day_dat in cleansed_dat]
    time_col =
        [DateTime(JSON_data[i]["calendar_date"]) for i in 1:length(JSON_data)]
    return hcat(time_col, data_col)
end
