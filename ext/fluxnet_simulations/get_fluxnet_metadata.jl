using DelimitedFiles
using ClimaLand

function parse_array_field(field)
    if field == "" || field == "NaN"
        return Float64[]
    elseif field isa Float64
        return isnan(field) ? Float64[] : [field]
    else
        return parse.(Float64, split(String(field), ";"))
    end
end


"""
    get_site_info(site_ID)

This function retrieves the site information from an artifact for a given fluxnet site ID.
    
Returns a Float64-valued namedtuple containing: 
- `lat`: Latitude of the site (deg)
- `long`: Longitude of the site (deg)
- `time_offset`: Time offset from UTC (hours). Note that the ClimaLand convention is for
    UTC = local_time + time_offset.
- `atmospheric_sensor_height` (Vector): Heights of atmospheric sensors (m)
"""
function get_site_info(site_ID; fluxnet2015_metadata_path = nothing)
    if fluxnet2015_metadata_path === nothing
        # Use the default path from ClimaLand.Artifacts
        fluxnet2015_data_path = ClimaLand.Artifacts.fluxnet2015_data_path()
        fluxnet2015_metadata_path =
            joinpath(fluxnet2015_data_path, "metadata_DD_clean.csv")
    end
    raw = readdlm(fluxnet2015_metadata_path, ',', Any)
    header = raw[1, :]
    data = raw[2:end, :]

    # Find the row matching site_ID
    row_idx = findfirst(row -> row[1] == site_ID, eachrow(data))
    if isnothing(row_idx)
        error("Site ID $site_ID not found in metadata.")
    end

    site_metadata = data[row_idx, :]

    varnames =
        ["latitude", "longitude", "utc_offset", "atmospheric_sensor_heights"]
    column_name_map = Dict(
        varname => findfirst(header[:] .== varname) for varname in varnames
    )

    nan_if_empty(x) =
        (x isa String && isempty(x)) ? NaN :
        (x isa AbstractFloat && isnan(x)) ? NaN : x

    for (varname, column) in column_name_map
        field = site_metadata[column]
        if field == "" || (field isa AbstractFloat && isnan(field))
            @warn "Field $(varname) is missing for site $site_ID"
        end
    end

    return (;
        lat = nan_if_empty(site_metadata[column_name_map["latitude"]]),
        long = nan_if_empty(site_metadata[column_name_map["longitude"]]),
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
        time_offset = Int64(
            -1 * nan_if_empty(site_metadata[column_name_map["utc_offset"]]),
=======
        time_offset = -1 * nan_if_empty(
            site_metadata[column_name_map["utc_offset"]],
>>>>>>> 4ccce0bef (for FluxnetSimulationsExt module, wrote access functions to get simulation info + parameters for 4 specific sites with hardcoded information and any other site with mapped info)
        ),
=======
        time_offset = Int64(-1 * nan_if_empty(
            site_metadata[column_name_map["utc_offset"]],
        )),
>>>>>>> a523e6aba (pain)
=======
        time_offset = Int64(
            -1 * nan_if_empty(site_metadata[column_name_map["utc_offset"]]),
        ),
>>>>>>> a199d603d (cleaning up for a branch clone)
        atmospheric_sensor_height = parse_array_field(
            site_metadata[column_name_map["atmospheric_sensor_heights"]],
        ),
    )
end


"""
    get_canopy_height(site_ID)

This function retrieves a canopy height from an artifact for a given fluxnet site ID.
    
Returns a Float64-valued namedtuple containing: 
- `canopy_height`: Canopy height of the site (m)
"""
function get_canopy_height(site_ID; fluxnet2015_metadata_path = nothing)
    if fluxnet2015_metadata_path === nothing
        # Use the default path from ClimaLand.Artifacts
        fluxnet2015_data_path = ClimaLand.Artifacts.fluxnet2015_data_path()
        fluxnet2015_metadata_path =
            joinpath(fluxnet2015_data_path, "metadata_DD_clean.csv")
    end

    raw = readdlm(fluxnet2015_metadata_path, ',', Any)
    header = raw[1, :]
    data = raw[2:end, :]

    site_info = get_site_info(site_ID; fluxnet2015_metadata_path)

    # Find the row matching site_ID
    row_idx = findfirst(row -> row[1] == site_ID, eachrow(data))
    if isnothing(row_idx)
        error("Site ID $site_ID not found in metadata.")
    end

    site_metadata = data[row_idx, :]

    varname = "canopy_height"

    col_idx = findfirst(header[:] .== varname)

    nan_if_empty(x) =
        (x isa String && isempty(x)) ? NaN :
        (x isa AbstractFloat && isnan(x)) ? NaN : x

    canopy_height = site_metadata[col_idx]
    if canopy_height == "" ||
       (canopy_height isa AbstractFloat && isnan(canopy_height))
        @warn "Field $(varname) is missing for site $site_ID"
    end

    return nan_if_empty(canopy_height)
end
