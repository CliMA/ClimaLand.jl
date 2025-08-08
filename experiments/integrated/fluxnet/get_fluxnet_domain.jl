using DelimitedFiles
using ClimaLand



function parse_array_field(field)
    @show field
    if field == "" || field == "NaN"
        return Float64[]
    elseif field isa Float64
        return isnan(field) ? Float64[] : [field]
    else
        return parse.(Float64, split(String(field), ";"))
    end
end


"""
    get_site_domain(site_ID)

This function retrieves the domain information for a given fluxnet site ID.
    
Returns a Float64-valued namedtuple containing: 
- `lat`: Latitude of the site (deg)
- `long`: Longitude of the site (deg)
- `time_offset`: Time offset from UTC (hours). Note that the ClimaLand convention is for
    UTC = local_time + time_offset. 
- `canopy_height`: Average height of the canopy (m)
- `atmospheric_sensor_height` (Vector): Heights of atmospheric sensors (m)
- `swc_depths` (Vector): Soil water content measurement depths (m)
- `ts_depths` (Vector): Soil temperature measurement depths (m)
"""
function get_site_domain(site_ID; fluxnet2015_metadata_path = nothing)
    if fluxnet2015_metadata_path === nothing
        # Use the default path from ClimaLand.Artifacts
        fluxnet2015_data_path = ClimaLand.Artifacts.fluxnet2015_data_path()
        fluxnet2015_metadata_path = joinpath(fluxnet2015_data_path, "metadata_DD_clean.csv")
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

    varnames = [
        "latitude", "longitude", "utc_offset", "canopy_height", "atmospheric_sensor_heights",
        "swc_depths", "ts_depths"
    ]
    column_name_map = Dict(
        varname => findfirst(header[:] .== varname) for varname in varnames
    )

    nan_if_empty(x) =
        (x isa String && isempty(x)) ? NaN :
        (x isa AbstractFloat && isnan(x)) ? NaN :
        x

    for (varname, column) in column_name_map
        field = site_metadata[column]
        if field == "" || (field isa AbstractFloat && isnan(field))
            @warn "Field $(varname) is missing for site $site_ID"
        end
    end
    
    return (;
        lat = nan_if_empty(site_metadata[column_name_map["latitude"]]),
        long = nan_if_empty(site_metadata[column_name_map["longitude"]]),
        time_offset = -1 * nan_if_empty(site_metadata[column_name_map["utc_offset"]]),
        h_canopy = nan_if_empty(site_metadata[column_name_map["canopy_height"]]),
        atmospheric_sensor_height = parse_array_field(site_metadata[column_name_map["atmospheric_sensor_heights"]]),
        swc_depths = parse_array_field(site_metadata[column_name_map["swc_depths"]]),
        ts_depths = parse_array_field(site_metadata[column_name_map["ts_depths"]]),
    )
end
