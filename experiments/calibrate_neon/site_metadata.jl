"""
NEON site metadata helper functions.
Loads site information (latitude, longitude, tower height) from CSV.
"""

using DelimitedFiles

const NEON_SITE_METADATA_CSV =
    "/kiwi-data/Data/groupMembers/evametz/ERA5/sitedata/NEON_Field_Site_Metadata_20260324.csv"

"""
    _neon_site_key(SITE_ID)

Convert NEON site ID to lookup key (uppercase, remove NEON prefix).
"""
function _neon_site_key(SITE_ID)
    return uppercase(replace(string(SITE_ID), "NEON_" => "", "NEON-" => ""))
end

"""
    _get_neon_site_metadata(SITE_ID)

Load site metadata (latitude, longitude, tower_height_m) from CSV.
Returns named tuple: (; lat, long, atmos_h)
"""
function _get_neon_site_metadata(SITE_ID)
    (data, columns) = DelimitedFiles.readdlm(NEON_SITE_METADATA_CSV, ','; header = true)
    header = vec(String.(columns))

    i_site = findfirst(==("site_id"), header)
    i_lat = findfirst(==("latitude"), header)
    i_long = findfirst(==("longitude"), header)
    i_h = findfirst(==("tower_height_m"), header)

    key = _neon_site_key(SITE_ID)
    row_idx = findfirst(i -> uppercase(string(data[i, i_site])) == key, axes(data, 1))

    lat = parse(Float64, string(data[row_idx, i_lat]))
    long = parse(Float64, string(data[row_idx, i_long]))
    atmos_h = parse(Float64, string(data[row_idx, i_h]))
    return (; lat, long, atmos_h)
end
