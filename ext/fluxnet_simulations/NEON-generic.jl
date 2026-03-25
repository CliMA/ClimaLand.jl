const NEON_SITE_METADATA_CSV =
    "/kiwi-data/Data/groupMembers/evametz/ERA5/sitedata/NEON_Field_Site_Metadata_20260324.csv"

function _neon_site_key(site_ID)
    return uppercase(replace(string(site_ID), "NEON_" => "", "NEON-" => ""))
end

function _get_neon_site_metadata(site_ID)
    (data, columns) = readdlm(NEON_SITE_METADATA_CSV, ','; header = true)
    header = vec(String.(columns))

    i_site = findfirst(==("site_id"), header)
    i_lat = findfirst(==("latitude"), header)
    i_long = findfirst(==("longitude"), header)
    i_h = findfirst(==("tower_height_m"), header)

    key = _neon_site_key(site_ID)
    row_idx = findfirst(i -> uppercase(string(data[i, i_site])) == key, axes(data, 1))

    lat = parse(Float64, string(data[row_idx, i_lat]))
    long = parse(Float64, string(data[row_idx, i_long]))
    atmos_h = parse(Float64, string(data[row_idx, i_h]))
    return (; lat, long, atmos_h)
end


"""
    get_domain_info(FT, ::Val{:US_Var}; dz_tuple = nothing,
        nelements = 14, zmin = FT(-0.5), zmax = FT(0))

Gets and returns primary domain information for the NEON CPER, which is a grassland,
NEON site. The values are provided as defaults, and can be overwritten by passing the corresponding
keyword arguments to this function.


"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{site_ID};
    #dz_bottom = FT(1.5),
    #dz_top = FT(0.1),
    #nelements = 20,
    #zmin = FT(-10),
    #zmax = FT(0),
    dz_bottom = FT(2), #FT(1.5),
    dz_top = FT(0.038),
    nelements = 24,
    zmin = FT(-6.2),
    zmax = FT(0),
)
    dz_tuple = (dz_bottom, dz_top)

    return (; dz_tuple, nelements, zmin, zmax)
end


function FluxnetSimulations.get_location(
    FT,
    ::Val{site_ID};
    time_offset = 0,
    lat = nothing,
    long = nothing,
) where {site_ID}
    metadata = _get_neon_site_metadata(site_ID)
    lat = isnothing(lat) ? FT(metadata.lat) : lat
    long = isnothing(long) ? FT(metadata.long) : long
    return (; time_offset, lat, long)
end

function FluxnetSimulations.get_fluxtower_height(
    FT,
    ::Val{site_ID};
    atmos_h = nothing,
) where {site_ID}
    metadata = _get_neon_site_metadata(site_ID)
    atmos_h = isnothing(atmos_h) ? FT(metadata.atmos_h) : atmos_h
    return (; atmos_h,)
end
