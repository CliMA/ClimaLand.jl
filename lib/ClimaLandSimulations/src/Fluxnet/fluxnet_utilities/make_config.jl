export make_config

"""
    make_config()

metadata for site:
- data_link, data file download link
- local_to_UTC, timezone (offset from UTC in hrs)
- atmos_h, height of sensor on flux tower
- lat, site latitude
- long, longitude
"""
function make_config(; data_link, local_to_UTC, atmos_h, lat, long)
    return (
        data_link = data_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,
    )
end

function make_config(site_ID)
    default_configs = Dict(
        "US-MOz" => ozark_default_configs,
        "US-Ha1" => harvard_default_configs,
        "US-NR1" => niwotridge_default_configs,
        "US-Var" => vairaranch_default_configs,
    )

    if haskey(default_configs, site_ID)
        make_config(; default_configs[site_ID]...)
    end
end
