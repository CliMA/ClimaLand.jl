export ozark_default_configs

"""
	function ozark_default_configs(;
		data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv",
		local_to_UTC = 7,
		atmos_h = FT(32),
		lat = FT(38.7441),
		long = FT(-92.2000),
		)
		return (
			data_link = data_link,
			local_to_UTC = local_to_UTC,
			atmos_h = atmos_h,
			lat = lat,
			long = long,
		)
	end

Named Tuple containing default configurations for the Ozark site. 
"""
function ozark_default_configs(;
    data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv",
    local_to_UTC = 7,
    atmos_h = FT(32),
    lat = FT(38.7441),
    long = FT(-92.2000),
)
    return (
        data_link = data_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,
    )
end
