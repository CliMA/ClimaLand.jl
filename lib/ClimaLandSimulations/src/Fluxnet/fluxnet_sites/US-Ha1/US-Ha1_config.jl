export harvard_default_configs

"""
	function harvard_default_configs(;
		data_link = "https://caltech.box.com/shared/static/xixaod6511cutz51ag81k1mtvy05hbol.csv",
		local_to_UTC = 5,
		lat = FT(42.5378), # degree
		long = FT(-72.1715), # degree
		atmos_h = FT(30),
		)
		return (
			data_link = data_link,
			local_to_UTC = local_to_UTC,
			atmos_h = atmos_h,
			lat = lat,
			long = long,	
		)
	end

Named Tuple containing default configurations for the Harvard site. 
"""
function harvard_default_configs(;
    data_link = "https://caltech.box.com/shared/static/xixaod6511cutz51ag81k1mtvy05hbol.csv",
    local_to_UTC = 5,
    lat = FT(42.5378), # degree
    long = FT(-72.1715), # degree
    atmos_h = FT(30),
)
    return (
        data_link = data_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,
    )
end
