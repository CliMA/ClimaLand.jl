export niwotridge_default_configs

"""
	function niwotridge_default_configs(;
		data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv",
		local_to_UTC = 7,
		lat = FT(40.0329), # degree,
		long = FT(-105.5464), # degree,
		atmos_h = FT(21.5),
		)
		return (
			data_link = data_link,
			local_to_UTC = local_to_UTC,
			atmos_h = atmos_h,
			lat = lat,
			long = long,	
		)
	end

Named Tuple containing default configurations for the Niwot Ridge site. 
"""
function niwotridge_default_configs(;
    data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv",
    local_to_UTC = 7,
    lat = FT(40.0329), # degree,
    long = FT(-105.5464), # degree,
    atmos_h = FT(21.5),
)
    return (
        data_link = data_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,
    )
end
