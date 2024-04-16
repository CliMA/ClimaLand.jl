export vairaranch_default_configs

"""
function vairaranch_default_configs(;
    data_link = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
    LAI_link = "https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv",
    atmos_h = FT(2),
    local_to_UTC = 8,
    lat = FT(38.4133), # degree
    long = FT(-120.9508), # degree
	)
	return (
	    data_link = data_link,
		LAI_link = LAI_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,	
	)
end

Named Tuple containing default configurations for the Vaira Ranch site. 
"""
function vairaranch_default_configs(;
    data_link = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
    LAI_link = "https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv",
    atmos_h = FT(2),
    local_to_UTC = 8,
    lat = FT(38.4133), # degree
    long = FT(-120.9508), # degree
)
    return (
        data_link = data_link,
        LAI_link = LAI_link,
        local_to_UTC = local_to_UTC,
        atmos_h = atmos_h,
        lat = lat,
        long = long,
    )
end
