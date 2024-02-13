export vairaranch_default_configs

"""
    vairaranch_default_configs = (
        data_link = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
        LAI_link = "https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv",
        atmos_h = FT(2),
        local_to_UTC = 8,
        lat = FT(38.4133), # degree
        long = FT(-120.9508), # degree
    )

Named Tuple containing default configurations for the Vaira Ranch site. 
"""
vairaranch_default_configs = (
    data_link = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
    LAI_link = "https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv",
    atmos_h = FT(2),
    local_to_UTC = 8,
    lat = FT(38.4133), # degree
    long = FT(-120.9508), # degree
)
