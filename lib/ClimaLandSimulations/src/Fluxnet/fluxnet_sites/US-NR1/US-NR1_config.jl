export niwotridge_default_configs

"""
    niwot_default_configs = (
        data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv",
        local_to_UTC = 7,
        lat = FT(40.0329), # degree,
        long = FT(-105.5464), # degree,
        atmos_h = FT(21.5),
    )

Named Tuple containing default configurations for the Niwot Ridge site. 
"""
niwotridge_default_configs = (
    data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv",
    local_to_UTC = 7,
    lat = FT(40.0329), # degree,
    long = FT(-105.5464), # degree,
    atmos_h = FT(21.5),
)
