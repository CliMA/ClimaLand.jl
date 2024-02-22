export niwotridge_default_args

"""
    niwotridge_default_args = (
        dz_bottom = FT(1.25),
        dz_top = FT(0.05),
        n_stem = Int64(1),
        n_leaf = Int64(1),
        h_leaf = FT(6.5), # m,
        h_stem = FT(7.5), # m,
        t0 = Float64(120 * 3600 * 24), # start mid year to avoid snow,
        dt = Float64(40),
        n = 45,
    )

Named Tuple containing default simulation parameter for the Niwot Ridge site. 
"""
niwotridge_default_args = (
    dz_bottom = FT(1.25),
    dz_top = FT(0.05),
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_leaf = FT(6.5), # m,
    h_stem = FT(7.5), # m,
    t0 = Float64(120 * 3600 * 24), # start mid year to avoid snow,
    dt = Float64(40),
    n = 45,
)
