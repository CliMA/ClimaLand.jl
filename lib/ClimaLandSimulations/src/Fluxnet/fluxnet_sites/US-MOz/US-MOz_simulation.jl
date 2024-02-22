export ozark_default_args

"""
    ozark_default_args = (
    dz_bottom = FT(1.5),
    dz_top = FT(0.025),
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(9),
    h_leaf = FT(9.5),
    t0 = Float64(120 * 3600 * 24), 
    dt = Float64(120),
    n = 15
   )

Named Tuple containing default simulation parameter for the Ozark site. 
"""
ozark_default_args = (
    dz_bottom = FT(1.5),
    dz_top = FT(0.025),
    n_stem = Int64(1),
    n_leaf = Int64(1),
    h_stem = FT(9),
    h_leaf = FT(9.5),
    t0 = Float64(120 * 3600 * 24),
    dt = Float64(120),
    n = 15,
)
