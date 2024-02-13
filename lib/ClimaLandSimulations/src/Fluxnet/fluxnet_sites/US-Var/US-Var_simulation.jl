export vairaranch_default_args

"""
	function vairaranch_default_args(;
		dz_bottom = FT(1.0),
		dz_top = FT(0.05),
		n_stem = Int64(0),
		n_leaf = Int64(1),
		h_leaf = FT(0.7), # m,
		h_stem = FT(0), # m,
		t0 = Float64(21 * 3600 * 24), # start day 21 of the year,
		dt = Float64(40),
		n = 45,
		)
		return (
			dz_bottom = dz_bottom,
			dz_top = dz_top,
			n_stem = n_stem,
			n_leaf = n_leaf,
			h_stem = h_stem,
			h_leaf = h_leaf,
			t0 = t0,
			dt = dt,
			n = n,	
		)
	end

Named Tuple containing default simulation parameter for the Vaira Ranch site. 
"""
function vairaranch_default_args(;
    dz_bottom = FT(1.0),
    dz_top = FT(0.05),
    n_stem = Int64(0),
    n_leaf = Int64(1),
    h_leaf = FT(0.7), # m,
    h_stem = FT(0), # m,
    t0 = Float64(21 * 3600 * 24), # start day 21 of the year,
    dt = Float64(40),
    n = 45,
)
    return (
        dz_bottom = dz_bottom,
        dz_top = dz_top,
        n_stem = n_stem,
        n_leaf = n_leaf,
        h_stem = h_stem,
        h_leaf = h_leaf,
        t0 = t0,
        dt = dt,
        n = n,
    )
end
