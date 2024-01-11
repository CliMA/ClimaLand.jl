# Domain setup
# For soil column
nelements = 10
zmin = FT(-2)
zmax = FT(0)
dz_bottom = FT(1.0)
dz_top = FT(0.04)

domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)
