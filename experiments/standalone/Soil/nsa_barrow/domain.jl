# Domain setup
# For soil column
nelements = 20
zmin = FT(-10)
zmax = FT(0)
dz_bottom = FT(1.0)
dz_top = FT(0.05)

domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)
