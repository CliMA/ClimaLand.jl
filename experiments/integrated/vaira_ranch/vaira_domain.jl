# Domain setup
# For soil column
nelements = 20
zmin = FT(-10)
zmax = FT(0)
dz_bottom = FT(1.0)
dz_top = FT(0.05)

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)
canopy_domain = ClimaLSM.Domains.obtain_surface_domain(land_domain)

# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
n_stem = Int64(0)
n_leaf = Int64(1)
h_leaf = FT(0.7) # m
h_canopy = h_leaf
compartment_midpoints = [h_leaf / 2]
compartment_surfaces = [zmax, h_leaf]
