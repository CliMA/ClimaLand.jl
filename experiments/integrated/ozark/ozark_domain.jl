# Domain setup
# For soil column
nelements = 10
zmin = FT(-2)
zmax = FT(0)
dz_bottom = FT(0.5)
dz_top = FT(0.025)

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)
canopy_domain = ClimaLSM.Domains.obtain_surface_domain(land_domain)

# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m, from Wang et al.
h_leaf = FT(9.5) # m from Wang et al.
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf]
