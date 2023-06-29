# Domain setup
# For soil column
nelements = 10
zmin = FT(-2)
zmax = FT(0)
land_domain =
    LSMSingleColumnDomain(; zlim = (zmin, zmax), nelements = nelements)

# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
n_stem = Int64(1)
n_leaf = Int64(1)
h_stem = FT(9) # m, from Wang et al.
h_leaf = FT(9.5) # m from Wang et al.
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf]
