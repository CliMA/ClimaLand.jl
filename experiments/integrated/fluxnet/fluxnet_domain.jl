"""This file sets up the domain to run Clima Land on a fluxtower site, using both
the generic parameters set in this file as well as the site-specific parameters
given in the {site-ID}_simulation.jl file in each site directory."""

# Domain setup
# For soil column
h_canopy = h_stem + h_leaf
compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]
