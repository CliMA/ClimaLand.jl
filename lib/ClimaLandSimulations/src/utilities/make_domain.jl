export make_domain

"""
    make_domain(setup;
        nelements = FT(20),
        zmin = FT(-10),
        zmax = FT(0)
        )

setup the soil and canopy domain of the simulation. Returns:
- h_canopy,
- compartment_midpoints, 
- compartment_surfaces,
- land_domain,
- canopy_domain,
"""
function make_domain(setup, FT; nelements = 20, zmin = FT(-10), zmax = FT(0))

    h_canopy = setup.h_stem + setup.h_leaf

    compartment_midpoints =
        setup.n_stem > 0 ? [setup.h_stem / 2, setup.h_stem + setup.h_leaf / 2] :
        [setup.h_leaf / 2]

    compartment_surfaces =
        setup.n_stem > 0 ? [zmax, setup.h_stem, h_canopy] : [zmax, setup.h_leaf]

    land_domain = Column(;
        zlim = (zmin, zmax),
        nelements = nelements,
        dz_tuple = (setup.dz_bottom, setup.dz_top),
    )

    canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

    return (
        nelements = nelements,
        zmin = zmin,
        zmax = zmax,
        h_canopy = h_canopy,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
        land_domain = land_domain,
        canopy_domain = canopy_domain,
    )
end
