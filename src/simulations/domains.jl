"""
    global_domain(
    FT;
    nelements = (101, 15),
    dz_tuple = (10.0, 0.05),
    depth = 50.0,
    npolynomial = 0,
)

Helper function to create a SphericalShell domain with (101,15) elements, a
depth of 50m, vertical layering ranging from 0.05m in depth at the surface to
10m at the bottom of the domain, with npolynomial = 0 and GL quadrature.

`npolynomial` determines the order of polynomial base to use for the spatial
discretization, which is correlated to the spatial resolution of the domain.

When `npolynomial` is zero, the element is equivalent to a single point. In this
case, the resolution of the model is sqrt((360*180)/(101*101*6)). The factor of 6 arises
because there are 101x101 elements per side of the cubed sphere, meaning 6*101*101 for the
entire globe.

When `npolynomial` is greater than 1, a Gauss-Legendre-Lobotto quadrature is
used, with `npolynomial + 1` points along the element. In this case, there are
always points two points on the boundaries for each direction with the other
points in the interior. These points are not equally spaced.

In practice, there is no reason to use `npolynomial` greater than 1 in the current
version of ClimaLand. To increase resolution, we recommend increasing the number of elements
rather than increasing the polynomial order.
"""
function global_domain(
    FT;
    apply_mask = true,
    mask_threshold = 0.5,
    nelements = (101, 15),
    dz_tuple = (10.0, 0.05),
    depth = 50.0,
    npolynomial = 0,
    context = ClimaComms.context(),
)
    if pkgversion(ClimaCore) < v"0.14.30" && apply_mask
        @warn "The land mask cannot be applied with ClimaCore < v0.14.30. Update ClimaCore for significant performance gains."
        apply_mask = false
    end

    radius = FT(6378.1e3)
    dz_tuple = FT.(dz_tuple)
    depth = FT(depth)
    domain = ClimaLand.Domains.SphericalShell(;
        radius,
        depth,
        nelements,
        npolynomial,
        dz_tuple,
        context,
    )
    if apply_mask
        surface_space = domain.space.surface # 2d space
        binary_mask = ClimaLand.landsea_mask(domain; threshold = mask_threshold)
        ClimaCore.Spaces.set_mask!(surface_space, binary_mask)
    end

    return domain
end
