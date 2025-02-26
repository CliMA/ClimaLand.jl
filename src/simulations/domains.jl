using ClimaLand
using ClimaUtilities.ClimaArtifacts
using ClimaComms
import Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
using ClimaCore: Spaces
"""
    global_domain(
    FT;
    apply_mask = true,
    mask_resolution = "60arcs",
    mask_threshold = 0.5,
    nelements = (101, 15),
    dz_tuple = (10.0, 0.05),
    depth = 50.0,
    npolynomial = 0,
)

Helper function to create a SphericalShell domain
with (101,15) elements, a depth of 50m, vertical
layering ranging from 0.05m in depth at the surface
to 10m at the bottom of the domain, with n_polynomial = 0
and GL quadrature.

The default continent dataset (60" resolution) is used with a threshold of 
50% land to create a land/sea mask.

`n_polynomial` determines the order of polynomial base to use for the spatial
discretization, which is correlated to the spatial resolution of the domain.

When `n_polynomial` is zero, the element is equivalent to a single point. In this
case, the resolution of the model is sqrt((360*180)/(101*101*6)). The factor of 6 arises
because there are 101x101 elements per side of the cubed sphere, meaning 6*101*101 for the
entire globe. 

When `n_polynomial` is greater than 1, a Gauss-Legendre-Lobotto quadrature
is used, with `n_polynomial + 1` points along the element. In this case, there are
always points two points on the boundaries for each direction with the other
points in the interior. These points are not equally spaced.

In practice, there is no reason to use `n_polynomial` greater than 1 in the current
version of ClimaLand. To increase resolution, we recommend increasing the number of elements
rather than increasing the polynomial order.
"""
function global_domain(
    FT;
    apply_mask = true,
    mask_resolution = "60arcs",
    mask_threshold = 0.5,
    nelements = (101, 15),
    dz_tuple = (10.0, 0.05),
    depth = 50.0,
    npolynomial = 0,
    context = nothing)
    radius = FT(6378.1e3)
    dz_tuple = FT.(dz_tuple)
    depth = FT(depth)
    domain = ClimaLand.Domains.SphericalShell(;
        radius,
        depth,
        nelements,
        npolynomial,
        dz_tuple,
    )
    if apply_mask
        surface_space = domain.space.surface # 2d space
        mask_path = ClimaLand.Artifacts.landseamask_file_path(;
            resolution = mask_resolution,
        )
        regridder_type = :InterpolationsRegridder
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        )
        mask = SpaceVaryingInput(
            mask_path,
            "landsea",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        binary_mask = apply_threshold.(mask, mask_threshold)
        # Here we would apply the mask to the domain
    end

    return domain
end

apply_threshold(field, value) =
    field > value ? eltype(field)(1) : eltype(field)(0)
