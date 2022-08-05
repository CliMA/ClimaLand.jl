### Long term plan: many of these domains will live in Models.jl repo, with
### the models.jl code as well. But LSM domains will live in ClimaLSM.
module Domains
using ClimaCore
using IntervalSets
using DocStringExtensions

import ClimaCore: Meshes, Spaces, Topologies, Geometry
import ClimaCore.Meshes: Uniform
### General type and methods all model domains will need

"""
    AbstractDomain{FT <:AbstractFloat}
An abstract type for domains.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end
Base.eltype(::AbstractDomain{FT}) where {FT} = FT
"""
    coordinates(domain::AbstractDomain)
Returns the coordinate field for the domain.
"""
function coordinates(domain::AbstractDomain)
    cc = ClimaCore.Fields.coordinate_field(domain.space)
    return cc
end

"""
    Point{FT} <: AbstractDomain{FT}
A domain for single column surface variables.
For models such as ponds, snow, roots, etc. Enables consistency 
in variable initialization across all domains.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Point{FT, S} <: AbstractDomain{FT}
    "Surface elevation relative to a reference (m)"
    z_sfc::FT
    "The associated ClimaCore Space"
    space::S
end

"""
    Point(;z_sfc::FT) where {FT}
Constructor for the `Point` domain using keyword arguments.
"""
function Point(; z_sfc::FT) where {FT}
    coord = ClimaCore.Geometry.ZPoint(z_sfc)
    space = ClimaCore.Spaces.PointSpace(coord)
    return Point{FT, typeof(space)}(z_sfc, space)
end

"""
    Column{FT} <: AbstractDomain{FT}
A struct holding the necessary information 
to construct a domain, a mesh, a center and face
space, etc. for use when a finite difference in
1D is suitable, as for a soil column model.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Column{FT, S} <: AbstractDomain{FT}
    "Domain interval limits, (zmin, zmax), in meters"
    zlim::Tuple{FT, FT}
    "Number of elements used to discretize the interval"
    nelements::Tuple{Int}
    "Boundary face identifiers"
    boundary_tags::Tuple{Symbol, Symbol}
    "The associated ClimaCore Space"
    space::S
end

"""
    Column(;
        zlim::Tuple{FT, FT},
        nelements::Int) where {FT}
Outer constructor for the `Column` type.
The `boundary_tags` field values are used to label the boundary faces 
at the top and bottom of the domain.
"""
function Column(; zlim::Tuple{FT, FT}, nelements::Int) where {FT}
    @assert zlim[1] < zlim[2]
    boundary_tags = (:bottom, :top)
    column = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint{FT}(zlim[1]),
        ClimaCore.Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = boundary_tags,
    )
    mesh = ClimaCore.Meshes.IntervalMesh(column; nelems = nelements)
    center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(mesh)
    return Column{FT, typeof(center_space)}(
        zlim,
        (nelements,),
        boundary_tags,
        center_space,
    )
end

"""
    Plane{FT} <: AbstractDomain{FT}
A struct holding the necessary information 
to construct a domain, a mesh, a 2d spectral
element space, and the resulting coordinate field.
Note that only periodic domains are currently supported.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Plane{FT, S} <: AbstractDomain{FT}
    "Domain interval limits along x axis, in meters"
    xlim::Tuple{FT, FT}
    "Domain interval limits along y axis, in meters"
    ylim::Tuple{FT, FT}
    "Number of elements to discretize interval, (nx, ny)"
    nelements::Tuple{Int, Int}
    "Flags for periodic boundaries; only true is supported"
    periodic::Tuple{Bool, Bool}
    "Polynomial order for both x and y"
    npolynomial::Int
    "The associated ClimaCore Space"
    space::S
end

"""
    Plane(;
        xlim::Tuple{FT,FT},
        ylim::Tuple{FT,FT},
        nelements::Tuple{Int,Int},
        periodic::Tuple{Bool,Bool},
        npolynomial::Int
        ) where {FT}
Outer constructor for the `Plane` domain, using keyword arguments.
"""
function Plane(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    nelements::Tuple{Int, Int},
    periodic::Tuple{Bool, Bool} = (true, true),
    npolynomial::Int,
) where {FT}
    @assert xlim[1] < xlim[2]
    @assert ylim[1] < ylim[2]
    @assert periodic == (true, true)
    domain_x = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.XPoint(xlim[1]),
        ClimaCore.Geometry.XPoint(xlim[2]);
        periodic = periodic[1],
    )
    domain_y = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.YPoint(ylim[1]),
        ClimaCore.Geometry.YPoint(ylim[2]);
        periodic = periodic[2],
    )
    plane = ClimaCore.Domains.RectangleDomain(domain_x, domain_y)

    mesh = ClimaCore.Meshes.RectilinearMesh(plane, nelements[1], nelements[2])
    grid_topology = ClimaCore.Topologies.Topology2D(mesh)
    if npolynomial == 0
        quad = ClimaCore.Spaces.Quadratures.GL{npolynomial + 1}()
    else
        quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    end
    space = ClimaCore.Spaces.SpectralElementSpace2D(grid_topology, quad)
    return Plane{FT, typeof(space)}(
        xlim,
        ylim,
        nelements,
        periodic,
        npolynomial,
        space,
    )
end

"""
    struct HybridBox{FT} <: AbstractDomain{FT}
        xlim::Tuple{FT, FT}
        ylim::Tuple{FT, FT}
        zlim::Tuple{FT, FT}
        nelements::Tuple{Int, Int, Int}
        npolynomial::Int
        periodic::Tuple{Bool, Bool}
    end
A struct holding the necessary information to construct a domain, a mesh, 
a 2d spectral element space (horizontal) x a 1d finite difference space
 (vertical), and the resulting coordinate field.
This domain is not periodic along the z-axis. Note that 
only periodic domains are supported
in the horizontal.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct HybridBox{FT, S} <: AbstractDomain{FT}
    "Domain interval limits along x axis, in meters"
    xlim::Tuple{FT, FT}
    "Domain interval limits along y axis, in meters"
    ylim::Tuple{FT, FT}
    "Domain interval limits along z axis, in meters"
    zlim::Tuple{FT, FT}
    "Number of elements to discretize interval, (nx, ny,nz)"
    nelements::Tuple{Int, Int, Int}
    " Polynomial order for the horizontal directions"
    npolynomial::Int
    "Flag indicating periodic boundaries in horizontal. only true is supported"
    periodic::Tuple{Bool, Bool}
    "The associated ClimaCore Space"
    space::S
end

"""
    HybridBox(;
        xlim::Tuple{FT, FT},
        ylim::Tuple{FT, FT},
        zlim::Tuple{FT, FT},
        nelements::Tuple{Int, Int, Int},
        npolynomial::Int,
        periodic = (true, true),
    ) where {FT}
Constructs the `HybridBox` domain
 with limits `xlim` `ylim` and `zlim` 
(where `xlim[1] < xlim[2]`,`ylim[1] < ylim[2]`, and `zlim[1] < zlim[2]`), 
`nelements` must be a tuple with three values, with the first 
value corresponding
to the x-axis, the second corresponding to the y-axis, and the third 
corresponding to the z-axis. The domain is periodic at the (xy) boundaries,
and the function space is of polynomial order `npolynomial` in the
horizontal directions.
"""
function HybridBox(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    zlim::Tuple{FT, FT},
    nelements::Tuple{Int, Int, Int},
    npolynomial::Int,
    periodic = (true, true),
) where {FT}
    @assert xlim[1] < xlim[2]
    @assert ylim[1] < ylim[2]
    @assert zlim[1] < zlim[2]
    @assert periodic == (true, true)
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(zlim[1]),
        ClimaCore.Geometry.ZPoint(zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = ClimaCore.Meshes.IntervalMesh(vertdomain, nelems = nelements[3])
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = Plane(;
        xlim = xlim,
        ylim = ylim,
        nelements = nelements[1:2],
        periodic = periodic,
        npolynomial = npolynomial,
    )
    horzspace = horzdomain.space

    hv_center_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
    return HybridBox{FT, typeof(hv_center_space)}(
        xlim,
        ylim,
        zlim,
        nelements,
        npolynomial,
        periodic,
        hv_center_space,
    )
end

"""
    struct SphericalShell{FT} <: AbstractDomain{FT}
        radius::FT
        height::FT
        nelements::Tuple{Int, Int}
        npolynomial::Int
    end
A struct holding the necessary information to construct a domain, a mesh, 
a 2d spectral element space (non-radial directions) 
x a 1d finite difference space (radial direction),
 and the resulting coordinate field.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalShell{FT, S} <: AbstractDomain{FT}
    "The radius of the shell"
    radius::FT
    "The radial extent of the shell"
    height::FT
    "The number of elements to be used in the non-radial and radial directions"
    nelements::Tuple{Int, Int}
    "The polynomial order to be used in the non-radial directions"
    npolynomial::Int
    "The associated ClimaCore Space"
    space::S
end

"""
    SphericalShell(;
        radius::FT,
        height::FT,
        nelements::Tuple{Int, Int},
        npolynomial::Int,
    ) where {FT}
Outer constructor for the `SphericalShell` domain, using keyword arguments.
"""
function SphericalShell(;
    radius::FT,
    height::FT,
    nelements::Tuple{Int, Int},
    npolynomial::Int,
) where {FT}
    @assert 0 < radius
    @assert 0 < height
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(0)),
        ClimaCore.Geometry.ZPoint(FT(height));
        boundary_tags = (:bottom, :top),
    )

    vertmesh = ClimaCore.Meshes.IntervalMesh(
        vertdomain,
        ClimaCore.Meshes.Uniform(),
        nelems = nelements[2],
    )
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements[1])
    horztopology = ClimaCore.Topologies.Topology2D(horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    hv_center_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
    return SphericalShell{FT, typeof(hv_center_space)}(
        radius,
        height,
        nelements,
        npolynomial,
        hv_center_space,
    )
end


"""
    struct SphericalSurface{FT} <: AbstractDomain{FT}
        radius::FT
        nelements::Tuple{Int, Int}
        npolynomial::Int
    end
A struct holding the necessary information to construct a domain, a mesh, 
a 2d spectral element space (non-radial directions) and the resulting coordinate field.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalSurface{FT, S} <: AbstractDomain{FT}
    "The radius of the surface"
    radius::FT
    "The number of elements to be used in the non-radial directions"
    nelements::Int
    "The polynomial order to be used in the non-radial directions"
    npolynomial::Int
    "The associated ClimaCore Space"
    space::S
end

"""
    SphericalSurface(;
        radius::FT,
        nelements::Int
        npolynomial::Int,
    ) where {FT}
Outer constructor for the `SphericalSurface` domain, using keyword arguments.
"""
function SphericalSurface(;
    radius::FT,
    nelements::Int,
    npolynomial::Int,
) where {FT}
    @assert 0 < radius
    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = Meshes.EquiangularCubedSphere(horzdomain, nelements)
    horztopology = Topologies.Topology2D(horzmesh)
    quad = Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = Spaces.SpectralElementSpace2D(horztopology, quad)
    return SphericalSurface{FT, typeof(horzspace)}(
        radius,
        nelements,
        npolynomial,
        horzspace,
    )
end


abstract type AbstractVegetationDomain{FT} <: AbstractDomain{FT} end
"""
   RootDomain{FT} <: AbstractVegetationDomain{FT}
Domain for a single bulk plant with roots ofvarying depths. The user needs
to specify the depths of the root tips as wel as the heights of the
compartments to be modeled within the plant. The compartment heights
are expected to be sorted in ascending order.
"""
struct RootDomain{FT} <: AbstractVegetationDomain{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "The height of the stem, leaf compartments, in meters"
    compartment_heights::Vector{FT}
end

function coordinates(domain::RootDomain{FT}) where {FT}
    return domain.compartment_heights
end

"""
    AbstractLSMDomain{FT} <: AbstractDomain{FT}
An abstract type for LSMDomains, which have two components: a surface
and a subsurface.
"""
abstract type AbstractLSMDomain{FT} <: AbstractDomain{FT} end

"""
    LSMSingleColumnDomain{FT} <: AbstractLSMDomain{FT}
A mixed domain, consisting of a column domain with z-coordinates at the
finite difference cell centers, and a point domain, with a single z
coordinate at the top boundary of the column domain.
For use in LSM modeling, where a subsurface finite difference space 
(for modeling soil hydrology and energy) and a surface space are both
needed.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LSMSingleColumnDomain{FT, S1, S2} <: AbstractLSMDomain{FT}
    "The subsurface Column domain"
    subsurface::Column{FT, S1}
    "The surface Point domain"
    surface::Point{FT, S2}
end

"""
    LSMSingleColumnDomain(;
        zlim::Tuple{FT, FT},
        nelements::Int,
    ) where {FT}
A constructor for the LSMSingleColumnDomain.
"""
function LSMSingleColumnDomain(; zlim::Tuple{FT, FT}, nelements::Int) where {FT}
    @assert zlim[1] < zlim[2]
    subsurface_domain = Column(; zlim = FT.(zlim), nelements = nelements)
    surface_space = obtain_surface_space(subsurface_domain.space)
    surface_domain = Point{FT, typeof(surface_space)}(zlim[2], surface_space)
    return LSMSingleColumnDomain{
        FT,
        typeof(subsurface_domain.space),
        typeof(surface_space),
    }(
        subsurface_domain,
        surface_domain,
    )
end


struct LSMMultiColumnDomain{FT, S1, S2} <: AbstractLSMDomain{FT}
    "The subsurface Column domain"
    subsurface::HybridBox{FT, S1}
    "The surface Point domain"
    surface::Plane{FT, S2}
end

"""
    LSMMultiColumnDomain(;
        xlim::Tuple{FT, FT}
        ylim::Tuple{FT, FT}
        zlim::Tuple{FT, FT}
        nelements::Tuple{Int, Int, Int}
        npolynomial::Int
        periodic::Tuple{Bool, Bool}
    ) where {FT}
A constructor for the LSMSingleColumnDomain.
"""
function LSMMultiColumnDomain(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    zlim::Tuple{FT, FT},
    nelements::Tuple{Int, Int, Int},
    npolynomial::Int,
    periodic::Tuple{Bool, Bool},
) where {FT}
    @assert xlim[1] < xlim[2]
    @assert ylim[1] < ylim[2]
    @assert zlim[1] < zlim[2]
    @assert periodic == (true, true)
    subsurface = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        nelements = nelements,
        npolynomial = npolynomial,
        periodic = periodic,
    )
    surface_space = obtain_surface_space(subsurface.space)
    surface = Plane{FT, typeof(surface_space)}(
        xlim,
        ylim,
        nelements[1:2],
        periodic,
        npolynomial,
        surface_space,
    )
    return LSMMultiColumnDomain{
        FT,
        typeof(subsurface.space),
        typeof(surface.space),
    }(
        subsurface,
        surface,
    )
end

"""
    coordinates(domain::AbstractLSMDomain{FT}) where {FT}
Returns the coordinates of the AbstractLSMDomain as a named tuple,
with keys of `subsurface` and `surface`.
"""
function coordinates(domain::AbstractLSMDomain{FT}) where {FT}
    return (
        subsurface = coordinates(domain.subsurface),
        surface = coordinates(domain.surface),
    )
end

"""
    obtain_face_space(cs::ClimaCore.Spaces.AbstractSpace)
Returns the face space, if applicable, for the center space `cs`.
"""
obtain_face_space(cs::ClimaCore.Spaces.AbstractSpace) =
    @error("No face space is defined for this space.")

"""
    obtain_surface_space(cs::ClimaCore.Spaces.AbstractSpace)
Returns the surface space, if applicable, for the center space `cs`.
"""
obtain_surface_space(cs::ClimaCore.Spaces.AbstractSpace) =
    @error("No surface space is defined for this space.")

"""
    obtain_face_space(cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace)
Returns the face space for the CenterExtrudedFiniteDifferenceSpace `cs`.
"""
function obtain_face_space(
    cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
)
    return ClimaCore.Spaces.FaceExtrudedFiniteDifferenceSpace(cs)
end

"""
    obtain_face_space(cs::ClimaCore.Spaces.CenterFiniteDifferenceSpace)
Returns the face space corresponding to the CenterFiniteDifferenceSpace `cs`.
"""
function obtain_face_space(cs::ClimaCore.Spaces.CenterFiniteDifferenceSpace)
    return ClimaCore.Spaces.FaceFiniteDifferenceSpace(cs)
end

"""
    obtain_surface_space(cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace)
Returns the horizontal space for the CenterExtrudedFiniteDifferenceSpace `cs`.
"""
function obtain_surface_space(
    cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
)
    return cs.horizontal_space
end

"""
    obtain_surface_space(cs::ClimaCore.Spaces.CenterFiniteDifferenceSpace)
Returns the top level of the face space corresponding to the CenterFiniteDifferenceSpace `cs`.
"""
function obtain_surface_space(cs::ClimaCore.Spaces.CenterFiniteDifferenceSpace)
    fs = obtain_face_space(cs)
    return ClimaCore.Spaces.level(
        fs,
        ClimaCore.Utilities.PlusHalf(ClimaCore.Spaces.nlevels(fs) - 1),
    )
end




export AbstractDomain, AbstractVegetationDomain, AbstractLSMDomain
export Column,
    Plane, HybridBox, RootDomain, Point, SphericalShell, SphericalSurface
export LSMSingleColumnDomain, LSMMultiColumnDomain
export coordinates, obtain_face_space, obtain_surface_space

end
