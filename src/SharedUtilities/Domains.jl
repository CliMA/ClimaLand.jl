### Long term plan: many of these domains will live in Models.jl repo, with
### the models.jl code as well. But LSM domains will live in ClimaLSM.
module Domains
using ClimaCore
using IntervalSets
using DocStringExtensions

import ClimaCore: Meshes, Spaces, Topologies, Geometry

### General type and methods all model domains will need

"""
    AbstractDomain{FT <:AbstractFloat}

An abstract type for domains.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end
Base.eltype(::AbstractDomain{FT}) where {FT} = FT
"""
    coordinates(domain::AbstractDomain)

Method which returns the coordinates appropriate for a given domain.

The coordinates can be Fields or Vectors.
"""
function coordinates(domain::AbstractDomain) end

"""
    Point{FT} <: AbstractDomain{FT}

A domain for single column surface variables.

For models such as ponds, snow, roots, etc. Enables consistency 
in variable initialization across all domains.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Point{FT} <: AbstractDomain{FT}
    "Surface elevation relative to a reference (m)"
    z_sfc::FT
end

"""
    Point(;z_sfc::FT) where {FT}

Constructor for the `Point` domain using keyword arguments.
"""
function Point(; z_sfc::FT) where {FT}
    return Point{FT}(z_sfc)
end


coordinates(domain::Point) = [domain.z_sfc]

"""
    Column{FT} <: AbstractDomain{FT}

A struct holding the necessary information 
to construct a domain, a mesh, a center and face
space, etc. for use when a finite difference in
1D is suitable, as for a soil column model.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Column{FT} <: AbstractDomain{FT}
    "Domain interval limits, (zmin, zmax), in meters"
    zlim::Tuple{FT, FT}
    "Number of elements used to discretize the interval"
    nelements::Tuple{Int}
    "Boundary face identifiers"
    boundary_tags::Tuple{Symbol, Symbol}
end

Base.ndims(::Column) = 1

Base.length(domain::Column) = domain.zlim[2] - domain.zlim[1]

Base.size(domain::Column) = length(domain)

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
    return Column{FT}(zlim, (nelements,), boundary_tags)
end

"""
    make_function_space(domain::Column)

Returns the center and face space of the column domain.
"""
function make_function_space(domain::Column{FT}) where {FT}
    column = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint{FT}(domain.zlim[1]),
        ClimaCore.Geometry.ZPoint{FT}(domain.zlim[2]);
        boundary_tags = domain.boundary_tags,
    )
    mesh = Meshes.IntervalMesh(column; nelems = domain.nelements[1])
    center_space = Spaces.CenterFiniteDifferenceSpace(mesh)
    face_space = Spaces.FaceFiniteDifferenceSpace(center_space)

    return center_space, face_space
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
struct Plane{FT} <: AbstractDomain{FT}
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
    return Plane{FT}(xlim, ylim, nelements, periodic, npolynomial)
end


"""
    make_function_space(domain::Plane)

Returns the 2d spectral element space of the
desired periodicity, nodal point type, and polynomial order.

Note that only periodic boundaries are supported.
"""
function make_function_space(domain::Plane{FT}) where {FT}
    domain_x = ClimaCore.Domains.IntervalDomain(
        Geometry.XPoint(domain.xlim[1]),
        Geometry.XPoint(domain.xlim[2]);
        periodic = domain.periodic[1],
    )
    domain_y = ClimaCore.Domains.IntervalDomain(
        Geometry.YPoint(domain.ylim[1]),
        Geometry.YPoint(domain.ylim[2]);
        periodic = domain.periodic[2],
    )
    plane = ClimaCore.Domains.RectangleDomain(domain_x, domain_y)

    mesh =
        Meshes.RectilinearMesh(plane, domain.nelements[1], domain.nelements[2])
    grid_topology = Topologies.Topology2D(mesh)
    if domain.npolynomial == 0
        quad = Spaces.Quadratures.GL{domain.npolynomial + 1}()
    else
        quad = Spaces.Quadratures.GLL{domain.npolynomial + 1}()
    end
    space = Spaces.SpectralElementSpace2D(grid_topology, quad)

    return space, nothing
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
struct HybridBox{FT} <: AbstractDomain{FT}
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
    return HybridBox{FT}(xlim, ylim, zlim, nelements, npolynomial, periodic)
end

"""
    make_function_space(domain::HybridBox)

Returns the extruded finite difference center and face
finite spaces of the
desired periodicity, nodal point type, and polynomial order
in the horizontal.

Note that only periodic boundaries are supported. 
"""
function make_function_space(domain::HybridBox{FT}) where {FT}
    vertdomain = ClimaCore.Domains.IntervalDomain(
        Geometry.ZPoint(domain.zlim[1]),
        Geometry.ZPoint(domain.zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = domain.nelements[3])
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = Plane{FT}(
        domain.xlim,
        domain.ylim,
        domain.nelements[1:2],
        domain.periodic,
        domain.npolynomial,
    )
    horzspace, _ = make_function_space(horzdomain)

    hv_center_space =
        Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)

    return hv_center_space, hv_face_space
end

"""
   coordinates(domain::Union{Column{FT}, Plane{FT}, HybridBox{FT}}) where {FT}

Returns the coordinate field for the domain.
"""
function coordinates(
    domain::Union{Column{FT}, Plane{FT}, HybridBox{FT}},
) where {FT}
    cs, _ = make_function_space(domain)
    cc = ClimaCore.Fields.coordinate_field(cs)
    return cc
end


### Example of component specific domain
"""
    AbstractVegetationDomain{FT} <: AbstractDomain{FT}

An abstract type for vegetation specific domains.
"""
abstract type AbstractVegetationDomain{FT} <: AbstractDomain{FT} end


"""
   RootDomain{FT} <: AbstractVegetationDomain{FT}

Domain for a single bulk plant with roots of varying depths. The user needs
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
    LSMSingleColumnDomain{FT} <: AbstractDomain{FT}

A mixed domain, consisting of a column domain with z-coordinates at the
finite difference cell centers, and a point domain, with a single z
coordinate at the top boundary of the column domain.

For use in LSM modeling, where a subsurface finite difference space 
(for modeling soil hydrology and energy) and a surface space are both
needed.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LSMSingleColumnDomain{FT} <: AbstractDomain{FT}
    "The subsurface Column domain"
    subsurface::Column{FT}
    "The surface Point domain"
    surface::Point{FT}
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
    surface_domain = Point{FT}(FT(zlim[2]))
    boundary_tags = (:bottom, :top)
    subsurface_domain = Column{FT}(FT.(zlim), (nelements,), boundary_tags)
    return LSMSingleColumnDomain{FT}(subsurface_domain, surface_domain)
end

"""
    coordinates(domain::LSMSingleColumnDomain{FT}) where {FT}

Returns the coordinates of the LSMSingleColumnDomain as a named tuple,
with keys of `subsurface` and `surface`.
"""
function coordinates(domain::LSMSingleColumnDomain{FT}) where {FT}
    return (
        subsurface = coordinates(domain.subsurface),
        surface = coordinates(domain.surface),
    )
end



export AbstractDomain, AbstractVegetationDomain
export Column, Plane, HybridBox, RootDomain, Point
export LSMSingleColumnDomain
export coordinates

end
