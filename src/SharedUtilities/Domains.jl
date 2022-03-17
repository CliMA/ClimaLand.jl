### Where will Domains.jl live eventually? Perhaps with Models.jl, as standalone
### models need to create domains of their own, and we dont want a copy of Domains
### in each component repo? Or maybe there does not need to be shared functionality,
### and each component can define their own domain as needed?
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
    nelements::Int32
    "Boundary face identifiers"
    boundary_tags::Tuple{Symbol, Symbol}
end

Base.ndims(::Column) = 1

Base.length(domain::Column) = domain.zlim[2] - domain.zlim[1]

Base.size(domain::Column) = length(domain)

"""
    function Column(FT::DataType = Float64; zlim, nelements)

Outer constructor for the `Column` type.
The `boundary_tags` field values are used to label the boundary faces 
at the top and bottom of the domain.
"""
function Column(FT::DataType = Float64; zlim, nelements)
    @assert zlim[1] < zlim[2]
    boundary_tags = (:bottom, :top)
    return Column{FT}(zlim, nelements, boundary_tags)
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
    mesh = Meshes.IntervalMesh(column; nelems = domain.nelements)
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
$(DocStringExtensions.FIELDS)
"""
struct Plane{FT} <: AbstractDomain{FT}
    xlim::Tuple{FT, FT}
    ylim::Tuple{FT, FT}
    nelements::Tuple{Int, Int}
    periodic::Tuple{Bool, Bool}
    npolynomial::Int
end

function Plane(
    FT::DataType = Float64;
    xlim,
    ylim,
    nelements,
    periodic,
    npolynomial,
)
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



# # # 3D hybrid domain
struct HybridBox{FT} <: AbstractDomain{FT}
    xlim::Tuple{FT, FT}
    ylim::Tuple{FT, FT}
    zlim::Tuple{FT, FT}
    nelements::Tuple{Int, Int, Int}
    npolynomial::Int
    periodic::Tuple{Bool, Bool}
end
"""
    HybridBox([FT = Float64]; xlim, ylim, zlim, nelements, npolynomial, periodic = (true, true).

Construct a domain of type `FT` that represents an xz-plane with limits `xlim` `ylim`
and `zlim` (where `xlim[1] < xlim[2]`,`ylim[1] < ylim[2]`, and `zlim[1] < zlim[2]`), `nelements`
elements of polynomial order `npolynomial`, (x,y)-axis periodicity = (true, true).

`nelements` must be a tuple with two values, with the first value corresponding
to the x-axis, the second corresponding to the y-axis, and the third corresponding to the z-axis. 

This domain is not periodic along the z-axis. Note that only periodic domains are supported
in the horizontal.
"""
function HybridBox(
    ::Type{FT} = Float64;
    xlim,
    ylim,
    zlim,
    nelements,
    npolynomial,
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
   coordinate(domain::Union{Column{FT}, Plane{FT}, HybridBox{FT}}) where {FT}

Returns the coordinate field for the domain.
"""
function coordinates(
    domain::Union{Column{FT}, Plane{FT}, HybridBox{FT}},
) where {FT}
    cs, _ = make_function_space(domain)
    cc = ClimaCore.Fields.coordinate_field(cs)
    return cc
end

export AbstractDomain, AbstractVegetationDomain
export Column, Plane, HybridBox, RootDomain
export coordinates

end
