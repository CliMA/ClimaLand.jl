module Domains
using ClimaCore
using DocStringExtensions

import ClimaCore: Meshes, Spaces, Topologies, Geometry

### General type and methods all model domains will need

"""
    AbstractDomain{FT <:AbstractFloat}

An abstract type for domains.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end

"""
    coordinates(domain::AbstractDomain)

Method which returns the coordinates appropriate for a given domain.
"""
function coordinates(domain::AbstractDomain) end

### Example of component specific domain
"""
    AbstractPlantDomain{FT} <: AbstractDomain{FT}

An abstract type for vegetation specific domains.
"""
abstract type AbstractPlantDomain{FT} <: AbstractDomain{FT} end


"""
   RootDomain{FT} <: AbstractPlantDomain{FT}

Domain for a single bulk plant with roots of varying depths. The user needs
to specify the depths of the root tips as wel as the heights of the
compartments to be modeled within the plant. The compartment heights
are expected to be sorted in ascending order.
"""
struct RootDomain{FT} <: AbstractPlantDomain{FT}
    "The depth of the root tips, in meters"
    root_depths::Vector{FT}
    "The height of the stem, leaf compartments, in meters"
    compartment_heights::Vector{FT}
end

function coordinates(domain::RootDomain{FT}) where {FT}
    return domain.compartment_heights
end


### Another example, but a ClimaCore FD Column space/domain may be general enough
### that we keep it as a concrete type of AbstractDomain, rather than an AbstractSoilDomain

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
    coordinates(domain::Column{FT}) where {FT}
        
Get the center coordinates from the `Column` domain. Face coordinates will be
added if necessary. Note that this function can only be called one time
for coordinates derived from ClimaCore.Spaces.

This is a required function for each concrete type of domain.
    """
function coordinates(domain::Column{FT}) where {FT}
    cs, _ = make_function_space(domain)
    cc = ClimaCore.Fields.coordinate_field(cs).z
    return cc
end

export AbstractDomain, AbstractPlantDomain
export Column, RootDomain
export coordinates


end
