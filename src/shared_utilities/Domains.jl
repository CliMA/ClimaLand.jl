### Long term plan: many of these domains will live in Models.jl repo, with
### the models.jl code as well. But LSM domains will live in ClimaLSM.
module Domains
using ClimaCore
using ClimaComms
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
For models such as ponds, snow, plant hydraulics, etc. Enables consistency 
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
    Point(;z_sfc::FT,
           comms = ClimaComms.SingletonCommsContext()
          ) where {FT}

Constructor for the `Point` domain using keyword arguments.

All other ClimaLSM domains rely on default comms set internally
by ClimaCore. However, the Point space is unique in this context,
and does not have the same default defined in ClimaCore. 
Because of this, we set the default here
in ClimaLSM. In long term, we will repeat the same for all ClimaLSM domains
and not rely on any internal defaults set in ClimaCore.
"""
function Point(;
    z_sfc::FT,
    comms = ClimaComms.SingletonCommsContext(),
) where {FT}
    coord = ClimaCore.Geometry.ZPoint(z_sfc)
    space = ClimaCore.Spaces.PointSpace(comms, coord)
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
    "Tuple for mesh stretching specifying *target* (dz_bottom, dz_top) (m). If nothing, no stretching is applied."
    dz_tuple::Union{Tuple{FT, FT}, Nothing}
    "Boundary face identifiers"
    boundary_tags::Tuple{Symbol, Symbol}
    "The associated ClimaCore Space"
    space::S
end

"""
    Column(;
        zlim::Tuple{FT, FT},
        nelements::Int,
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing) where {FT}

Outer constructor for the `Column` type. 

Using `ClimaCore` tools, the coordinate
mesh can be stretched such that the top of the domain has finer resolution
than the bottom of the domain. In order to activate this, a tuple with the
target dz_bottom and dz_top should be passed via keyword argument. The default is
uniform spacing. Please note that in order to use this feature, ClimaCore requires that
the elements of zlim be <=0. Additionally, the dz_tuple you supply may not be compatible
with the domain boundaries in some cases, in which case you may need to choose
different values.

The `boundary_tags` field values are used to label the boundary faces 
at the top and bottom of the domain.
"""
function Column(;
    zlim::Tuple{FT, FT},
    nelements::Int,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert zlim[1] < zlim[2]
    boundary_tags = (:bottom, :top)
    column = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint{FT}(zlim[1]),
        ClimaCore.Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = boundary_tags,
    )
    if dz_tuple isa Nothing
        mesh = ClimaCore.Meshes.IntervalMesh(column; nelems = nelements)
    else
        @assert zlim[2] <= 0
        mesh = ClimaCore.Meshes.IntervalMesh(
            column,
            ClimaCore.Meshes.GeneralizedExponentialStretching{FT}(
                dz_tuple[1],
                dz_tuple[2],
            );
            nelems = nelements,
            reverse_mode = true,
        )
    end

    center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(mesh)
    return Column{FT, typeof(center_space)}(
        zlim,
        (nelements,),
        dz_tuple,
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
        dz_tuple::Union{Tuple{FT, FT}, Nothing}
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
    "Tuple for mesh stretching specifying *target* (dz_bottom, dz_top) (m). If nothing, no stretching is applied."
    dz_tuple::Union{Tuple{FT, FT}, Nothing}
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
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
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

Using `ClimaCore` tools, the coordinate
mesh can be stretched such that the top of the domain has finer resolution
than the bottom of the domain. In order to activate this, a tuple with the
target dz_bottom and dz_top should be passed via keyword argument. The default is
uniform spacing. Please note that in order to use this feature, ClimaCore requires that
the elements of zlim be <=0. Additionally, the dz_tuple you supply may not be compatible
with the domain boundaries in some cases, in which case you may need to choose
different values.
"""
function HybridBox(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    zlim::Tuple{FT, FT},
    nelements::Tuple{Int, Int, Int},
    npolynomial::Int,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
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
    if dz_tuple isa Nothing
        vertmesh =
            ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = nelements[3])
    else
        @assert zlim[2] <= 0
        vertmesh = ClimaCore.Meshes.IntervalMesh(
            vertdomain,
            ClimaCore.Meshes.GeneralizedExponentialStretching{FT}(
                dz_tuple[1],
                dz_tuple[2],
            );
            nelems = nelements[3],
            reverse_mode = true,
        )
    end
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
        dz_tuple,
        nelements,
        npolynomial,
        periodic,
        hv_center_space,
    )
end

"""
    struct SphericalShell{FT} <: AbstractDomain{FT}
        radius::FT
        depth::FT
        dz_tuple::Union{Tuple{FT, FT}, Nothing}
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
    depth::FT
    "Tuple for mesh stretching specifying *target* (dz_bottom, dz_top) (m). If nothing, no stretching is applied."
    dz_tuple::Union{Tuple{FT, FT}, Nothing}
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
        depth::FT,
        nelements::Tuple{Int, Int},
        npolynomial::Int,
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    ) where {FT}
Outer constructor for the `SphericalShell` domain, using keyword arguments.

Using `ClimaCore` tools, the coordinate
mesh can be stretched such that the top of the domain has finer resolution
than the bottom of the domain. In order to activate this, a tuple with the
target dz_bottom and dz_top should be passed via keyword argument. The default is
uniform spacing. Please note that the dz_tuple you supply may not be compatible
with the depth/nelements chosen, in which case you may need to choose
different values.
"""
function SphericalShell(;
    radius::FT,
    depth::FT,
    nelements::Tuple{Int, Int},
    npolynomial::Int,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert 0 < radius
    @assert 0 < depth
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(-depth)),
        ClimaCore.Geometry.ZPoint(FT(0));
        boundary_tags = (:bottom, :top),
    )
    if dz_tuple isa Nothing
        vertmesh =
            ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = nelements[2])
    else
        vertmesh = ClimaCore.Meshes.IntervalMesh(
            vertdomain,
            ClimaCore.Meshes.GeneralizedExponentialStretching{FT}(
                dz_tuple[1],
                dz_tuple[2],
            );
            nelems = nelements[2],
            reverse_mode = true,
        )
    end
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
        depth,
        dz_tuple,
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
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    ) where {FT}
A constructor for the LSMSingleColumnDomain.
"""
function LSMSingleColumnDomain(;
    zlim::Tuple{FT, FT},
    nelements::Int,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert zlim[1] < zlim[2]
    subsurface_domain =
        Column(; zlim = FT.(zlim), nelements = nelements, dz_tuple = dz_tuple)
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

"""
    LSMMultiColumnDomain{FT} <: AbstractLSMDomain{FT}

A mixed domain, consisting of a hybdrid box domain with z-coordinates at the
finite difference cell centers, and a plane domain, coinciding with the surface
of the box.

For use in LSM modeling, where a subsurface finite difference space 
(for modeling soil hydrology and energy) and a surface space are both
needed.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LSMMultiColumnDomain{FT, S1, S2} <: AbstractLSMDomain{FT}
    "The subsurface box domain"
    subsurface::HybridBox{FT, S1}
    "The surface plane domain"
    surface::Plane{FT, S2}
end

"""
    LSMMultiColumnDomain(;
        xlim::Tuple{FT, FT}
        ylim::Tuple{FT, FT}
        zlim::Tuple{FT, FT}
        nelements::Tuple{Int, Int, Int}
        npolynomial::Int
        periodic::Tuple{Bool, Bool} = (true, true),
        dz_tuple::Union{Tuple{FT,FT}, Nothing} = nothing,
    ) where {FT}
A constructor for the LSMMultiColumnDomain.
"""
function LSMMultiColumnDomain(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    zlim::Tuple{FT, FT},
    nelements::Tuple{Int, Int, Int},
    npolynomial::Int,
    periodic::Tuple{Bool, Bool} = (true, true),
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert xlim[1] < xlim[2]
    @assert ylim[1] < ylim[2]
    @assert zlim[1] < zlim[2]
    @assert periodic == (true, true)
    subsurface = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        dz_tuple = dz_tuple,
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
    LSMSphericalShellDomain{FT} <: AbstractLSMDomain{FT}

A mixed domain, consisting of a spherical shell domain with z-coordinates at the
finite difference cell centers, and a spherical surface domain, coinciding with 
the surface of the shell.

For use in LSM modeling, where a subsurface finite difference space 
(for modeling soil hydrology and energy) and a surface space are both
needed.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LSMSphericalShellDomain{FT, S1, S2} <: AbstractLSMDomain{FT}
    "The subsurface shell domain"
    subsurface::SphericalShell{FT, S1}
    "The surface domain"
    surface::SphericalSurface{FT, S2}
end

"""
    LSMSphericalShellDomain(;
        radius::FT,
        depth::FT,
        nelements::Tuple{Int, Int},
        npolynomial::Int,
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    ) where {FT}
A constructor for the LSMSphericalShellDomain.
"""
function LSMSphericalShellDomain(;
    radius::FT,
    depth::FT,
    nelements::Tuple{Int, Int},
    npolynomial::Int,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert radius > FT(0)
    @assert depth > FT(0)
    subsurface = SphericalShell(;
        radius = radius,
        depth = depth,
        dz_tuple = dz_tuple,
        nelements = nelements,
        npolynomial = npolynomial,
    )
    surface_space = obtain_surface_space(subsurface.space)
    surface = SphericalSurface{FT, typeof(surface_space)}(
        radius,
        nelements[1],
        npolynomial,
        surface_space,
    )
    return LSMSphericalShellDomain{
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

"""
    top_center_to_surface(center_field::ClimaCore.Fields.Field)

Creates and returns a ClimaCore.Fields.Field defined on the space
corresponding to the surface of the space on which `center_field` 
is defined, with values equal to the those at the level of the top 
center.

For example, given a `center_field` defined on 1D center finite difference space,
this would return a field defined on the Point space of the surface of
the column. The value would be the value of the oroginal `center_field` 
at the topmost location. Given a `center_field` defined on a 3D
extruded center finite difference space, this would return a 2D field
corresponding to the surface, with values equal to the topmost level.
"""
function top_center_to_surface(center_field::ClimaCore.Fields.Field)
    cs = axes(center_field)
    face_space = obtain_face_space(cs)
    N = ClimaCore.Spaces.nlevels(face_space)
    interp_c2f = ClimaCore.Operators.InterpolateC2F(
        top = ClimaCore.Operators.Extrapolate(),
        bottom = ClimaCore.Operators.Extrapolate(),
    )
    surface_space = obtain_surface_space(cs)
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaCore.Fields.level(
                interp_c2f.(center_field),
                ClimaCore.Utilities.PlusHalf(N - 1),
            ),
        ),
        surface_space,
    )
end


export AbstractDomain, AbstractLSMDomain
export Column, Plane, HybridBox, Point, SphericalShell, SphericalSurface
export LSMSingleColumnDomain, LSMMultiColumnDomain, LSMSphericalShellDomain
export coordinates,
    obtain_face_space, obtain_surface_space, top_center_to_surface

end
