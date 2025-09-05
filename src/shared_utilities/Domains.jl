module Domains

import ..compat_add_mask, ..compat_set_mask!
import ..Artifacts.landseamask_file_path

using ClimaCore
using ClimaComms
using DocStringExtensions
using Interpolations
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import StaticArrays: SMatrix

import ClimaCore: Meshes, Spaces, Topologies, Geometry
import ClimaCore.Meshes: Uniform
### General type and methods all model domains will need

"""
    AbstractDomain{FT <:AbstractFloat}

An abstract type for domains.

The domain structs typically hold information
regarding the bounds of the domain, the boundary condition type
(periodic or not), and the spatial discretization.

Additionally, the domain struct holds the relevant spaces for that
domain. For example, a 3D domain holds the center space (in terms of
finite difference - the space corresponding to the centers of each
element), and the top face space where surface fluxes are computed.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end
Base.eltype(::AbstractDomain{FT}) where {FT} = FT

ClimaComms.context(domain::AbstractDomain) =
    ClimaComms.context(first(domain.space))
ClimaComms.device(domain::AbstractDomain) =
    ClimaComms.device(first(domain.space))

"""
    coordinates(domain::AbstractDomain)

Returns the coordinate fields for the domain as a NamedTuple.

The returned coordinates are stored with keys :surface, :subsurface, e.g.
as relevant for the domain.
"""
function coordinates(domain::AbstractDomain)
    domain_names = propertynames(domain.space)
    domain_spaces = getproperty.(Ref(domain.space), domain_names)
    return NamedTuple{domain_names}(
        ClimaCore.Fields.coordinate_field.(domain_spaces),
    )
end

"""
    Point{FT} <: AbstractDomain{FT}

A domain for single column surface variables.
For models such as ponds, snow, plant hydraulics, etc. Enables consistency
in variable initialization across all domains.

`space` is a NamedTuple holding the surface space (in this case,
the Point space).
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Point{FT, NT <: NamedTuple} <: AbstractDomain{FT}
    "Surface elevation relative to a reference (m)"
    z_sfc::FT
    "A NamedTuple of associated ClimaCore spaces: in this case, the Point (surface) space"
    space::NT
end

"""
    Point(; z_sfc::FT,
        longlat::Union{Tuple{FT, FT}, Nothing} = nothing,
        comms = ClimaComms.SingletonCommsContext()
        ) where {FT}

Constructor for the `Point` domain using keyword arguments.

If `longlat` is provided, the domain's space contains those coordinates.
This should be used for simulations at a point that require reading in spatial data
from a file defined on a lat/long grid, such as CLM data.
Note that constructing a Point with lat/long requires first constructing a 3D HybridBox
domain centered at `longlat`, then extracting a point from the intermediate 3D space.

The latitude and longitude of the returned domain can be accessed as follows:
- lat = ClimaLand.Domains.get_lat(domain.space.surface)
- long = ClimaLand.Domains.get_long(domain.space.surface)
"""
function Point(;
    z_sfc::FT,
    longlat::Union{Tuple{FT, FT}, Nothing} = nothing,
    context = ClimaComms.SingletonCommsContext(),
) where {FT}
    if isnothing(longlat)
        coord = ClimaCore.Geometry.ZPoint(z_sfc)
        space = (; surface = ClimaCore.Spaces.PointSpace(context, coord))
    else
        long, lat = longlat
        zlim = FT.((z_sfc - 1, z_sfc))
        nelements = 2
        box_domain = HybridBox(;
            xlim = FT.((long, long)),
            ylim = FT.((lat, lat)),
            zlim = zlim,
            longlat = longlat,
            nelements = (1, 1, nelements),
        )
        # Extract a point at the surface from the 3D space
        column_space =
            ClimaCore.Spaces.column(box_domain.space.subsurface_face, 1, 1, 1)
        point_space = ClimaCore.Spaces.level(
            column_space,
            nelements + ClimaCore.Utilities.half,
        )
        space = (; surface = point_space)
    end
    return Point{FT, typeof(space)}(z_sfc, space)
end


"""
    Column{FT} <: AbstractDomain{FT}
A struct holding the necessary information
to construct a domain, a mesh, a center and face
space, etc. for use when a finite difference in
1D is suitable, as for a soil column model.

`space` is a NamedTuple holding the surface space (in this case,
the top face space) and the center space for the subsurface.
These are stored using the keys :surface and :subsurface.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Column{FT, NT1 <: NamedTuple, NT2 <: NamedTuple} <: AbstractDomain{FT}
    "Domain interval limits, (zmin, zmax), in meters"
    zlim::Tuple{FT, FT}
    "Number of elements used to discretize the interval"
    nelements::Tuple{Int}
    "Tuple for mesh stretching specifying *target* (dz_bottom, dz_top) (m). If nothing, no stretching is applied."
    dz_tuple::Union{Tuple{FT, FT}, Nothing}
    "Boundary face identifiers"
    boundary_names::Tuple{Symbol, Symbol}
    "A NamedTuple of associated ClimaCore spaces: in this case, the surface space and subsurface center space"
    space::NT1
    "Fields and field data associated with the coordinates of the domain that are useful to store"
    fields::NT2
end

"""
    Column(;
        zlim::Tuple{FT, FT},
        nelements::Int,
        device = ClimaComms.device(),
        longlat::Union{Tuple{FT, FT}, Nothing} = nothing,
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing)

Outer constructor for the 1D `Column` domain.

Using `ClimaCore` tools, the coordinate mesh can be stretched such that the top
of the domain has finer resolution than the bottom of the domain. To activate this,
a tuple with the target dz_bottom and dz_top should be passed via keyword argument.
The default is uniform spacing. Please note that in order to use this feature, ClimaCore
requires that the elements of zlim be <=0. Additionally, the dz_tuple you supply may not
be compatible with the domain boundaries in some cases, in which case you may need to choose
different values.

The `boundary_names` field values are used to label the boundary faces
at the top and bottom of the domain.

If `longlat` is provided, the column is created at those coordinates.
This should be used for simulations on a column that require reading in spatial data
from a file defined on a lat/long grid, such as CLM data.
Note that constructing a Column with lat/long requires first constructing a 3D HybridBox
domain centered at `longlat`, then extracting a column from the intermediate 3D space.

The latitude and longitude of the returned domain can be accessed as follows:
- lat = get_lat(domain.space.surface)
- long = get_long(domain.space.surface)
"""
function Column(;
    zlim::Tuple{FT, FT},
    nelements::Int,
    device = ClimaComms.device(),
    longlat::Union{Tuple{FT, FT}, Nothing} = nothing,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
) where {FT}
    @assert zlim[1] < zlim[2]
    boundary_names = (:bottom, :top)

    if isnothing(longlat)
        column = ClimaCore.Domains.IntervalDomain(
            ClimaCore.Geometry.ZPoint{FT}(zlim[1]),
            ClimaCore.Geometry.ZPoint{FT}(zlim[2]);
            boundary_names = boundary_names,
        )
        if isnothing(dz_tuple)
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
        subsurface_space =
            ClimaCore.Spaces.CenterFiniteDifferenceSpace(device, mesh)
    else
        # Set limits at the longlat coordinates
        long, lat = longlat
        box_domain = HybridBox(;
            xlim = FT.((long, long)),
            ylim = FT.((lat, lat)),
            zlim = zlim,
            longlat = longlat,
            nelements = (1, 1, nelements),
            dz_tuple = dz_tuple,
        )
        # Extract a column from the 3D space
        colidx = ClimaCore.Grids.ColumnIndex((1, 1), 1)
        subsurface_space =
            ClimaCore.Spaces.column(box_domain.space.subsurface, colidx)
    end

    surface_space = obtain_surface_space(subsurface_space)
    subsurface_face_space = ClimaCore.Spaces.face_space(subsurface_space)
    space = (;
        surface = surface_space,
        subsurface = subsurface_space,
        subsurface_face = subsurface_face_space,
    )
    fields = get_additional_coordinate_field_data(subsurface_space)
    return Column{FT, typeof(space), typeof(fields)}(
        zlim,
        (nelements,),
        dz_tuple,
        boundary_names,
        space,
        fields,
    )
end


"""
    Plane{FT} <: AbstractDomain{FT}

A struct holding the necessary information
to construct a domain, a mesh, a 2d spectral
element space, and the resulting coordinate field.

When `longlat` is not nothing, the plane is assumed to be centered around these
coordinates. In this case, the curvature of the Earth is not accounted for.

`longlat` are in degrees, with longitude going from -180 to 180.

> :warning: Only independent columns are supported! (No lateral flow).

`space` is a NamedTuple holding the surface space (in this case,
the entire Plane space).
# Fields
$(DocStringExtensions.FIELDS)
"""
struct Plane{FT, NT <: NamedTuple} <: AbstractDomain{FT}
    "Domain interval limits along x axis, in meters or degrees (if `latlong != nothing`)"
    xlim::Tuple{FT, FT}
    "Domain interval limits along y axis, in meters or degrees (if `latlong != nothing`)"
    ylim::Tuple{FT, FT}
    "When not `nothing`, a Tuple that contains the center `long` and `lat` (in degrees, with `long`
    from -180 to 180)."
    longlat::Union{Nothing, Tuple{FT, FT}}
    "Number of elements to discretize interval, (nx, ny)"
    nelements::Tuple{Int, Int}
    "Flags for periodic boundaries. Only periodic or no lateral flow is supported."
    periodic::Tuple{Bool, Bool}
    "Polynomial order for both x and y"
    npolynomial::Int
    "A NamedTuple of associated ClimaCore spaces: in this case, the surface(Plane) space"
    space::NT
end

"""
    Plane(;
        xlim::Tuple{FT,FT},
        ylim::Tuple{FT,FT},
        nelements::Tuple{Int,Int},
        periodic::Tuple{Bool,Bool},
        npolynomial::Int,
        longlat = nothing,
        context = ClimaComms.context(),
        radius_earth = FT(6.378e6)
        ) where {FT}

Outer constructor for the `Plane` domain, using keyword arguments.

When not `nothing`, `longlat` is a Tuple that specifies the center of the plane on the Earth.
Example: `longlat = (FT(-118.14452), FT(34.14778))`. `longlat` is in degrees, with long from
-180 to 180. `longlat` is then converted to meters using `radius_earth`. In this, curvature
is not accounted for. We compute `xlim` and `ylim` as

```julia
xlim = (long - xlim[1] / (2π*radius_earth) * 360, long + xlim[2] / (2π*radius_earth) * 360)
```
"""
function Plane(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    nelements::Tuple{Int, Int},
    longlat = nothing,
    periodic::Tuple{Bool, Bool} = isnothing(longlat) ? (true, true) :
                                  (false, false),
    npolynomial::Int = 0,
    context = ClimaComms.context(),
    radius_earth = 6.378e6,
) where {FT}
    if isnothing(longlat)
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
    else
        @assert periodic == (false, false)
        radius_earth = FT(radius_earth)
        long, lat = longlat
        dxlim = abs.(xlim) # long bounds
        dylim = abs.(ylim) # lat bounds
        # Now make x refer to lat, and y refer to long,
        # for compatibility with ClimaCore
        xlim = (
            lat - dylim[1] / FT(2π * radius_earth) * 360,
            lat + dylim[2] / FT(2π * radius_earth) * 360,
        )
        ylim = (
            long - dxlim[1] / FT(2π * radius_earth) * 360,
            long + dxlim[2] / FT(2π * radius_earth) * 360,
        )

        @assert xlim[1] < xlim[2]
        @assert ylim[1] < ylim[2]

        # NOTE: We have LatLong instead of the other way because of ClimaCore
        domain_x = ClimaCore.Domains.IntervalDomain(
            ClimaCore.Geometry.LatPoint(xlim[1]),
            ClimaCore.Geometry.LatPoint(xlim[2]);
            boundary_names = (:north, :south),
        )
        domain_y = ClimaCore.Domains.IntervalDomain(
            ClimaCore.Geometry.LongPoint(ylim[1]),
            ClimaCore.Geometry.LongPoint(ylim[2]);
            boundary_names = (:west, :east),
        )
    end
    plane = ClimaCore.Domains.RectangleDomain(domain_x, domain_y)

    mesh = ClimaCore.Meshes.RectilinearMesh(plane, nelements[1], nelements[2])
    grid_topology = ClimaCore.Topologies.Topology2D(context, mesh)
    if npolynomial == 0
        quad = ClimaCore.Spaces.Quadratures.GL{npolynomial + 1}()
    else
        quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    end
    space = ClimaCore.Spaces.SpectralElementSpace2D(
        grid_topology,
        quad;
        compat_add_mask()...,
    )
    compat_set_mask!(space)
    space = (; surface = space)
    return Plane{FT, typeof(space)}(
        xlim,
        ylim,
        longlat,
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
        longlat::Union{Nothing, Tuple{FT, FT}},
        dz_tuple::Union{Tuple{FT, FT}, Nothing}
        nelements::Tuple{Int, Int, Int}
        npolynomial::Int
        periodic::Tuple{Bool, Bool}
    end

A struct holding the necessary information to construct a domain, a mesh, a 2d spectral
element space (horizontal) x a 1d finite difference space (vertical), and the resulting
coordinate field. This domain is not periodic along the z-axis. Note that no-flow boundary
conditions are supported in the horizontal.

When `longlat` is not `nothing`, assume that the box describes a region on the
globe centered around the `long` and `lat`.

`space` is a NamedTuple holding the surface space (in this case,
the top face space) and the center space for the subsurface.
These are stored using the keys :surface and :subsurface.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct HybridBox{FT, NT1 <: NamedTuple, NT2 <: NamedTuple} <: AbstractDomain{FT}
    "Domain interval limits along x axis, in meters or degrees (if `latlong != nothing`)"
    xlim::Tuple{FT, FT}
    "Domain interval limits along y axis, in meters or degrees (if `latlong != nothing`)"
    ylim::Tuple{FT, FT}
    "Domain interval limits along z axis, in meters"
    zlim::Tuple{FT, FT}
    "When not `nothing`, a Tuple that contains the center `long` and `lat`."
    longlat::Union{Nothing, Tuple{FT, FT}}
    "Tuple for mesh stretching specifying *target* (dz_bottom, dz_top) (m). If nothing, no stretching is applied."
    dz_tuple::Union{Tuple{FT, FT}, Nothing}
    "Number of elements to discretize interval, (nx, ny,nz)"
    nelements::Tuple{Int, Int, Int}
    " Polynomial order for the horizontal directions"
    npolynomial::Int
    "Flag indicating periodic boundaries in horizontal"
    periodic::Tuple{Bool, Bool}
    "A NamedTuple of associated ClimaCore spaces: in this case, the surface space and subsurface center space"
    space::NT1
    "Fields and field data associated with the coordinates of the domain that are useful to store"
    fields::NT2
end

"""
    HybridBox(;
    xlim::Tuple{FT, FT},
    ylim::Tuple{FT, FT},
    zlim::Tuple{FT, FT},
    nelements::Tuple{Int, Int, Int},
    npolynomial::Int = 0,
    device = ClimaComms.device(),
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    longlat = nothing,
    periodic::Tuple{Bool, Bool} = (true, true))

Constructs the `HybridBox` domain
 with limits `xlim` `ylim` and `zlim
(where `xlim[1] < xlim[2]`,`ylim[1] < ylim[2]`, and `zlim[1] < zlim[2]`),
`nelements` must be a tuple with three values, with the first
value corresponding
to the x-axis, the second corresponding to the y-axis, and the third
corresponding to the z-axis. The domain is periodic at the (xy) boundaries,
and the function space is of polynomial order `npolynomial` in the
horizontal directions.

When `longlat` is not `nothing`, assume that the box describes a region on the
globe centered around the `long` and `lat`. In this case, the radius of Earth is
assumed to be `6.378e6` meters and curvature is not included. We compute `xlim` and `ylim` as

```julia
xlim = (long - xlim[1] / (2π*radius_earth) * 360, long + xlim[2] / (2π*radius_earth) * 360)
```

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
    npolynomial::Int = 0,
    device = ClimaComms.device(),
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    longlat = nothing,
    periodic::Tuple{Bool, Bool} = isnothing(longlat) ? (true, true) :
                                  (false, false),
) where {FT}
    @assert zlim[1] < zlim[2]
    if isnothing(longlat)
        @assert xlim[1] < xlim[2]
        @assert ylim[1] < ylim[2]
        @assert periodic == (true, true)
    end
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(zlim[1]),
        ClimaCore.Geometry.ZPoint(zlim[2]);
        boundary_names = (:bottom, :top),
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
    vert_center_space =
        ClimaCore.Spaces.CenterFiniteDifferenceSpace(device, vertmesh)

    horzdomain = Plane(;
        xlim = xlim,
        ylim = ylim,
        longlat,
        nelements = nelements[1:2],
        periodic = periodic,
        npolynomial = npolynomial,
    )
    horzspace = horzdomain.space.surface

    subsurface_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )

    surface_space = obtain_surface_space(subsurface_space)
    subsurface_face_space = ClimaCore.Spaces.face_space(subsurface_space)
    space = (;
        surface = surface_space,
        subsurface = subsurface_space,
        subsurface_face = subsurface_face_space,
    )
    fields = get_additional_coordinate_field_data(subsurface_space)
    return HybridBox{FT, typeof(space), typeof(fields)}(
        horzdomain.xlim,
        horzdomain.ylim,
        zlim,
        longlat,
        dz_tuple,
        nelements,
        npolynomial,
        periodic,
        space,
        fields,
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

`space` is a NamedTuple holding the surface space (in this case,
the top face space) and the center space for the subsurface.
These are stored using the keys :surface and :subsurface.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalShell{FT, NT1 <: NamedTuple, NT2 <: NamedTuple} <:
       AbstractDomain{FT}
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
    "A NamedTuple of associated ClimaCore spaces: in this case, the surface space and subsurface center space"
    space::NT1
    "Fields and field data associated with the coordinates of the domain that are useful to store"
    fields::NT2
end

"""
    SphericalShell(;
        radius::FT,
        depth::FT,
        nelements::Tuple{Int, Int},
        npolynomial::Int,
        dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
        context = ClimaComms.context(),
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
    npolynomial::Int = 0,
    dz_tuple::Union{Tuple{FT, FT}, Nothing} = nothing,
    context = ClimaComms.context(),
) where {FT}
    @assert 0 < radius
    @assert 0 < depth
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(-depth)),
        ClimaCore.Geometry.ZPoint(FT(0));
        boundary_names = (:bottom, :top),
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
    device = ClimaComms.device(context)
    vert_center_space =
        ClimaCore.Spaces.CenterFiniteDifferenceSpace(device, vertmesh)

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements[1])
    horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
    if npolynomial == 0
        quad = ClimaCore.Spaces.Quadratures.GL{npolynomial + 1}()
    else
        quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    end
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(
        horztopology,
        quad;
        compat_add_mask()...,
    )
    compat_set_mask!(horzspace)
    subsurface_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
    surface_space = obtain_surface_space(subsurface_space)
    subsurface_face_space = ClimaCore.Spaces.face_space(subsurface_space)
    space = (;
        surface = surface_space,
        subsurface = subsurface_space,
        subsurface_face = subsurface_face_space,
    )
    fields = get_additional_coordinate_field_data(subsurface_space)
    return SphericalShell{FT, typeof(space), typeof(fields)}(
        radius,
        depth,
        dz_tuple,
        nelements,
        npolynomial,
        space,
        fields,
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

`space` is a NamedTuple holding the surface space (in this case,
the entire SphericalSurface space).
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalSurface{FT, NT <: NamedTuple} <: AbstractDomain{FT}
    "The radius of the surface"
    radius::FT
    "The number of elements to be used in the non-radial directions"
    nelements::Int
    "The polynomial order to be used in the non-radial directions"
    npolynomial::Int
    "A NamedTuple of associated ClimaCore spaces: in this case, the surface (SphericalSurface) space"
    space::NT
end

"""
    SphericalSurface(;
        radius::FT,
        nelements::Int
        npolynomial::Int,
        context = ClimaComms.context(),
    ) where {FT}
Outer constructor for the `SphericalSurface` domain, using keyword arguments.
"""
function SphericalSurface(;
    radius::FT,
    nelements::Int,
    npolynomial::Int = 0,
    context = ClimaComms.context(),
) where {FT}
    @assert 0 < radius
    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = Meshes.EquiangularCubedSphere(horzdomain, nelements)
    horztopology = Topologies.Topology2D(context, horzmesh)
    if npolynomial == 0
        quad = ClimaCore.Spaces.Quadratures.GL{npolynomial + 1}()
    else
        quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    end
    horzspace =
        Spaces.SpectralElementSpace2D(horztopology, quad; compat_add_mask()...)
    compat_set_mask!(horzspace)
    space = (; surface = horzspace)
    return SphericalSurface{FT, typeof(space)}(
        radius,
        nelements,
        npolynomial,
        space,
    )
end


"""
    obtain_surface_domain(d::AbstractDomain) where {FT}

Default method throwing an error; any domain with a corresponding
domain should define a new method of this function.
"""
obtain_surface_domain(d::AbstractDomain) =
    @error("No surface domain is defined for this domain.")

"""
    obtain_surface_domain(c::Column{FT}) where {FT}

Returns the Point domain corresponding to the top face (surface) of the
Column domain `c`.
"""
function obtain_surface_domain(c::Column{FT}) where {FT}
    surface_domain = Point{FT, typeof((; surface = c.space.surface))}(
        c.zlim[2],
        (; surface = c.space.surface),
    )
    return surface_domain
end

"""
    obtain_surface_domain(b::HybridBox{FT}) where {FT}

Returns the Plane domain corresponding to the top face (surface) of the
HybridBox domain `b`.
"""
function obtain_surface_domain(b::HybridBox{FT}) where {FT}
    surface_domain = Plane{FT, typeof((; surface = b.space.surface))}(
        b.xlim,
        b.ylim,
        b.longlat,
        b.nelements[1:2],
        b.periodic,
        b.npolynomial,
        (; surface = b.space.surface),
    )
    return surface_domain
end


"""
    obtain_surface_domain(s::SphericalShell{FT}) where {FT}

Returns the SphericalSurface domain corresponding to the top face
(surface) of the SphericalShell domain `s`.
"""
function obtain_surface_domain(s::SphericalShell{FT}) where {FT}
    surface_domain =
        SphericalSurface{FT, typeof((; surface = s.space.surface))}(
            s.radius,
            s.nelements[1],
            s.npolynomial,
            (; surface = s.space.surface),
        )
    return surface_domain
end

"""
    obtain_surface_space(cs::ClimaCore.Spaces.AbstractSpace)

Returns the surface space, if applicable, for the center space `cs`.
"""
obtain_surface_space(cs::ClimaCore.Spaces.AbstractSpace) =
    @error("No surface space is defined for this space.")


"""
    obtain_surface_space(cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace)

Returns the horizontal space for the CenterExtrudedFiniteDifferenceSpace `cs`.
"""
function obtain_surface_space(
    cs::ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
)
    return Spaces.horizontal_space(cs)
end

"""
    obtain_surface_space(cs::ClimaCore.Spaces.FiniteDifferenceSpace)

Returns the top level of the face space corresponding to the input space `cs`.
"""
function obtain_surface_space(cs::ClimaCore.Spaces.FiniteDifferenceSpace)
    fs = ClimaCore.Spaces.face_space(cs)
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
    center_space = axes(center_field)
    N = ClimaCore.Spaces.nlevels(center_space)
    surface_space = obtain_surface_space(center_space)
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(
            ClimaCore.Fields.level(center_field, N),
        ),
        surface_space,
    )
end

"""
    top_center_to_surface(val)

When `val` is a scalar (e.g. a single float or struct), returns `val`.
"""
top_center_to_surface(val) = val

"""
    linear_interpolation_to_surface!(sfc_field, center_field, z, Δz_top)

Linearly interpolate the center field `center_field` to the surface
defined by the top face coordinate of `z` with a center to face distance
`Δz_top` in the first layer; updates the `sfc_field`
on the surface (face) space in place.
"""
function linear_interpolation_to_surface!(sfc_field, center_field, z, Δz_top)
    surface_space = axes(sfc_field)
    Δz_top = ClimaCore.Fields.field_values(Δz_top)
    nz = Spaces.nlevels(axes(center_field))
    f1 = ClimaCore.Fields.field_values(ClimaCore.Fields.level(center_field, nz))
    f2 = ClimaCore.Fields.field_values(
        ClimaCore.Fields.level(center_field, nz - 1),
    )
    z1 = ClimaCore.Fields.field_values(ClimaCore.Fields.level(z, nz))
    z2 = ClimaCore.Fields.field_values(ClimaCore.Fields.level(z, nz - 1))
    ClimaCore.Fields.field_values(sfc_field) .=
        @. (f1 - f2) / (z1 - z2) * (Δz_top + z1 - z2) + f2
end

"""
    bottom_center_to_surface(center_field::ClimaCore.Fields.Field)

Creates and returns a ClimaCore.Fields.Field defined on the space
corresponding to the bottom surface of the space on which `center_field`
is defined, with values equal to the those at the level of the bottom
center.

For example, given a `center_field` defined on 1D center finite difference space,
this would return a field defined on the Point space of the bottom surface of
the column. The value would be the value of the original `center_field`
at the bottommost location. Given a `center_field` defined on a 3D
extruded center finite difference space, this would return a 2D field
corresponding to the bottom surface, with values equal to the bottommost level.
"""
function bottom_center_to_surface(center_field::ClimaCore.Fields.Field)
    center_space = axes(center_field)
    surface_space = obtain_surface_space(center_space)
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(ClimaCore.Fields.level(center_field, 1)),
        surface_space,
    )
end

"""
    bottom_center_to_surface(val)

When `val` is a scalar (e.g. a single float or struct), returns `val`.
"""
bottom_center_to_surface(val) = val

"""
    get_Δz(z::ClimaCore.Fields.Field)

A function to return a tuple containing the distance between the top boundary
and its closest center, and the bottom boundary and its closest center,
both as Fields. It also returns the widths of each layer as a field.
"""
function get_Δz(z::ClimaCore.Fields.Field)
    # Extract the differences between levels of the face space
    fs = ClimaCore.Spaces.face_space(axes(z))
    z_face = ClimaCore.Fields.coordinate_field(fs).z
    Δz_face = ClimaCore.Fields.Δz_field(z_face)
    Δz_top = ClimaCore.Fields.level(
        Δz_face,
        ClimaCore.Utilities.PlusHalf(ClimaCore.Spaces.nlevels(fs) - 1),
    )
    Δz_bottom = ClimaCore.Fields.level(Δz_face, ClimaCore.Utilities.PlusHalf(0))

    #Layer widths:
    Δz_center = ClimaCore.Fields.Δz_field(z)
    return Δz_top ./ 2, Δz_bottom ./ 2, Δz_center
end

"""
    get_lat(surface_space::ClimaCore.Spaces.PointSpace)
    get_lat(subsurface_space::ClimaCore.Spaces.FiniteDifferenceSpace)

Returns the latitude of provided surface space as a Field defined on the space.
If the space does not have latitude information, an error is thrown.

This function is not implemented for subsurface spaces, since we want latitude
as a 2D Field, not 3D.
"""
function get_lat(
    surface_space::Union{
        ClimaCore.Spaces.PointSpace,
        ClimaCore.Spaces.SpectralElementSpace2D,
    },
)
    if hasproperty(ClimaCore.Fields.coordinate_field(surface_space), :lat)
        return ClimaCore.Fields.coordinate_field(surface_space).lat
    else
        error("Surface space does not have latitude information")
    end
end
get_lat(
    _::Union{
        ClimaCore.Spaces.FiniteDifferenceSpace,
        ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
    },
) = error("`get_lat` is not implemented for subsurface spaces")

"""
    get_long(surface_space::ClimaCore.Spaces.PointSpace)
    get_long(subsurface_space::ClimaCore.Spaces.FiniteDifferenceSpace)

Returns the longitude of the provided surface space as a Field defined on the space.
If the space does not have longitude information, an error is thrown.

This function is not implemented for subsurface spaces, since we want longitude
as a 2D Field, not 3D.
"""
function get_long(
    surface_space::Union{
        ClimaCore.Spaces.PointSpace,
        ClimaCore.Spaces.SpectralElementSpace2D,
    },
)
    if hasproperty(ClimaCore.Fields.coordinate_field(surface_space), :long)
        return ClimaCore.Fields.coordinate_field(surface_space).long
    else
        error("Surface space does not have longitude information")
    end
end
get_long(
    _::Union{
        ClimaCore.Spaces.FiniteDifferenceSpace,
        ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
    },
) = error("`get_long` is not implemented for subsurface spaces")

"""
    top_face_to_surface(face_field::ClimaCore.Fields.Field, surface_space)

Creates and returns a ClimaCore.Fields.Field defined on the space
corresponding to the surface of the space on which `face_field`
is defined, with values equal to the those at the level of the top
face.

Given a `face_field` defined on a 3D
extruded face finite difference space, this would return a 2D field
corresponding to the surface, with values equal to the topmost level.
 """
function top_face_to_surface(face_field::ClimaCore.Fields.Field, surface_space)
    face_space = axes(face_field)
    N = ClimaCore.Spaces.nlevels(face_space)
    sfc_level =
        ClimaCore.Fields.level(face_field, ClimaCore.Utilities.PlusHalf(N - 1))
    # Project onto surface space
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(sfc_level),
        surface_space,
    )
end

"""
    get_additional_coordinate_field_data(subsurface_space)

A helper function which returns additional fields and field data corresponding to ClimaLand
domains which have a subsurface_space (Column, HybridBox, SphericalShell).
The fields are the center coordinates of the subsurface space, the spacing between
the top center and top surface and bottom center and bottom surface, as well as the
field corresponding to the surface height z and layer widths. The field data are the
depth of the domain and the minimum top layer thickness over the entire domain.

We allocate these once, upon domain construction, so that they are accessible
during the simulation.
"""
function get_additional_coordinate_field_data(subsurface_space)
    surface_space = obtain_surface_space(subsurface_space)
    z = ClimaCore.Fields.coordinate_field(subsurface_space).z
    Δz_top, Δz_bottom, Δz = get_Δz(z)
    face_space = ClimaCore.Spaces.face_space(subsurface_space)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    z_sfc = top_face_to_surface(z_face, surface_space)
    d = depth(subsurface_space)
    Δz_min = minimum(Δz)
    fields = (;
        z = z,
        Δz_top = Δz_top,
        Δz_bottom = Δz_bottom,
        z_sfc = z_sfc,
        depth = d,
        Δz = Δz,
        Δz_min = Δz_min,
    )
    return fields
end

"""
    depth(space::Union{ClimaCore.Spaces.CenterFiniteDifferenceSpace,
                       ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace})

Returns the depth of the domain as a scalar. Note that these functions
will need to be modified upon the introduction of
- topography at surface
- depth to bedrock (topography at bottom)

Since the depth will be a field in this case, it should be allocated and
stored in domain.fields, which is why we store it there now even though it is not a field.
"""
function depth(
    space::Union{
        ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.CenterFiniteDifferenceSpace,
    },
)
    zmin, zmax = extrema(
        ClimaCore.Fields.coordinate_field(ClimaCore.Spaces.face_space(space)).z,
    )
    return zmax - zmin
end

"""
    horizontal_resolution_degrees(domain::AbstractDomain)

Return a tuple with the approximate resolution on the domain in degrees along
the two directions.

For boxes and planes, the order is `(latitude, longitude)`.

Examples
=======

```jldoctest
julia> using ClimaLand

julia> domain = ClimaLand.Domains.SphericalShell(;
        radius = 6300e3,
        depth = 15.,
        nelements = (10, 3),
        dz_tuple = ((1.0, 0.05)));

julia> ClimaLand.Domains.average_horizontal_resolution_degrees(domain)
(9.0, 9.0)

julia> domain = ClimaLand.Domains.Plane(;
        xlim = (50000.0, 80000.),
        ylim = (30000.0, 40000.),
        longlat = (-118.14452, 34.14778),
        nelements = (20, 3));

julia> round.(ClimaLand.Domains.average_horizontal_resolution_degrees(domain), digits = 3)
(0.031, 0.389)
```
"""
function average_horizontal_resolution_degrees(
    domain::Union{SphericalShell, SphericalSurface},
)
    num_elements_lat = num_elements_lon = first(domain.nelements)
    quad = ClimaCore.Spaces.quadrature_style(domain.space.surface)
    num_points_per_element =
        ClimaCore.Quadratures.unique_degrees_of_freedom(quad)

    # There are 2 full cubed-sphere panels along latitudes, and 4 along
    # longitudes
    return (
        180 / (2 * num_points_per_element * num_elements_lat),
        360 / (4 * num_points_per_element * num_elements_lon),
    )
end

function average_horizontal_resolution_degrees(domain::Union{Plane, HybridBox})
    isnothing(domain.longlat) &&
        error("Can only compute resolution with domains in latlong")

    # x is Lat, y is Lon because ClimaCore defines LatLongPoints
    num_elements_lat, num_elements_lon = domain.nelements[1:2]
    delta_lat = domain.xlim[2] - domain.xlim[1]
    delta_lon = domain.ylim[2] - domain.ylim[1]

    quad = ClimaCore.Spaces.quadrature_style(domain.space.surface)
    num_points_per_element =
        ClimaCore.Quadratures.unique_degrees_of_freedom(quad)

    return (
        delta_lat / (num_points_per_element * num_elements_lat),
        delta_lon / (num_points_per_element * num_elements_lon),
    )
end

apply_threshold(field, value) =
    field > value ? eltype(field)(1) : eltype(field)(0)

"""
    landsea_mask(
        surface_space;
        filepath = landseamask_file_path(;
                 context = ClimaComms.context(surface_space),
                 ),
        threshold = 0.5,
        regridder_type = :InterpolationsRegridder,
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
       interpolation_method = Interpolations.Constant()
    )

Reads in the default Clima land/sea mask, regrids to the `surface_space`, and
treats any point with a land fraction < threshold as ocean.

Note that by default we use nearest-neighbor interpolation, in which case
the threshold does nothing since the land-sea mask is natively 0/1, and only
becomes a non-integer number when linearly interpolating. You can change to linear
interpolation by passing `interpolation_method = Interpolations.Linear()`.
"""
function landsea_mask(
    surface_space;
    filepath = landseamask_file_path(;
        context = ClimaComms.context(surface_space),
    ),
    threshold = 0.5,
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant(),
)
    mask = SpaceVaryingInput(
        filepath,
        "landsea",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    binary_mask = apply_threshold.(mask, threshold)
    return binary_mask
end

# Points and Columns do not have a horizontal dim, so a horizontal mask cannot be applied
landsea_mask(domain::Union{Point, Column}; kwargs...) = nothing

function landsea_mask(
    domain::Union{SphericalShell, SphericalSurface, HybridBox, Plane};
    kwargs...,
)
    # HybridBox and Plane domains might not have longlat, which is needed for the mask
    (hasproperty(domain, :longlat) && isnothing(domain.longlat)) &&
        return nothing
    # average_horizontal_resolution_degrees returns a tuple with the resolution
    # along the two directions, so we take the minimum
    resolution_degrees = minimum(average_horizontal_resolution_degrees(domain))
    if resolution_degrees > 1
        resolution = "1deg"
    else
        resolution_arcsec = 3600resolution_degrees
        # Pick the landsea mask at 60 arcseconds if the nodal distance is more than
        # 120", otherwise pick the higher resolution
        resolution = resolution_arcsec > 120 ? "60arcs" : "30arcs"
    end

    filepath = landseamask_file_path(;
        resolution,
        context = ClimaComms.context(domain.space.surface),
    )
    return landsea_mask(domain.space.surface; filepath, kwargs...)
end


"""
    global_domain(
    FT;
    apply_mask = true,
    mask_threshold = 0.5,
    nelements = (101, 15),
    dz_tuple = (10.0, 0.05),
    depth = 50.0,
    npolynomial = 0,
    context = ClimaComms.context(),
    filepath = landseamask_file_path(;context),
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant()
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
    filepath = landseamask_file_path(; context),
    regridder_type = :InterpolationsRegridder,
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    ),
    interpolation_method = Interpolations.Constant(),
)
    if pkgversion(ClimaCore) < v"0.14.30" && apply_mask
        @warn "The land mask cannot be applied with ClimaCore < v0.14.30. Update ClimaCore for significant performance gains."
        apply_mask = false
    end

    radius = FT(6378.1e3)
    dz_tuple = FT.(dz_tuple)
    depth = FT(depth)
    domain = SphericalShell(;
        radius,
        depth,
        nelements,
        npolynomial,
        dz_tuple,
        context,
    )
    if apply_mask
        surface_space = domain.space.surface # 2d space
        binary_mask = landsea_mask(
            domain;
            threshold = mask_threshold,
            regridder_type,
            extrapolation_bc,
            interpolation_method,
            filepath,
        )
        Spaces.set_mask!(surface_space, binary_mask)
    end

    return domain
end


"""
    use_lowres_clm(space)

Returns true if the simulation space is closer in resolution to the low
resolution  CLM data (at 0.9x1.25 degree lat/long) compared with the
high resolution CLM data ( 0.125x0.125 degree lat/long).

If the simulation space is closer to the low resolution CLM data,
that data will be used, and vice versa. The high resolution data
takes longer to process and so using lower resolution data when it
suffices is desirable.

If the surface space is a point, `use_lowres_clm` always returns true.
"""
function use_lowres_clm(
    surface_space::ClimaCore.Spaces.AbstractSpectralElementSpace,
)
    node_scale = ClimaCore.Spaces.node_horizontal_length_scale(surface_space)
    surface_mesh = ClimaCore.Spaces.topology(surface_space).mesh
    if surface_mesh isa ClimaCore.Meshes.AbstractCubedSphere
        # in this case, node_scale is in meters
        sphere_radius = surface_mesh.domain.radius
        horizontal_length_scale(lat_res, long_res) = sqrt(
            4 * pi * sphere_radius^2 / ((360 / long_res) * (180 / lat_res)),
        )
        highres_scale = horizontal_length_scale(0.125, 0.125)
        lowres_scale = horizontal_length_scale(0.9, 1.25)
    elseif surface_mesh isa ClimaCore.Meshes.RectilinearMesh
        # in this case, node_scale is in degrees
        highres_scale = 0.125
        lowres_scale = sqrt(1.25 * 0.9)
    else
        return false
    end
    return abs(lowres_scale - node_scale) < abs(highres_scale - node_scale)
end

use_lowres_clm(surface_space::ClimaCore.Spaces.PointSpace) = true


export AbstractDomain
export Column, Plane, HybridBox, Point, SphericalShell, SphericalSurface
export coordinates,
    obtain_surface_space,
    top_center_to_surface,
    bottom_center_to_surface,
    top_face_to_surface,
    obtain_surface_domain,
    linear_interpolation_to_surface!,
    get_Δz,
    average_horizontal_resolution_degrees,
    use_lowres_clm,
    global_domain
end
