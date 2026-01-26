"""
    AbstractZSamplingMethod

The `AbstractZInterpolationMethod` defines how points along the vertical axis should be
sampled.

In other words, if a column is defined between 0 and 100 and the target number of points is
50, how should those 50 points be chosen?

Available methods are:
- `LevelMethod`: just use the grid levels
- `FakePressureLevelsMethod`: linearly spaced in (very) approximate atmospheric pressure levels
"""
abstract type AbstractZSamplingMethod end

"""
    LevelsMethod

Do not perform interpolation on `z`, use directly the grid levels instead.
"""
struct LevelsMethod <: AbstractZSamplingMethod end

"""
    FakePressureLevelsMethod

Linearly sample points from `z_min` to `z_max` in pressure levels assuming a very simplified hydrostatic balance model.

Pressure is approximated with

p ~ p₀ exp(-z/H)

H is assumed to be 7000 m, which is a good scale height for the Earth atmosphere.
"""
struct FakePressureLevelsMethod <: AbstractZSamplingMethod end

"""
    add_dimension!(nc::NCDatasets.NCDataset,
                   name::String,
                   points;
                   kwargs...)

Add dimension identified by `name` in the given `nc` file and fill it with the given
`points`.
"""
function add_dimension!(
    nc::NCDatasets.NCDataset,
    name::String,
    points;
    kwargs...,
)
    FT = eltype(points)

    NCDatasets.defDim(nc, name, size(points)[end])

    dim = NCDatasets.defVar(nc, name, FT, (name,))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end

    dim[:] = points

    return nothing
end

"""
    dimension_exists(
        nc::NCDatasets.NCDataset,
        name::String,
        expected_size::Tuple,
        )

Return whether the given dimension exists in the given dataset, and if yes, it has the same
size as `expected_size`.
"""
function dimension_exists(
    nc::NCDatasets.NCDataset,
    name::String,
    expected_size::Tuple,
)
    if haskey(nc, name)
        if size(nc[name]) != expected_size
            file_path = NCDatasets.path(nc)
            error(
                "Incompatible $name dimension already exists in file $file_path",
            )
        else
            return true
        end
    else
        return false
    end
end

"""
    add_time_maybe!(nc::NCDatasets.NCDataset,
                    float_type::Type{FT};
                    kwargs...) where {FT}

Add the `time` dimension (with infinite size) to the given NetCDF file if not already there.
Optionally, add all the keyword arguments as attributes.
"""
function add_time_maybe!(
    nc::NCDatasets.NCDataset,
    float_type::Type{FT};
    kwargs...,
) where {FT}

    # If we already have time, do nothing
    haskey(nc, "time") && return nothing

    NCDatasets.defDim(nc, "time", Inf)
    dim = NCDatasets.defVar(nc, "time", FT, ("time",))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end
    return nothing
end

"""
    add_date_maybe!(nc::NCDatasets.NCDataset,
                    kwargs...)

Add the `date` dimension (with infinite size) to the given NetCDF file if not already there.
Optionally, add all the keyword arguments as attributes.
"""
function add_date_maybe!(nc::NCDatasets.NCDataset; kwargs...)
    haskey(nc, "date") && return nothing
    !haskey(nc, "time") &&
        error("`add_date_maybe!` is meant to be called after `add_time_maybe!`")
    dim = NCDatasets.defVar(nc, "date", Float64, ("time",))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end
end

"""
    add_time_bounds_maybe!(nc::NCDatasets.NCDataset,
                           float_type::Type{FT};
                           kwargs...)

Add the `time_bnds` dimension (with infinite size) to the given NetCDF file if
not already there. Optionally, add all the keyword arguments as attributes.

This function is meant to be called after `add_time_maybe!`.
"""
function add_time_bounds_maybe!(
    nc::NCDatasets.NCDataset,
    float_type::Type{FT};
    kwargs...,
) where {FT}
    haskey(nc, "time_bnds") && return nothing
    !haskey(nc, "time") && error(
        "Time dimension does not exist. Call add_time_maybe! before calling this function",
    )
    !("nv" in NCDatasets.dimnames(nc)) && NCDatasets.defDim(nc, "nv", 2) # number of vertices
    dim = NCDatasets.defVar(nc, "time_bnds", FT, ("nv", "time"))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end
    return nothing
end

"""
    add_date_bounds_maybe!(nc::NCDatasets.NCDataset; kwargs...)

Add the `date_bnds` dimension (with infinite size) to the given NetCDF file if
not already there. Optionally, add all the keyword arguments as attributes.

This function is meant to be called after `add_time_maybe!`.
"""
function add_date_bounds_maybe!(nc::NCDatasets.NCDataset; kwargs...)
    haskey(nc, "date_bnds") && return nothing
    !haskey(nc, "date") && error(
        "Time dimension does not exist. Call add_time_maybe! before calling this function",
    )
    !("nv" in NCDatasets.dimnames(nc)) && NCDatasets.defDim(nc, "nv", 2) # number of vertices
    dim = NCDatasets.defVar(nc, "date_bnds", Float64, ("nv", "time"))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end
    return nothing
end

"""
    add_space_coordinates_maybe!(nc::NCDatasets.NCDataset,
                                 space::Spaces.AbstractSpace,
                                 num_points;
                                 names)

Add dimensions relevant to the `space` to the given `nc` NetCDF file. The range is
automatically determined and the number of points is set with `num_points`, which has to be
an iterable of size N, where N is the number of dimensions of the space. For instance, 3 for
a cubed sphere, 2 for a surface, 1 for a column.

The function returns an array with the names of the relevant dimensions. (We want arrays
because we want to preserve the order to match the one in num_points).

In some cases, the names are adjustable passing the keyword `names`.
"""
function add_space_coordinates_maybe! end

"""
    target_coordinates(space, num_points)

Return the range of interpolation coordinates. The range is automatically determined and the
number of points is set with `num_points`, which has to be an iterable of size N, where N is
the number of dimensions of the space. For instance, 3 for a cubed sphere, 2 for a surface,
1 for a column.
"""
function target_coordinates(space, num_points) end

function target_coordinates(
    space::S,
    num_points,
    z_sampling_method::FakePressureLevelsMethod,
) where {
    S <:
    Union{Spaces.CenterFiniteDifferenceSpace, Spaces.FaceFiniteDifferenceSpace},
}
    # Exponentially spaced with base e
    #
    # We mimic something that looks like pressure levels
    #
    # p ~ p₀ exp(-z/H)
    #
    # We assume H to be 7000, which is a good scale height for the Earth atmosphere
    H_EARTH = 7000

    num_points_z = last(num_points)
    FT = Spaces.undertype(space)
    topology = Spaces.topology(space)
    vert_domain = topology.mesh.domain
    z_min, z_max = FT(vert_domain.coord_min.z), FT(vert_domain.coord_max.z)
    # We floor z_min to avoid having to deal with the singular value z = 0.
    z_min = max(z_min, 100)
    exp_z_min = exp(-z_min / H_EARTH)
    exp_z_max = exp(-z_max / H_EARTH)
    return collect(-H_EARTH * log.(range(exp_z_min, exp_z_max, num_points_z)))
end

function target_coordinates(
    space::S,
    num_points,
    z_sampling_method::LevelsMethod,
) where {
    S <:
    Union{Spaces.CenterFiniteDifferenceSpace, Spaces.FaceFiniteDifferenceSpace},
}
    cspace = Spaces.space(space, Grids.CellCenter())
    return Array(parent(Fields.coordinate_field(cspace).z))[:, 1]
end

# Column
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.FiniteDifferenceSpace,
    num_points_z;
    z_sampling_method,
    names = ("z",),
    interpolated_physical_z = nothing, # Not needed here, but needed for consistency of
    # interface and dispatch
)
    name, _... = names
    z_dimension_exists = dimension_exists(nc, name, num_points_z)

    if !z_dimension_exists
        zpts = target_coordinates(space, num_points_z, z_sampling_method)
        add_dimension!(nc, name, zpts, units = "m", axis = "Z")
    end
    return [name]
end

# PointSpace
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.PointSpace,
    num_points_z;
    z_sampling_method,
    names = (),
    interpolated_physical_z = nothing, # Not needed here, but needed for consistency of
    # interface and dispatch
)
    return []
end

add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.AbstractSpectralElementSpace,
    num_points;
) = add_space_coordinates_maybe!(
    nc,
    space,
    num_points,
    Meshes.domain(Spaces.topology(space));
)


# For the horizontal space, we also have to look at the domain, so we define another set of
# functions that dispatches over the domain
target_coordinates(space::Spaces.AbstractSpectralElementSpace, num_points) =
    target_coordinates(space, num_points, Meshes.domain(Spaces.topology(space)))

# Box
function target_coordinates(
    space::Spaces.SpectralElementSpace2D,
    num_points,
    domain::Domains.RectangleDomain,
)
    if islatlonbox(domain)
        # ClimaCore assumes LatLon, but we really want LongLat, so we need to flip
        # This function return Lat - Long points
        num_points_y, num_points_x = num_points
    else
        num_points_x, num_points_y = num_points
    end

    xmin = Geometry.tofloat(domain.interval1.coord_min)
    xmax = Geometry.tofloat(domain.interval1.coord_max)
    ymin = Geometry.tofloat(domain.interval2.coord_min)
    ymax = Geometry.tofloat(domain.interval2.coord_max)
    # Case of box with one single point
    if num_points_x == 1
        xpts = [(xmax + xmin) / 2]
    else
        xpts = collect(range(xmin, xmax, num_points_x))
    end
    if num_points_y == 1
        ypts = [(ymax + ymin) / 2]
    else
        ypts = collect(range(ymin, ymax, num_points_y))
    end
    return (xpts, ypts)
end

# Plane
function target_coordinates(
    space::Spaces.SpectralElementSpace1D,
    num_points,
    domain::Domains.IntervalDomain,
)
    num_points_x, _... = num_points
    FT = Spaces.undertype(space)
    xmin = Geometry.tofloat(domain.coord_min)
    xmax = Geometry.tofloat(domain.coord_max)
    xpts = collect(range(xmin, xmax, num_points_x))
    return (xpts)
end

# Cubed sphere
function target_coordinates(
    space::Spaces.SpectralElementSpace2D,
    num_points,
    ::Domains.SphereDomain,
)
    num_points_long, num_points_lat = num_points
    FT = Spaces.undertype(space)
    longpts = collect(range(FT(-180), FT(180), num_points_long))
    latpts = collect(range(FT(-90), FT(90), num_points_lat))

    return (longpts, latpts)
end

islatlonbox(space::Spaces.PointSpace) = false
islatlonbox(space::Spaces.FiniteDifferenceSpace) = false
islatlonbox(space::Domains.AbstractDomain) = false
function islatlonbox(space::Spaces.AbstractSpace)
    return islatlonbox(
        Meshes.domain(Spaces.topology(Spaces.horizontal_space(space))),
    )
end

# Box
function islatlonbox(domain::Domains.RectangleDomain)
    return domain.interval1.coord_max isa Geometry.LatPoint &&
           domain.interval2.coord_max isa Geometry.LongPoint
end


function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace2D,
    num_points,
    domain::Domains.RectangleDomain;
    names = islatlonbox(domain) ? ("lon", "lat") : ("x", "y"),
)
    name1, name2 = names
    num_points1, num_points2 = num_points

    maybe_reverse = identity
    if islatlonbox(domain)
        more_attribs = (
            Dict(
                :units => "degrees_east",
                :axis => "X",
                :standard_name => "longitude",
                :long_name => "Longitude",
            ),
            Dict(
                :units => "degrees_north",
                :axis => "Y",
                :standard_name => "latitude",
                :long_name => "Latitude",
            ),
        )
        # ClimaCore assumes LatLon, but we really want LongLat, so we need to flip
        maybe_reverse = reverse
    else
        more_attribs = (
            Dict(:units => "m", :axis => "X"),
            Dict(:units => "m", :axis => "Y"),
        )
    end

    dim1_exists = dimension_exists(nc, name1, (num_points1,))
    dim2_exists = dimension_exists(nc, name2, (num_points2,))

    if !dim1_exists && !dim2_exists
        pts1, pts2 = maybe_reverse(target_coordinates(space, num_points))
        add_dimension!(nc, name1, pts1; more_attribs[1]...)
        add_dimension!(nc, name2, pts2; more_attribs[2]...)
    end

    return [name1, name2]
end

# Plane
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace1D,
    num_points,
    ::Domains.IntervalDomain;
    names = ("x",),
)
    xname, _... = names
    num_points_x, = num_points
    x_dimension_exists = dimension_exists(nc, xname, (num_points_x,))

    if !x_dimension_exists
        xpts = target_coordinates(space, num_points)
        add_dimension!(nc, xname, xpts; units = "m", axis = "X")
    end
    return [xname]
end

# Cubed sphere
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace2D,
    num_points,
    ::Domains.SphereDomain;
    names = ("lon", "lat"),
)
    longname, latname = names
    num_points_long, num_points_lat = num_points

    long_dimension_exists = dimension_exists(nc, longname, (num_points_long,))
    lat_dimension_exists = dimension_exists(nc, latname, (num_points_lat,))

    if !long_dimension_exists && !lat_dimension_exists
        longpts, latpts = target_coordinates(space, num_points)
        add_dimension!(
            nc,
            longname,
            longpts;
            units = "degrees_east",
            axis = "X",
            standard_name = "longitude",
            long_name = "Longitude",
        )
        add_dimension!(
            nc,
            latname,
            latpts;
            units = "degrees_north",
            axis = "Y",
            standard_name = "latitude",
            long_name = "Latitude",
        )
    end

    return [longname, latname]
end

# General hybrid space. This calls both the vertical and horizontal add_space_coordinates_maybe!
# and combines the resulting dictionaries
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.ExtrudedFiniteDifferenceSpace,
    num_points;
    z_sampling_method,
    interpolated_physical_z = nothing,
)

    hdims_names = vdims_names = []

    num_points_horiz..., num_points_vertic = num_points

    # Being an Extruded space, we can assume that we have an horizontal and a vertical space.
    # We can also assume that the vertical space has dimension 1
    horizontal_space = Spaces.horizontal_space(space)

    hdims_names =
        add_space_coordinates_maybe!(nc, horizontal_space, num_points_horiz)

    vertical_space = Spaces.FiniteDifferenceSpace(
        Spaces.grid(space).vertical_grid,
        Spaces.staggering(space),
    )

    if Spaces.grid(space).hypsography isa Grids.Flat
        vdims_names = add_space_coordinates_maybe!(
            nc,
            vertical_space,
            (num_points_vertic,);
            z_sampling_method,
        )
    else
        vdims_names = add_space_coordinates_maybe!(
            nc,
            vertical_space,
            (num_points_vertic,),
            interpolated_physical_z;
            z_sampling_method,
            names = ("z_reference",),
            depending_on_dimensions = hdims_names,
        )
    end

    return (hdims_names..., vdims_names...)
end

# Ignore the interpolated_physical_z/z_sampling_method keywords in the general
# case (we only case about the specialized one for extruded spaces)
add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space,
    num_points;
    interpolated_physical_z = nothing,
    z_sampling_method = nothing,
) = add_space_coordinates_maybe!(nc::NCDatasets.NCDataset, space, num_points)

# Elevation with topography

# `depending_on_dimensions` identifies the dimensions upon which the current one depends on
# (excluding itself). In pretty much all cases, the dimensions depend only on themselves
# (e.g., `lat` is a variable only defined on the latitudes.), and `depending_on_dimensions`
# should be an empty tuple. The only case in which this is not what happens is with `z` with
# topography. With topography, the altitude will depend on the spatial coordinates. So,
# `depending_on_dimensions` might be `("lon", "lat)`, or similar.
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.FiniteDifferenceSpace,
    num_points,
    interpolated_physical_z;
    names = ("z_reference",),
    z_sampling_method,
    depending_on_dimensions,
)
    name, _... = names

    # Add z_reference
    z_reference_dimension_dimension_exists =
        dimension_exists(nc, name, num_points)

    if !z_reference_dimension_dimension_exists
        reference_altitudes =
            target_coordinates(space, num_points, z_sampling_method)
        add_dimension!(nc, name, reference_altitudes; units = "m", axis = "Z")
    end

    # We also have to add an extra variable with the physical altitudes
    physical_name = "z_physical"
    z_physical_dimension_dimension_exists =
        dimension_exists(nc, physical_name, size(interpolated_physical_z))

    if !z_physical_dimension_dimension_exists
        FT = eltype(interpolated_physical_z)
        dim = NCDatasets.defVar(
            nc,
            physical_name,
            FT,
            (depending_on_dimensions..., name),
        )
        dim.attrib["units"] = "m"
        if length(depending_on_dimensions) == 2
            dim[:, :, :] = interpolated_physical_z
        elseif length(depending_on_dimensions) == 1
            dim[:, :] = interpolated_physical_z
        else
            error("Error in calculating z_physical")
        end
    end
    # We do not output this name because it is not an axis

    return [name]
end

# General hybrid space. This calls both the vertical and horizontal add_space_coordinates_maybe!
# and combines the resulting dictionaries
function target_coordinates(
    space::Spaces.ExtrudedFiniteDifferenceSpace,
    num_points,
    z_sampling_method,
)

    hcoords = vcoords = ()

    num_points_horiz..., num_points_vertic = num_points

    hcoords =
        target_coordinates(Spaces.horizontal_space(space), num_points_horiz)

    vertical_space = Spaces.FiniteDifferenceSpace(
        Spaces.grid(space).vertical_grid,
        Spaces.staggering(space),
    )
    vcoords =
        target_coordinates(vertical_space, num_points_vertic, z_sampling_method)

    hcoords == vcoords == () && error("Found empty space")

    return hcoords, vcoords
end

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace2D,
    domain::Domains.SphereDomain,
    hpts,
)
    # Notice LatLong not LongLat!
    return [Geometry.LatLongPoint(hc2, hc1) for hc1 in hpts[1], hc2 in hpts[2]]
end


# Workaround for https://github.com/CliMA/ClimaCore.jl/issues/1936
#
# Note that here we are assuming that lat is first
Geometry._coordinate(pt::Geometry.LatLongPoint, ::Val{1}) =
    Geometry.LatPoint(pt.lat)
Geometry._coordinate(pt::Geometry.LatLongPoint, ::Val{2}) =
    Geometry.LongPoint(pt.long)

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace2D,
    domain::Domains.RectangleDomain,
    hpts,
)
    XYPointType = typeof(
        Geometry.product_coordinates(
            domain.interval1.coord_max,
            domain.interval2.coord_max,
        ),
    )

    return [XYPointType(hc1, hc2) for hc1 in hpts[1], hc2 in hpts[2]]
end

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace1D,
    domain::Domains.IntervalDomain{PointType},
    hpts,
) where {PointType}
    return [PointType(hc1) for hc1 in hpts]
end

"""
    hcoords_from_horizontal_space(space, domain, hpts)

Prepare the matrix of horizontal coordinates with the correct type according to the given `space`
and `domain` (e.g., `ClimaCore.Geometry.LatLongPoint`s).
"""
function hcoords_from_horizontal_space(space, domain, hpts) end


"""
    default_num_points(space)

Return a tuple with number of points that are optimally suited to interpolate the given
`space`.

"Optimally suited" here means approximately the same as the number of points as the given
`space`.
"""
function default_num_points(space::Spaces.ExtrudedFiniteDifferenceSpace)
    horizontal_space = Spaces.horizontal_space(space)
    num_horz = default_num_points(horizontal_space)

    vertical_space = Spaces.FiniteDifferenceSpace(
        Spaces.grid(space).vertical_grid,
        Spaces.staggering(space),
    )
    num_vert = default_num_points(vertical_space)
    return (num_horz..., num_vert...)
end

# 2D sphere
function default_num_points(space::Spaces.CubedSphereSpectralElementSpace2D)
    # A cubed sphere has 4 panels to cover the range of longitudes, each panel has
    # `num_elements_per_panel` elements, each with `unique_degrees_of_freedom` points. Same
    # for latitudes, except that we need 2 panels to cover from 0 to 180.

    unique_degrees_of_freedom = ClimaCore.Quadratures.unique_degrees_of_freedom(
        Grids.quadrature_style(space),
    )
    num_elements_per_panel = ClimaCore.Meshes.n_elements_per_panel_direction(
        Spaces.topology(space).mesh,
    )
    num_lat = 2 * num_elements_per_panel * unique_degrees_of_freedom
    num_lon = 2num_lat
    return (num_lon, num_lat)
end

# TODO: Maybe move to ClimaCore?
const RectilinearSpectralElementSpace1D = Spaces.SpectralElementSpace1D{
    <:Grids.SpectralElementGrid1D{<:ClimaCore.Topologies.IntervalTopology},
}

# 1D box
function default_num_points(space::RectilinearSpectralElementSpace1D)
    unique_degrees_of_freedom = ClimaCore.Quadratures.unique_degrees_of_freedom(
        Grids.quadrature_style(space),
    )
    return (
        unique_degrees_of_freedom *
        ClimaCore.Meshes.nelements(Spaces.topology(space).mesh),
    )
end

# 2D box
function default_num_points(space::Spaces.RectilinearSpectralElementSpace2D)
    unique_degrees_of_freedom = ClimaCore.Quadratures.unique_degrees_of_freedom(
        Grids.quadrature_style(space),
    )

    return (
        unique_degrees_of_freedom *
        ClimaCore.Meshes.nelements(Spaces.topology(space).mesh.intervalmesh1),
        unique_degrees_of_freedom *
        ClimaCore.Meshes.nelements(Spaces.topology(space).mesh.intervalmesh2),
    )
end

# Column
function default_num_points(space::Spaces.FiniteDifferenceSpace)
    # We always want the center space for interpolation
    cspace = Spaces.center_space(space)
    return (Spaces.nlevels(cspace),)
end
