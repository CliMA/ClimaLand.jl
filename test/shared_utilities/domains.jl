import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using NCDatasets
using Test
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
using ClimaLand
using ClimaLand: Domains
using Interpolations
using ClimaLand.Domains:
    Column,
    HybridBox,
    Plane,
    Point,
    SphericalShell,
    SphericalSurface,
    coordinates,
    obtain_surface_space,
    obtain_surface_domain,
    get_Δz,
    top_face_to_surface,
    average_horizontal_resolution_degrees,
    get_lat,
    get_long

FT = Float32
@testset "Clima Core Domains, FT = $FT" begin
    zmin = FT(-1.0)
    zmax = FT(0.0)
    xlim = FT.((0.0, 10.0))
    ylim = FT.((0.0, 1.0))
    zlim = FT.((zmin, zmax))
    nelements = (1, 1, 10)
    radius = FT(100)
    depth = FT(30)
    n_elements_sphere = (6, 20)

    # NOTE: Here we set npoly_sphere to 3, instead of 0 to test that npoly != 0
    # works as expected
    npoly_sphere = 3
    # Spherical Shell
    shell = SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = n_elements_sphere,
        npolynomial = npoly_sphere,
    )
    @test shell.fields.Δz_min == depth / 20
    @test shell.fields.depth == depth
    @test shell.fields.z ==
          ClimaCore.Fields.coordinate_field(shell.space.subsurface).z
    face_space = ClimaCore.Spaces.face_space(shell.space.subsurface)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    @test shell.fields.z_sfc == top_face_to_surface(z_face, shell.space.surface)
    Δz_top, Δz_bottom, Δz = get_Δz(shell.fields.z)
    @test shell.fields.Δz_top == Δz_top
    @test shell.fields.Δz_bottom == Δz_bottom
    @test shell.radius == radius
    @test shell.depth == depth
    @test shell.nelements == n_elements_sphere
    @test shell.npolynomial == npoly_sphere
    shell_coords = coordinates(shell).subsurface
    @test eltype(shell_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
    @test typeof(shell_coords) <: ClimaCore.Fields.Field
    @test typeof(ClimaCore.Spaces.horizontal_space(shell.space.subsurface)) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test typeof(shell.space.subsurface) <:
          ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace
    @test typeof(shell.space.surface) <: ClimaCore.Spaces.SpectralElementSpace2D
    @test obtain_surface_space(shell.space.subsurface) == shell.space.surface
    @test obtain_surface_domain(shell) ==
          SphericalSurface{FT, typeof((; surface = shell.space.surface))}(
        radius,
        n_elements_sphere[1],
        npoly_sphere,
        (; surface = shell.space.surface),
    )
    @test ClimaComms.context(shell) == ClimaComms.context()
    @test ClimaComms.device(shell) == ClimaComms.device()
    shell_stretch = SphericalShell(;
        radius = radius,
        depth = FT(1.0),
        dz_tuple = FT.((0.3, 0.03)),
        nelements = (6, 10),
    )
    shell_coords_stretch = coordinates(shell_stretch).subsurface
    dz =
        Array(parent(shell_coords_stretch.z))[:, end, end, end, end][2:end] .-
        Array(parent(shell_coords_stretch.z))[:, end, end, end, end][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2
    @test shell.fields.Δz_min == minimum(shell.fields.Δz)

    @test get_lat(shell.space.surface) ==
          ClimaCore.Fields.coordinate_field(shell.space.surface).lat
    @test get_long(shell.space.surface) ==
          ClimaCore.Fields.coordinate_field(shell.space.surface).long
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        shell.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        shell.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        shell.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        shell.space.subsurface_face,
    )


    # Spherical Surface
    shell_surface =
        SphericalSurface(; radius = radius, nelements = n_elements_sphere[1])
    @test shell_surface.radius == radius
    @test shell_surface.nelements == n_elements_sphere[1]
    shell_surface_coords = coordinates(shell_surface).surface
    @test eltype(shell_surface_coords) == ClimaCore.Geometry.LatLongPoint{FT}
    @test typeof(shell_surface_coords) <: ClimaCore.Fields.Field
    @test typeof(shell_surface.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D

    @test ClimaComms.context(shell_surface) == ClimaComms.context()
    @test ClimaComms.device(shell_surface) == ClimaComms.device()

    @test average_horizontal_resolution_degrees(shell_surface) ==
          (180 / (2 * n_elements_sphere[1]), 180 / (2n_elements_sphere[1]))
    @test get_lat(shell_surface.space.surface) ==
          ClimaCore.Fields.coordinate_field(shell_surface.space.surface).lat
    @test get_long(shell_surface.space.surface) ==
          ClimaCore.Fields.coordinate_field(shell_surface.space.surface).long


    # HybridBox
    box = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        nelements = nelements,
    )
    @test ClimaComms.context(box) == ClimaComms.context()
    @test ClimaComms.device(box) == ClimaComms.device()
    @test box.fields.depth == zlim[2] - zlim[1]
    @test box.fields.z ==
          ClimaCore.Fields.coordinate_field(box.space.subsurface).z
    face_space = ClimaCore.Spaces.face_space(box.space.subsurface)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    @test box.fields.z_sfc == top_face_to_surface(z_face, box.space.surface)
    Δz_top, Δz_bottom, Δz = get_Δz(box.fields.z)
    @test box.fields.Δz_top == Δz_top
    @test box.fields.Δz_bottom == Δz_bottom
    box_coords = coordinates(box).subsurface
    @test eltype(box_coords) == ClimaCore.Geometry.XYZPoint{FT}
    @test typeof(box_coords) <: ClimaCore.Fields.Field
    @test box.xlim == FT.(xlim)
    @test box.ylim == FT.(ylim)
    @test box.zlim == FT.(zlim)
    @test box.nelements == nelements
    @test box.npolynomial == 0
    @test box.periodic == (true, true)
    @test typeof(ClimaCore.Spaces.horizontal_space(box.space.subsurface)) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test typeof(box.space.subsurface) <:
          ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace
    @test typeof(box.space.surface) <: ClimaCore.Spaces.SpectralElementSpace2D
    @test obtain_surface_space(box.space.subsurface) == box.space.surface
    @test obtain_surface_domain(box) ==
          Plane{FT, typeof((; surface = box.space.surface))}(
        xlim,
        ylim,
        nothing,
        nelements[1:2],
        (true, true),
        0,
        (; surface = box.space.surface),
    )

    @test_throws "Surface space does not have latitude information" get_lat(
        box.space.surface,
    )
    @test_throws "Surface space does not have longitude information" get_long(
        box.space.surface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        box.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        box.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        box.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        box.space.subsurface_face,
    )

    stretch_box = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        dz_tuple = FT.((0.3, 0.03)),
        nelements = nelements,
    )
    box_coords_stretch = coordinates(stretch_box).subsurface
    dz =
        Array(parent(box_coords_stretch.z))[:][2:end] .-
        Array(parent(box_coords_stretch.z))[:][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2


    # Plane
    xy_plane = Plane(;
        xlim = xlim,
        ylim = ylim,
        nelements = nelements[1:2],
        periodic = (true, true),
    )
    @test ClimaComms.context(xy_plane) == ClimaComms.context()
    @test ClimaComms.device(xy_plane) == ClimaComms.device()
    plane_coords = coordinates(xy_plane).surface
    @test eltype(plane_coords) == ClimaCore.Geometry.XYPoint{FT}
    @test typeof(plane_coords) <: ClimaCore.Fields.Field
    @test xy_plane.xlim == FT.(xlim)
    @test xy_plane.ylim == FT.(ylim)
    @test xy_plane.nelements == nelements[1:2]
    @test xy_plane.npolynomial == 0
    @test xy_plane.periodic == (true, true)
    @test typeof(xy_plane.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D

    @test_throws ErrorException average_horizontal_resolution_degrees(xy_plane)
    @test_throws "Surface space does not have latitude information" get_lat(
        xy_plane.space.surface,
    )
    @test_throws "Surface space does not have longitude information" get_long(
        xy_plane.space.surface,
    )


    # Plane latlong
    dxlim = (FT(50_000), FT(80_000))
    dylim = (FT(30_000), FT(40_000))
    longlat = (FT(-118.14452), FT(34.14778))
    radius_earth = FT(6.378e6)
    xlim_longlat = (
        longlat[2] - dylim[1] / FT(2π * radius_earth) * 360,
        longlat[2] + dylim[2] / FT(2π * radius_earth) * 360,
    )
    ylim_longlat = (
        longlat[1] - dxlim[1] / FT(2π * radius_earth) * 360,
        longlat[1] + dxlim[2] / FT(2π * radius_earth) * 360,
    )

    longlat_plane =
        Plane(; xlim = dxlim, ylim = dylim, longlat, nelements = nelements[1:2])
    plane_coords = coordinates(longlat_plane).surface
    @test eltype(plane_coords) == ClimaCore.Geometry.LatLongPoint{FT}
    @test typeof(plane_coords) <: ClimaCore.Fields.Field
    @test longlat_plane.xlim == FT.(xlim_longlat)
    @test longlat_plane.ylim == FT.(ylim_longlat)
    @test longlat_plane.nelements == nelements[1:2]
    @test longlat_plane.npolynomial == 0
    @test longlat_plane.periodic == (false, false)
    @test typeof(longlat_plane.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D

    expected_resolution_x = (xlim_longlat[2] - xlim_longlat[1]) / nelements[1]
    expected_resolution_y = (ylim_longlat[2] - ylim_longlat[1]) / nelements[2]

    @test average_horizontal_resolution_degrees(longlat_plane) ==
          (expected_resolution_x, expected_resolution_y)

    @test Array(parent(get_lat(longlat_plane.space.surface)))[1, 1, 1, 1] ≈
          (xlim_longlat[1] + xlim_longlat[2]) / 2
    @test Array(parent(get_long(longlat_plane.space.surface)))[1, 1, 1, 1] ≈
          (ylim_longlat[1] + ylim_longlat[2]) / 2

    # Box latlong
    longlat_box = HybridBox(;
        xlim = dxlim,
        ylim = dylim,
        zlim = zlim,
        longlat,
        nelements = nelements,
    )
    @test longlat_box.fields.depth == zlim[2] - zlim[1]
    @test longlat_box.fields.z ==
          ClimaCore.Fields.coordinate_field(longlat_box.space.subsurface).z
    face_space = ClimaCore.Spaces.face_space(longlat_box.space.subsurface)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    @test longlat_box.fields.z_sfc ==
          top_face_to_surface(z_face, longlat_box.space.surface)
    Δz_top, Δz_bottom = get_Δz(longlat_box.fields.z)
    @test longlat_box.fields.Δz_top == Δz_top
    @test longlat_box.fields.Δz_bottom == Δz_bottom
    longlat_box_coords = coordinates(longlat_box).subsurface
    @test eltype(longlat_box_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
    @test typeof(longlat_box_coords) <: ClimaCore.Fields.Field
    @test longlat_box.xlim == FT.(xlim_longlat)
    @test longlat_box.ylim == FT.(ylim_longlat)
    @test longlat_box.zlim == FT.(zlim)
    @test longlat_box.nelements == nelements
    @test longlat_box.npolynomial == 0
    @test longlat_box.periodic == (false, false)
    @test typeof(
        ClimaCore.Spaces.horizontal_space(longlat_box.space.subsurface),
    ) <: ClimaCore.Spaces.SpectralElementSpace2D
    @test typeof(longlat_box.space.subsurface) <:
          ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace
    @test typeof(longlat_box.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test obtain_surface_space(longlat_box.space.subsurface) ==
          longlat_box.space.surface
    @test obtain_surface_domain(longlat_box) ==
          Plane{FT, typeof((; surface = longlat_box.space.surface))}(
        xlim_longlat,
        ylim_longlat,
        longlat,
        nelements[1:2],
        (false, false),
        0,
        (; surface = longlat_box.space.surface),
    )

    @test average_horizontal_resolution_degrees(longlat_box) ==
          (expected_resolution_x, expected_resolution_y)

    @test Array(parent(get_lat(longlat_box.space.surface)))[1, 1, 1, 1] ≈
          (xlim_longlat[1] + xlim_longlat[2]) / 2
    @test Array(parent(get_long(longlat_box.space.surface)))[1, 1, 1, 1] ≈
          (ylim_longlat[1] + ylim_longlat[2]) / 2
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_box.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_box.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_box.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_box.space.subsurface_face,
    )

    # Column

    z_column = Column(; zlim = zlim, nelements = nelements[3])
    @test ClimaComms.context(z_column) == ClimaComms.context()
    @test ClimaComms.device(z_column) == ClimaComms.device()
    @test z_column.fields.z ==
          ClimaCore.Fields.coordinate_field(z_column.space.subsurface).z
    @test z_column.fields.depth == zlim[2] - zlim[1]
    face_space = ClimaCore.Spaces.face_space(z_column.space.subsurface)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    @test z_column.fields.z_sfc ==
          top_face_to_surface(z_face, z_column.space.surface)
    Δz_top, Δz_bottom, Δz = get_Δz(z_column.fields.z)
    z = ClimaCore.Fields.coordinate_field(z_column.space.subsurface).z
    @test z_column.fields.z == z
    @test z_column.fields.Δz_top == Δz_top
    @test z_column.fields.Δz_bottom == Δz_bottom
    column_coords = coordinates(z_column).subsurface
    @test z_column.zlim == FT.(zlim)
    @test z_column.nelements[1] == nelements[3]
    @test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
    @test typeof(column_coords) <: ClimaCore.Fields.Field
    @test typeof(z_column.space.subsurface) <:
          ClimaCore.Spaces.CenterFiniteDifferenceSpace
    @test typeof(z_column.space.surface) <: ClimaCore.Spaces.PointSpace
    @test any(
        Array(parent(column_coords))[:][2:end] .-
        Array(parent(column_coords))[:][1:(end - 1)] .≈
        (zmax - zmin) / nelements[3],
    )
    @test obtain_surface_space(z_column.space.subsurface) ==
          z_column.space.surface
    @test obtain_surface_domain(z_column) ==
          Point{FT, typeof((; surface = z_column.space.surface))}(
        zlim[2],
        (; surface = z_column.space.surface),
    )

    @test_throws "Surface space does not have latitude information" get_lat(
        z_column.space.surface,
    )
    @test_throws "Surface space does not have longitude information" get_long(
        z_column.space.surface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        z_column.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        z_column.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        z_column.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        z_column.space.subsurface_face,
    )

    z_column_stretch =
        Column(; zlim = zlim, nelements = 10, dz_tuple = FT.((0.3, 0.03)))
    column_coords = coordinates(z_column_stretch).subsurface
    @test z_column_stretch.zlim == FT.(zlim)
    dz =
        Array(parent(column_coords))[:][2:end] .-
        Array(parent(column_coords))[:][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2
end

@testset "Point Domain, FT = $FT" begin
    zmin = FT(1.0)
    point = Point(; z_sfc = zmin)
    @test ClimaComms.context(point) == ClimaComms.context()
    @test ClimaComms.device(point) == ClimaComms.device()
    @test point.z_sfc == zmin
    point_space = point.space.surface
    @test point_space isa ClimaCore.Spaces.PointSpace
    coords = coordinates(point).surface
    @test coords isa ClimaCore.Fields.Field
    @test eltype(coords) == ClimaCore.Geometry.ZPoint{FT}
    @test Array(parent(point.space.surface.local_geometry.coordinates.z)) ==
          [zmin]
    @test_throws "Surface space does not have latitude information" get_lat(
        point_space,
    )
    @test_throws "Surface space does not have longitude information" get_long(
        point_space,
    )

    # Point with lat/long
    lat = FT(34.14778)
    long = FT(-118.14452)
    longlat_point = Point(z_sfc = zmin, longlat = (long, lat))
    @test ClimaComms.context(longlat_point) == ClimaComms.context()
    @test ClimaComms.device(longlat_point) == ClimaComms.device()
    @test longlat_point.z_sfc == zmin
    point_space = longlat_point.space.surface
    @test point_space isa ClimaCore.Spaces.PointSpace
    coords = coordinates(longlat_point).surface
    @test coords isa ClimaCore.Fields.Field
    @test eltype(coords) == ClimaCore.Geometry.LatLongZPoint{FT}
    @test Array(
        parent(longlat_point.space.surface.local_geometry.coordinates.z),
    ) == [zmin]
    @test Array(
        parent(longlat_point.space.surface.local_geometry.coordinates.lat),
    ) == [lat]
    @test Array(
        parent(longlat_point.space.surface.local_geometry.coordinates.long),
    ) == [long]

    @test get_lat(point_space) == ClimaCore.Fields.zeros(point_space) .+ lat
    @test get_long(point_space) == ClimaCore.Fields.zeros(point_space) .+ long
end


@testset "Column Domain with lat/lon, FT = $FT" begin
    lat = FT(34.14778)
    long = FT(-118.14452)
    zlim = FT.((-1, 0))
    nelements = 10

    # Construct a Column with longlat, no vertical stretching
    longlat_column = Column(; zlim, nelements, longlat = (long, lat))
    @test ClimaComms.context(longlat_column) == ClimaComms.context()
    @test ClimaComms.device(longlat_column) == ClimaComms.device()
    @test longlat_column.fields.z ==
          ClimaCore.Fields.coordinate_field(longlat_column.space.subsurface).z
    @test longlat_column.fields.depth == zlim[2] - zlim[1]
    face_space = ClimaCore.Spaces.face_space(longlat_column.space.subsurface)
    z_face = ClimaCore.Fields.coordinate_field(face_space).z
    @test longlat_column.fields.z_sfc ==
          top_face_to_surface(z_face, longlat_column.space.surface)
    Δz_top, Δz_bottom, Δz = get_Δz(longlat_column.fields.z)
    z = ClimaCore.Fields.coordinate_field(longlat_column.space.subsurface).z
    @test longlat_column.fields.z == z
    @test longlat_column.fields.Δz_top == Δz_top
    @test longlat_column.fields.Δz_bottom == Δz_bottom
    column_coords = coordinates(longlat_column).subsurface
    @test longlat_column.zlim == FT.(zlim)
    @test longlat_column.nelements[1] == nelements
    @test eltype(column_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
    @test typeof(column_coords) <: ClimaCore.Fields.Field
    @test typeof(longlat_column.space.subsurface) <:
          ClimaCore.Spaces.CenterFiniteDifferenceSpace
    @test typeof(longlat_column.space.surface) <: ClimaCore.Spaces.PointSpace
    @test any(
        Array(parent(column_coords))[:][2:end] .-
        Array(parent(column_coords))[:][1:(end - 1)] .≈
        (zlim[2] - zlim[1]) / nelements,
    )
    @test obtain_surface_space(longlat_column.space.subsurface) ==
          longlat_column.space.surface
    @test obtain_surface_domain(longlat_column) ==
          Point{FT, typeof((; surface = longlat_column.space.surface))}(
        zlim[2],
        (; surface = longlat_column.space.surface),
    )
    @test get_lat(longlat_column.space.surface) ==
          ClimaCore.Fields.zeros(longlat_column.space.surface) .+ lat
    @test get_long(longlat_column.space.surface) ==
          ClimaCore.Fields.zeros(longlat_column.space.surface) .+ long
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_column.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_column.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_column.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_column.space.subsurface_face,
    )

    # Construct Column with longlat and vertical stretching
    longlat_column_stretch = Column(;
        zlim,
        nelements,
        longlat = (long, lat),
        dz_tuple = FT.((0.3, 0.03)),
    )
    column_coords = coordinates(longlat_column_stretch).subsurface
    @test longlat_column_stretch.zlim == FT.(zlim)
    dz =
        Array(parent(column_coords.z))[:][2:end] .-
        Array(parent(column_coords.z))[:][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2

    @test get_lat(longlat_column_stretch.space.surface) ==
          ClimaCore.Fields.zeros(longlat_column_stretch.space.surface) .+ lat
    @test get_long(longlat_column_stretch.space.surface) ==
          ClimaCore.Fields.zeros(longlat_column_stretch.space.surface) .+ long
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_column_stretch.space.subsurface,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_column_stretch.space.subsurface,
    )
    @test_throws "`get_lat` is not implemented for subsurface spaces" get_lat(
        longlat_column_stretch.space.subsurface_face,
    )
    @test_throws "`get_long` is not implemented for subsurface spaces" get_long(
        longlat_column_stretch.space.subsurface_face,
    )
end

@testset "use_clm_lowres" begin
    FT = Float64
    # test `use_lowres_clm` on the sphere domain, and then again with a sphere domain with 2x
    # the horizontal resolution. Then repeat with plane domains.
    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (101, 15),
        dz_tuple = FT.((10.0, 0.05)),
    )
    surface_space = domain.space.surface
    @test ClimaLand.Domains.use_lowres_clm(surface_space)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (202, 2),
    )
    surface_space = domain.space.surface
    @test !ClimaLand.Domains.use_lowres_clm(surface_space)
    domain = ClimaLand.Domains.Plane(;
        longlat = (-117.0, 34.0),
        xlim = (0.0, FT(2e6)),
        ylim = (0.0, FT(2e6)),
        nelements = (10, 10),
    )
    surface_space = domain.space.surface
    @test ClimaLand.Domains.use_lowres_clm(surface_space)
    domain = ClimaLand.Domains.Plane(;
        longlat = (-117.0, 34.0),
        xlim = (0.0, FT(2e5)),
        ylim = (0.0, FT(2e5)),
        nelements = (10, 10),
    )
    surface_space = domain.space.surface
    @test !ClimaLand.Domains.use_lowres_clm(surface_space)
end

@testset "read in spatial data at a point" begin
    # Manually access the rooting depth from the CLM data to compare with
    context = ClimaComms.context()
    clm_artifact_path =
        ClimaLand.Artifacts.clm_data_folder_path(; context, lowres = true)
    filepath = joinpath(clm_artifact_path, "vegetation_properties_map.nc")

    NCDataset(filepath) do ncd
        lats = ncd["lat"][:]
        longs = ncd["lon"][:]
        lat_ind = 2
        long_ind = 10
        lat = FT(lats[lat_ind])
        long = FT(longs[long_ind])

        rooting_depth_nc = ncd["rooting_depth"][long_ind, lat_ind]

        # Construct a Point domain at a longlat that appears in the dataset (not testing interpolation)
        zmin = FT(0)
        nelements = 10
        longlat = (long, lat)
        zlim = FT.((-1, 0))
        column = ClimaLand.Domains.Column(; zlim, nelements, longlat)
        surface_space = column.space.surface

        # Get the rooting depth at this lat/lon using ClimaUtilities
        rooting_depth_pt = ClimaLand.Canopy.clm_rooting_depth(surface_space)

        @test all(Array(parent(rooting_depth_pt)) .== FT(rooting_depth_nc))
    end
end

@testset "read in spatial data at a column" begin
    # Manually access 3D soil parameter data to compare with
    context = ClimaComms.context()

    soil_params_artifact_path =
        ClimaLand.Artifacts.soil_params_artifact_folder_path(; context)
    filepath = joinpath(
        soil_params_artifact_path,
        "residual_map_gupta_etal2020_1.0x1.0x4.nc",
    )

    NCDataset(filepath) do ncd
        lats = ncd["lat"][:]
        longs = ncd["lon"][:]
        zs = ncd["z"][:]
        lat_ind = 90
        long_ind = 100
        lat = FT(lats[lat_ind])
        long = FT(longs[long_ind])

        θ_r_nc = ncd["θ_r"][long_ind, lat_ind, :]

        # Construct a Column domain at a longlat that appears in the dataset
        # The dz_tuple is chosen so that the depths of the column will map to
        # the 4 layers in the data using nearest neighbor interpolation
        nelements = length(zs)
        longlat = (long, lat)
        zlim = FT.((-1, 0))
        column = ClimaLand.Domains.Column(;
            zlim,
            nelements,
            longlat,
            dz_tuple = (FT(0.35), FT(0.25)),
        )
        subsurface_space = column.space.subsurface

        # Get the residual water content at this lat/lon using ClimaUtilities
        # use nearest neighbor interpolation to map the column to data in the depth
        # dimension
        (; θ_r) = ClimaLand.Soil.soil_vangenuchten_parameters(
            subsurface_space,
            FT;
            interpolation_method = Interpolations.Constant(),
        )

        @test all(Array(parent(θ_r)) .== FT.(θ_r_nc[:]))
    end
end
