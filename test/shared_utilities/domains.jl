using ClimaCore
using Test
using ClimaLSM
using ClimaLSM: Domains
using ClimaLSM.Domains:
    Column,
    HybridBox,
    Plane,
    Point,
    SphericalShell,
    SphericalSurface,
    coordinates,
    obtain_surface_space,
    obtain_face_space,
    obtain_surface_domain

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
    npoly_sphere = 3
    # Spherical Shell
    shell = SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = n_elements_sphere,
        npolynomial = npoly_sphere,
    )
    @test shell.radius == radius
    @test shell.depth == depth
    @test shell.nelements == n_elements_sphere
    @test shell.npolynomial == npoly_sphere
    shell_coords = coordinates(shell).subsurface
    @test eltype(shell_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
    @test typeof(shell_coords) <: ClimaCore.Fields.Field
    @test typeof(shell.space.subsurface.horizontal_space) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test typeof(shell.space.subsurface) <:
          ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace
    @test typeof(shell.space.surface) <: ClimaCore.Spaces.SpectralElementSpace2D
    @test obtain_surface_space(shell.space.subsurface) == shell.space.surface
    @test obtain_surface_domain(shell) == SphericalSurface{FT}(
        radius,
        n_elements_sphere[1],
        npoly_sphere,
        (; surface = shell.space.surface),
    )
    shell_stretch = SphericalShell(;
        radius = radius,
        depth = FT(1.0),
        dz_tuple = FT.((0.3, 0.03)),
        nelements = (6, 10),
        npolynomial = 3,
    )
    shell_coords_stretch = coordinates(shell_stretch).subsurface
    dz =
        parent(shell_coords_stretch.z)[:, 1, 4, 1, 216][2:end] .-
        parent(shell_coords_stretch.z)[:, 1, 4, 1, 216][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2

    # Spherical Surface
    shell_surface = SphericalSurface(;
        radius = radius,
        nelements = n_elements_sphere[1],
        npolynomial = npoly_sphere,
    )
    @test shell_surface.radius == radius
    @test shell_surface.nelements == n_elements_sphere[1]
    @test shell_surface.npolynomial == npoly_sphere
    shell_surface_coords = coordinates(shell_surface).surface
    @test eltype(shell_surface_coords) == ClimaCore.Geometry.LatLongPoint{FT}
    @test typeof(shell_surface_coords) <: ClimaCore.Fields.Field
    @test typeof(shell_surface.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D


    # HybridBox
    xyz_column_box = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        nelements = nelements,
        npolynomial = 0,
    )
    box_coords = coordinates(xyz_column_box).subsurface
    @test eltype(box_coords) == ClimaCore.Geometry.XYZPoint{FT}
    @test typeof(box_coords) <: ClimaCore.Fields.Field
    @test xyz_column_box.xlim == FT.(xlim)
    @test xyz_column_box.ylim == FT.(ylim)
    @test xyz_column_box.zlim == FT.(zlim)
    @test xyz_column_box.nelements == nelements
    @test xyz_column_box.npolynomial == 0
    @test xyz_column_box.periodic == (true, true)
    @test typeof(xyz_column_box.space.subsurface.horizontal_space) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test typeof(xyz_column_box.space.subsurface) <:
          ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace
    @test typeof(xyz_column_box.space.surface) <:
          ClimaCore.Spaces.SpectralElementSpace2D
    @test obtain_surface_space(xyz_column_box.space.subsurface) ==
          xyz_column_box.space.surface
    @test obtain_surface_domain(xyz_column_box) == Plane{FT}(
        xlim,
        ylim,
        nelements[1:2],
        (true, true),
        0,
        (; surface = xyz_column_box.space.surface),
    )

    xyz_stretch_column_box = HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        dz_tuple = FT.((0.3, 0.03)),
        nelements = nelements,
        npolynomial = 0,
    )
    box_coords_stretch = coordinates(xyz_stretch_column_box).subsurface
    dz =
        parent(box_coords_stretch.z)[:][2:end] .-
        parent(box_coords_stretch.z)[:][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2

    # Plane
    xy_plane = Plane(;
        xlim = xlim,
        ylim = ylim,
        nelements = nelements[1:2],
        periodic = (true, true),
        npolynomial = 0,
    )
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


    # Column

    z_column = Column(; zlim = zlim, nelements = nelements[3])
    column_coords = coordinates(z_column).subsurface
    @test z_column.zlim == FT.(zlim)
    @test z_column.nelements[1] == nelements[3]
    @test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
    @test typeof(column_coords) <: ClimaCore.Fields.Field
    @test typeof(z_column.space.subsurface) <:
          ClimaCore.Spaces.CenterFiniteDifferenceSpace
    @test typeof(z_column.space.surface) <: ClimaCore.Spaces.PointSpace
    @test any(
        parent(column_coords)[:][2:end] .-
        parent(column_coords)[:][1:(end - 1)] .≈ (zmax - zmin) / nelements[3],
    )
    @test obtain_surface_space(z_column.space.subsurface) ==
          z_column.space.surface
    @test obtain_surface_domain(z_column) ==
          Point{FT}(zlim[2], (; surface = z_column.space.surface))
    z_column_stretch =
        Column(; zlim = zlim, nelements = 10, dz_tuple = FT.((0.3, 0.03)))
    column_coords = coordinates(z_column_stretch).subsurface
    @test z_column_stretch.zlim == FT.(zlim)
    dz =
        parent(column_coords)[:][2:end] .- parent(column_coords)[:][1:(end - 1)]
    @test abs(dz[1] - 0.3) < 1e-1
    @test abs(dz[end] - 0.03) < 1e-2
end

@testset "Point Domain, FT = $FT" begin
    zmin = FT(1.0)
    point = Point(; z_sfc = zmin)
    @test point.z_sfc == zmin
    point_space = point.space.surface
    @test point_space isa ClimaCore.Spaces.PointSpace
    coords = coordinates(point).surface
    @test coords isa ClimaCore.Fields.Field
    @test eltype(coords) == ClimaCore.Geometry.ZPoint{FT}
    @test ClimaCore.Fields.field_values(coords)[] ==
          ClimaCore.Geometry.ZPoint(zmin)
end
