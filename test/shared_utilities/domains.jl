using ClimaCore
using Test
using ClimaLSM
using ClimaLSM: Domains
using ClimaLSM.Domains:
    Column,
    HybridBox,
    Plane,
    Point,
    LSMSingleColumnDomain,
    LSMMultiColumnDomain,
    LSMSphericalShellDomain,
    SphericalShell,
    SphericalSurface,
    coordinates,
    obtain_surface_space,
    obtain_face_space

TestFloatTypes = (Float32, Float64)

@testset "Clima Core Domains" begin
    for FT in TestFloatTypes
        zmin = FT(-1.0)
        zmax = FT(0.0)
        xlim = FT.((0.0, 10.0))
        ylim = FT.((0.0, 1.0))
        zlim = FT.((zmin, zmax))
        nelements = (1, 1, 10)
        shell = SphericalShell(;
            radius = FT(100.0),
            depth = FT(30.0),
            nelements = (6, 20),
            npolynomial = 3,
        )
        @test shell.radius == FT(100)
        @test shell.depth == FT(30)
        @test shell.nelements == (6, 20)
        @test shell.npolynomial == 3
        shell_coords = coordinates(shell)
        @test eltype(shell_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
        @test typeof(shell_coords) <: ClimaCore.Fields.Field
        @test typeof(shell.space.horizontal_space) <:
              ClimaCore.Spaces.SpectralElementSpace2D
        @test typeof(shell.space) <:
              ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace

        shell_stretch = SphericalShell(;
            radius = FT(100.0),
            depth = FT(1.0),
            dz_tuple = FT.((0.3, 0.03)),
            nelements = (6, 10),
            npolynomial = 3,
        )
        shell_coords_stretch = coordinates(shell_stretch)
        dz =
            parent(shell_coords_stretch.z)[:, 1, 4, 1, 216][2:end] .-
            parent(shell_coords_stretch.z)[:, 1, 4, 1, 216][1:(end - 1)]
        @test abs(dz[1] - 0.3) < 1e-1
        @test abs(dz[end] - 0.03) < 1e-2


        shell_surface = SphericalSurface(;
            radius = FT(100.0),
            nelements = 6,
            npolynomial = 3,
        )
        @test shell_surface.radius == FT(100)
        @test shell_surface.nelements == 6
        @test shell_surface.npolynomial == 3
        shell_surface_coords = coordinates(shell_surface)
        @test eltype(shell_surface_coords) ==
              ClimaCore.Geometry.LatLongPoint{FT}
        @test typeof(shell_surface_coords) <: ClimaCore.Fields.Field
        @test typeof(shell_surface.space) <:
              ClimaCore.Spaces.SpectralElementSpace2D

        xy_plane = Plane(;
            xlim = xlim,
            ylim = ylim,
            nelements = nelements[1:2],
            periodic = (true, true),
            npolynomial = 0,
        )

        xyz_column_box = HybridBox(;
            xlim = xlim,
            ylim = ylim,
            zlim = zlim,
            nelements = nelements,
            npolynomial = 0,
        )
        box_coords = coordinates(xyz_column_box)
        @test eltype(box_coords) == ClimaCore.Geometry.XYZPoint{FT}
        @test typeof(box_coords) <: ClimaCore.Fields.Field
        @test xyz_column_box.xlim == FT.(xlim)
        @test xyz_column_box.ylim == FT.(ylim)
        @test xyz_column_box.zlim == FT.(zlim)
        @test xyz_column_box.nelements == nelements
        @test xyz_column_box.npolynomial == 0
        @test xyz_column_box.periodic == (true, true)
        @test typeof(xyz_column_box.space.horizontal_space) <:
              ClimaCore.Spaces.SpectralElementSpace2D
        @test typeof(xyz_column_box.space) <:
              ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace

        xyz_stretch_column_box = HybridBox(;
            xlim = xlim,
            ylim = ylim,
            zlim = zlim,
            dz_tuple = FT.((0.3, 0.03)),
            nelements = nelements,
            npolynomial = 0,
        )
        box_coords_stretch = coordinates(xyz_stretch_column_box)
        dz =
            parent(box_coords_stretch.z)[:][2:end] .-
            parent(box_coords_stretch.z)[:][1:(end - 1)]
        @test abs(dz[1] - 0.3) < 1e-1
        @test abs(dz[end] - 0.03) < 1e-2
        xy_plane = Plane(;
            xlim = xlim,
            ylim = ylim,
            nelements = nelements[1:2],
            periodic = (true, true),
            npolynomial = 0,
        )
        plane_coords = coordinates(xy_plane)
        @test eltype(plane_coords) == ClimaCore.Geometry.XYPoint{FT}
        @test typeof(plane_coords) <: ClimaCore.Fields.Field
        @test xy_plane.xlim == FT.(xlim)
        @test xy_plane.ylim == FT.(ylim)
        @test xy_plane.nelements == nelements[1:2]
        @test xy_plane.npolynomial == 0
        @test xy_plane.periodic == (true, true)
        @test typeof(xy_plane.space) <: ClimaCore.Spaces.SpectralElementSpace2D

        z_column = Column(; zlim = zlim, nelements = nelements[3])
        column_coords = coordinates(z_column)
        @test z_column.zlim == FT.(zlim)
        @test z_column.nelements[1] == nelements[3]
        @test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
        @test typeof(column_coords) <: ClimaCore.Fields.Field
        @test typeof(z_column.space) <:
              ClimaCore.Spaces.CenterFiniteDifferenceSpace
        @test any(
            parent(column_coords)[:][2:end] .-
            parent(column_coords)[:][1:(end - 1)] .≈
            (zmax - zmin) / nelements[3],
        )

        z_column_stretch =
            Column(; zlim = zlim, nelements = 10, dz_tuple = FT.((0.3, 0.03)))
        column_coords = coordinates(z_column_stretch)
        @test z_column_stretch.zlim == FT.(zlim)
        dz =
            parent(column_coords)[:][2:end] .-
            parent(column_coords)[:][1:(end - 1)]
        @test abs(dz[1] - 0.3) < 1e-1
        @test abs(dz[end] - 0.03) < 1e-2
    end
end


@testset "Point Domain" begin
    for FT in TestFloatTypes
        zmin = FT(1.0)
        point = Point(; z_sfc = zmin)
        @test point.z_sfc == zmin
        point_space = point.space
        @test point_space isa ClimaCore.Spaces.PointSpace
        coords = coordinates(point)
        @test coords isa ClimaCore.Fields.Field
        @test eltype(coords) == ClimaCore.Geometry.ZPoint{FT}
        @test ClimaCore.Fields.field_values(coords)[] ==
              ClimaCore.Geometry.ZPoint(zmin)
    end
end


@testset "LSMSingleColumnDomain" begin
    for FT in TestFloatTypes
        zmin = FT(-1.0)
        zmax = FT(0.0)
        zlim = FT.((zmin, zmax))
        nelements = 5
        domain = LSMSingleColumnDomain(;
            zlim = zlim,
            nelements = nelements,
            dz_tuple = FT.((0.3, 0.03)),
        )
        point = domain.surface
        column = domain.subsurface

        @test typeof(column) == Column{FT, typeof(column.space)}
        @test typeof(point) == Point{FT, typeof(point.space)}
        @test coordinates(point) === coordinates(domain).surface
        @test coordinates(domain).subsurface === coordinates(column)
        @test obtain_surface_space(column.space) === point.space
    end
end


@testset "LSMMultiColumnDomain" begin
    for FT in TestFloatTypes
        zmin = FT(-1.0)
        zmax = FT(0.0)
        zlim = FT.((zmin, zmax))
        xlim = FT.((0.0, 10.0))
        ylim = FT.((0.0, 1.0))
        zlim = FT.((zmin, zmax))
        nelements = (1, 1, 5)
        npolynomial = 2
        domain = LSMMultiColumnDomain(;
            xlim = xlim,
            ylim = ylim,
            zlim = zlim,
            nelements = nelements,
            periodic = (true, true),
            npolynomial = npolynomial,
            dz_tuple = FT.((0.3, 0.03)),
        )
        plane = domain.surface
        box = domain.subsurface

        @test typeof(box) == HybridBox{FT, typeof(box.space)}
        @test typeof(plane) == Plane{FT, typeof(plane.space)}
        @test coordinates(plane) === coordinates(domain).surface
        @test coordinates(domain).subsurface === coordinates(box)
        @test obtain_surface_space(box.space) === plane.space
    end
end



@testset "LSMSphericalShellDomain" begin
    for FT in TestFloatTypes
        domain = LSMSphericalShellDomain(;
            radius = FT(100.0),
            depth = FT(30.0),
            nelements = (6, 20),
            npolynomial = 3,
            dz_tuple = FT.((5.0, 0.3)),
        )
        surf = domain.surface
        shell = domain.subsurface

        @test typeof(shell) == SphericalShell{FT, typeof(shell.space)}
        @test typeof(surf) == SphericalSurface{FT, typeof(surf.space)}
        @test coordinates(surf) === coordinates(domain).surface
        @test coordinates(domain).subsurface === coordinates(shell)
        @test obtain_surface_space(shell.space) === surf.space
    end
end
