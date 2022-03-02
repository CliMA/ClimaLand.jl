using ClimaCore
using Test
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM: Domains
using ClimaLSM.Domains:
    Column,
    PlantHydraulicsDomain,
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
        zmin = FT(1.0)
        zmax = FT(2.0)
        xlim = FT.((0.0, 10.0))
        ylim = FT.((0.0, 1.0))
        zlim = FT.((zmin, zmax))
        nelements = (1, 1, 5)
        shell = SphericalShell(;
            radius = FT(100.0),
            height = FT(30.0),
            nelements = (6, 20),
            npolynomial = 3,
        )
        @test shell.radius == FT(100)
        @test shell.height == FT(30)
        @test shell.nelements == (6, 20)
        @test shell.npolynomial == 3
        shell_coords = coordinates(shell)
        @test eltype(shell_coords) == ClimaCore.Geometry.LatLongZPoint{FT}
        @test typeof(shell_coords) <: ClimaCore.Fields.Field
        @test typeof(shell.space.horizontal_space) <:
              ClimaCore.Spaces.SpectralElementSpace2D
        @test typeof(shell.space) <:
              ClimaCore.Spaces.CenterExtrudedFiniteDifferenceSpace


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
        zmin = FT(1.0)
        zmax = FT(2.0)
        zlim = FT.((zmin, zmax))
        nelements = 5
        domain = LSMSingleColumnDomain(; zlim = zlim, nelements = nelements)
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
        zmin = FT(1.0)
        zmax = FT(2.0)
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
            height = FT(30.0),
            nelements = (6, 20),
            npolynomial = 3,
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


@testset "PlantHydraulicsDomain" begin
    for FT in TestFloatTypes
        Δz = 1.0
        roots = [FT(1.0), FT(2.0), FT(3.0)]
        stem = Int64(5)
        leaves = Int64(4)
        comp_points = Vector(FT(0.5):FT(1.0):FT(8.5))
        top_of_compartments = Vector(FT(0.0):FT(1.0):FT(9.0))
        comp_labels =
            [:stem, :stem, :stem, :stem, :stem, :leaf, :leaf, :leaf, :leaf]
        test_tuple = (root = 1, stem = 2, leaf = 3)
        @test test_tuple[:root] == 1

        testing = PlantHydraulicsDomain(roots, stem, leaves, FT(Δz))
        @test testing.root_depths == [1.0, 2.0, 3.0]
        @test testing.n_stem == 5
        @test testing.n_leaf == 4
        @test testing.compartment_midpoints == comp_points
        @test testing.compartment_surfaces == top_of_compartments
        @test testing.compartment_labels == comp_labels
        # Check that root depth array is monotonic and increasing 
        @test(
            (
                (
                    testing.root_depths[i + 1] - testing.root_depths[i] for
                    i in 1:(length(testing.root_depths) - 1)
                ) .> FT(0)
            ) == Bool.(ones(length(testing.root_depths) - 1))
        )

        for i in 1:5
            @test test_tuple[testing.compartment_labels[i]] == 2
        end
        for i in 6:9
            @test test_tuple[testing.compartment_labels[i]] == 3
        end

        testing2 = PlantHydraulicsDomain(
            roots,
            stem,
            leaves,
            comp_points,
            top_of_compartments,
        )
        @test testing2.root_depths == [1.0, 2.0, 3.0]
        @test testing2.n_stem == 5
        @test testing2.n_leaf == 4
        @test testing2.compartment_midpoints ==
              [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
        @test testing2.compartment_surfaces ==
              [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        @test testing2.compartment_labels == comp_labels
        for i in 1:5
            @test test_tuple[testing2.compartment_labels[i]] == 2
        end
        for i in 6:9
            @test test_tuple[testing2.compartment_labels[i]] == 3
        end

        coords = coordinates(testing)
        @test coords == comp_points
    end
end
