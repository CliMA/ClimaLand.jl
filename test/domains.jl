using ClimaCore
using Test
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM: Domains
using ClimaLSM.Domains:
    Column,
    RootDomain,
    HybridBox,
    Plane,
    Point,
    LSMSingleColumnDomain,
    SphericalShell,
    coordinates

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

        z_column = Column(; zlim = zlim, nelements = nelements[3])
        column_coords = coordinates(z_column)
        @test z_column.zlim == FT.(zlim)
        @test z_column.nelements[1] == nelements[3]
        @test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
        @test typeof(column_coords) <: ClimaCore.Fields.Field
    end
end


@testset "Point Domain" begin
    for FT in TestFloatTypes
        zmin = FT(1.0)
        point = Point(; z_sfc = zmin)
        @test point.z_sfc == zmin
        point_space, _ = Domains.make_function_space(point)
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

        @test typeof(column) == Column{FT}
        @test typeof(point) == Point{FT}
        @test parent(coordinates(point)) == parent(coordinates(domain).surface)
        @test parent(coordinates(domain).subsurface) ==
              parent(coordinates(column))
    end
end

TestFloatTypes = (Float32, Float64)
@testset "RootDomain" begin
    for FT in TestFloatTypes
        d_z = 1.0
        roots = [1.0, 2.0, 3.0]
        stem = Int64(5)
        leaves = Int64(4)
        comp_heights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        comp_labels =
            [:stem, :stem, :stem, :stem, :stem, :leaf, :leaf, :leaf, :leaf]
        test_tuple = (root = 1, stem = 2, leaf = 3)
        @test test_tuple[:root] == 1

        testing = RootDomain(roots, stem, leaves, d_z)
        @test testing.root_depths == [1.0, 2.0, 3.0]
        @test testing.n_stem == 5
        @test testing.n_leaf == 4
        @test testing.compartment_heights == comp_heights
        @test testing.compartment_labels == comp_labels
        for i in 1:5
            @test test_tuple[testing.compartment_labels[i]] == 2
        end
        for i in 6:9
            @test test_tuple[testing.compartment_labels[i]] == 3
        end

        testing2 = RootDomain(roots, stem, leaves, comp_heights)
        @test testing2.root_depths == [1.0, 2.0, 3.0]
        @test testing2.n_stem == 5
        @test testing2.n_leaf == 4
        @test testing2.compartment_heights ==
              [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        @test testing2.compartment_labels == comp_labels
        for i in 1:5
            @test test_tuple[testing2.compartment_labels[i]] == 2
        end
        for i in 6:9
            @test test_tuple[testing2.compartment_labels[i]] == 3
        end

        testing3 = RootDomain(roots, stem, leaves, comp_heights, comp_labels)
        @test testing3.root_depths == [1.0, 2.0, 3.0]
        @test testing3.n_stem == 5
        @test testing3.n_leaf == 4
        @test testing3.compartment_heights ==
              [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        @test testing3.compartment_labels ==
              [:stem, :stem, :stem, :stem, :stem, :leaf, :leaf, :leaf, :leaf]
        for i in 1:5
            @test test_tuple[testing3.compartment_labels[i]] == 2
        end
        for i in 6:9
            @test test_tuple[testing3.compartment_labels[i]] == 3
        end

        coords = coordinates(testing)
        @test coords == comp_heights
    end
end
