using ClimaCore
using Test
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains:
    Column, RootDomain, HybridBox, Plane, Point, LSMSingleColumnDomain
using ClimaLSM.Domains: coordinates

FT = Float64
zmin = FT(1.0)
zmax = FT(2.0)
xlim = (0.0, 10.0)
ylim = (0.0, 1.0)
zlim = (zmin, zmax)
nelements = (1, 1, 5)

@testset "Box Domain" begin
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
end

@testset "Plane Domain" begin
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
end

@testset "Column Domain" begin
    z_column = Column(; zlim = zlim, nelements = nelements[3])
    column_coords = coordinates(z_column)
    @test z_column.zlim == FT.(zlim)
    @test z_column.nelements[1] == nelements[3]
    @test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
    @test typeof(column_coords) <: ClimaCore.Fields.Field
end

@testset "Point Domain" begin
    point = Point(; z_sfc = zlim[1])
    @test point.z_sfc == zlim[1]
    @test coordinates(point) == [zlim[1]]
end


@testset "LSMSingleColumnDomain" begin
    domain = LSMSingleColumnDomain(; zlim = zlim, nelements = nelements[3])
    point = domain.surface
    column = domain.subsurface

    @test typeof(column) == Column{FT}
    @test typeof(point) == Point{FT}
    @test coordinates(point) == coordinates(domain).surface
    @test parent(coordinates(domain).subsurface) == parent(coordinates(column))
end
