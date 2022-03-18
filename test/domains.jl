using ClimaLSM.Domains: coordinates
zmin = FT(1.0)
zmax = FT(2.0)
xlim = (0.0, 10.0)
ylim = (0.0, 1.0)
zlim = (zmin, zmax)
nelements = (1, 1, 5)
xyz_column_box = HybridBox(
    Float64;
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

xy_plane = Plane{FT}(xlim, ylim, nelements[1:2], (true, true), 0)
plane_coords = coordinates(xy_plane)
@test eltype(plane_coords) == ClimaCore.Geometry.XYPoint{FT}
@test typeof(plane_coords) <: ClimaCore.Fields.Field
@test xy_plane.xlim == FT.(xlim)
@test xy_plane.ylim == FT.(ylim)
@test xy_plane.nelements == nelements[1:2]
@test xy_plane.npolynomial == 0
@test xy_plane.periodic == (true, true)


z_column = Column(FT, zlim = zlim, nelements = nelements[3])
column_coords = coordinates(z_column)
@test z_column.zlim == FT.(zlim)
@test z_column.nelements == nelements[3]
@test eltype(column_coords) == ClimaCore.Geometry.ZPoint{FT}
@test typeof(column_coords) <: ClimaCore.Fields.Field
