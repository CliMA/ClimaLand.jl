using ClimaLSM.Regridder: regrid_netcdf_to_field
using ClimaLSM.FileReader: PrescribedDataStatic
using ClimaLSM.Bucket:
    BulkAlbedoStatic,
    bareground_albedo_dataset_path,
    set_initial_parameter_field!
using ClimaLSM.Domains: SphericalSurface
import ClimaLSM
using ClimaComms
using ClimaCore
using Test


@testset "Spatially varying map - regrid to field" begin
    FT = Float32
    path = bareground_albedo_dataset_path()
    regrid_dirpath =
        joinpath(pkgdir(ClimaLSM), "test/standalone/Bucket/regridder_tmpfiles")
    mkpath(regrid_dirpath)

    varname = "sw_alb"
    albedo = PrescribedDataStatic(path, regrid_dirpath, varname)

    surface_domain =
        SphericalSurface(; radius = FT(1), nelements = 2, npolynomial = 3)
    boundary_space = surface_domain.space.surface
    comms = boundary_space.topology.context
    field = regrid_netcdf_to_field(
        FT,
        regrid_dirpath,
        comms,
        path,
        varname,
        boundary_space,
    )
    @test axes(field) == boundary_space


    p = (; :bucket => (; :α_sfc => ClimaCore.Fields.zeros(boundary_space)))
    set_initial_parameter_field!(
        BulkAlbedoStatic{FT}(FT(0.08), albedo),
        p,
        ClimaCore.Fields.coordinate_field(boundary_space),
    )
    @test p.bucket.α_sfc == field
    rm(regrid_dirpath, recursive = true)

end
