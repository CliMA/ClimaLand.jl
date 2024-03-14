using Test

import ClimaLand
import ClimaLand: Regridders

@testset "InterpolationsRegridder" begin

    lon, lat = collect(0.0:1:360), collect(-90.0:1:90)
    dimensions = (lon, lat)
    sized = (361, 181)
    data_lat = zeros(sized)
    data_lon = zeros(sized)
    for i in 1:length(lon)
        data_lat[i, :] .= lat
    end
    for i in 1:length(lat)
        data_lon[:, i] .= lon
    end

    for FT in (Float32, Float64)
        target_space =
            ClimaLand.Domains.SphericalShell(;
                radius = FT(6731e3),
                depth = FT(1.0),
                nelements = (4, 4),
                npolynomial = 4,
            ).space.surface

        reg = Regridders.InterpolationsRegridder(target_space)

        regridded_lat = Regridders.regrid(reg, data_lat, dimensions)
        regridded_lon = Regridders.regrid(reg, data_lon, dimensions)

        coordinates = ClimaCore.Fields.coordinate_field(target_space)

        # Compute max err
        err_lat = coordinates.lat .- regridded_lat
        err_lon = coordinates.long .- regridded_lon

        @test maximum(err_lat) < 1e-5
        @test maximum(err_lon) < 1e-5
    end

    @test_throws ErrorException Regridders.InterpolationsRegridder(
        ClimaLand.Domains.HybridBox(;
            xlim = (-1.0, 1.0),
            ylim = (-1.0, 1.0),
            zlim = (-1.0, 1.0),
            nelements = (3, 3, 3),
            npolynomial = 0,
        ).space.surface,
    )
end
