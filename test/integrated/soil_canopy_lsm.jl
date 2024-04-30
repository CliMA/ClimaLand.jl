using Test
import ClimaComms
pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy

@testset "PrognosticSoil" begin
    FT = Float64
    α_PAR = 0.2
    α_NIR = 0.4
    soil_driver = ClimaLand.PrognosticSoil{FT}(α_PAR, α_NIR)
    @test Canopy.ground_albedo_PAR(soil_driver, nothing, nothing, nothing) ==
          0.2
    @test Canopy.ground_albedo_NIR(soil_driver, nothing, nothing, nothing) ==
          0.4
end
