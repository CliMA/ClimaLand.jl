using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy

@testset "PrognosticSoil" begin
    FT = Float64
    α_PAR = 0.2
    α_NIR = 0.4
    soil_driver = ClimaLand.PrognosticSoilConditions{FT}(α_PAR, α_NIR)
    prognostic_land_components = (:canopy, :soil, :soilco2)
    @test Canopy.ground_albedo_PAR(
        Val(prognostic_land_components),
        soil_driver,
        nothing,
        nothing,
        nothing,
    ) == 0.2
    @test Canopy.ground_albedo_NIR(
        Val(prognostic_land_components),
        soil_driver,
        nothing,
        nothing,
        nothing,
    ) == 0.4

end
