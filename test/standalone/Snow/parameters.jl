using Test
using ClimaLand
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float64;

@testset "Scalar Parameters" begin
    Δt = 450.0
    default_params = ClimaLand.Snow.SnowParameters(FT, Δt)
    @test default_params.ϵ_snow ==
          LP.get_default_parameter(FT, :snow_emissivity)
    @test default_params.α_snow == ClimaLand.Snow.ConstantAlbedoModel(
        LP.get_default_parameter(FT, :snow_albedo),
    )
    @test default_params.z_0b ==
          LP.get_default_parameter(FT, :snow_scalar_roughness_length)

    α_snow = ClimaLand.Snow.ConstantAlbedoModel(FT(0.7))
    ϵ_snow = FT(0.999)
    overwritten_params = ClimaLand.Snow.SnowParameters(FT, Δt; α_snow, ϵ_snow)
    @test overwritten_params.ϵ_snow == ϵ_snow
    @test overwritten_params.α_snow == α_snow
end
