using Test

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float64;

@testset "Spatially varying soil parameters" begin
    regridder_type = :InterpolationsRegridder
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    context = ClimaComms.context()
    ClimaComms.init(context)

    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (101, 15),
        dz_tuple = FT.((10.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    soil_albedo_params = ClimaLand.Soil.clm_soil_albedo_parameters(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    composition_params = ClimaLand.Soil.soil_composition_parameters(
        subsurface_space,
        FT;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    vg_params = ClimaLand.Soil.soil_vangenuchten_parameters(
        subsurface_space,
        FT;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    for p in soil_albedo_params
        @test axes(p) == surface_space
    end
    for p in [composition_params..., vg_params...]
        @test axes(p) == subsurface_space
    end

    hcm = vg_params.hydrology_cm
    @test eltype(hcm) == ClimaLand.Soil.vanGenuchten{FT}
end

# this tests if the upstream ClimaArtifacts generation for the low res (1 deg x 1 deg)
# van Genuchten param maps is creating depth-varying maps as expected
@testset "Depth varying van Genuchten parameters" begin
    regridder_type = :InterpolationsRegridder
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )

    # construct a column domain at the US-Ha1 and check that the vG params are depth varying
    domain = ClimaLand.Domains.Column(;
        zlim = (FT(-1.0), FT(0.0)),
        nelements = 10,
        dz_tuple = (FT(0.2), FT(0.05)),
        longlat = (FT(-72.1715), FT(42.5378)),
    )

    vg_params = ClimaLand.Soil.soil_vangenuchten_parameters(
        domain.space.subsurface,
        FT;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )

    # check that vG params are not all the same
    for p in vg_params
        # this gives a nelements x (n) matrix
        vals = parent(p)
        @test size(vals, 1) > 1 && any(vals .!= @view vals[1:1, :])
    end
end

@testset "Scalar Parameters" begin
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = vg_α, n = vg_n)
    θ_r = FT(0.1)
    S_c = hcm.S_c

    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    default_params = ClimaLand.Soil.EnergyHydrologyParameters(
        toml_dict;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm = hcm,
        K_sat,
        S_s,
        θ_r,
    )
    @test default_params.emissivity ==
          LP.get_default_parameter(FT, :emissivity_bare_soil)
    @test default_params.z_0m ==
          LP.get_default_parameter(FT, :soil_momentum_roughness_length)
    @test default_params.z_0b ==
          LP.get_default_parameter(FT, :soil_scalar_roughness_length)

    overwritten_params = ClimaLand.Soil.EnergyHydrologyParameters(
        toml_dict;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm = hcm,
        K_sat,
        S_s,
        θ_r,
        emissivity = FT(0.8),
        z_0m = FT(0.1),
        z_0b = FT(0.1),
    )
    @test overwritten_params.emissivity == FT(0.8)
    @test overwritten_params.z_0m == FT(0.1)
    @test overwritten_params.z_0b == FT(0.1)

end
