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

    earth_param_set = LP.LandParameters(FT)
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
