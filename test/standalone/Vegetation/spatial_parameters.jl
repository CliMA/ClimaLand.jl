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
@testset "Spatially varying canopy parameters" begin

    regridder_type = :InterpolationsRegridder
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    context = ClimaComms.context()
    ClimaComms.init(context)

    toml_dict = LP.create_toml_dict(FT)
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

    g1 = ClimaLand.Canopy.clm_medlyn_g1(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    @test axes(g1) == surface_space
    rooting_depth = ClimaLand.Canopy.clm_rooting_depth(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    @test axes(rooting_depth) == surface_space
    (; is_c3, Vcmax25) = ClimaLand.Canopy.clm_photosynthesis_parameters(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )
    @test axes(is_c3) == surface_space
    @test axes(Vcmax25) == surface_space
    (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf) =
        ClimaLand.Canopy.clm_canopy_radiation_parameters(
            surface_space;
            regridder_type = regridder_type,
            extrapolation_bc = extrapolation_bc,
        )
    @test axes(Ω) == surface_space
    @test axes(G_Function.χl) == surface_space
    @test axes(α_PAR_leaf) == surface_space
    @test axes(α_NIR_leaf) == surface_space
    @test axes(τ_PAR_leaf) == surface_space
    @test axes(τ_NIR_leaf) == surface_space
end

@testset "Spatially varying canopy height" begin
    regridder_type = :InterpolationsRegridder
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    context = ClimaComms.context()
    ClimaComms.init(context)

    toml_dict = LP.create_toml_dict(FT)
    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = (101, 15),
        dz_tuple = FT.((10.0, 0.05)),
    )
    surface_space = domain.space.surface

    # Test clm_canopy_height function (reads raw CLM data)
    raw_height = ClimaLand.Canopy.clm_canopy_height(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
    )

    # Check that field is on the correct space
    @test axes(raw_height) == surface_space

    # Check that heights are positive
    @test all(Array(parent(raw_height)) .>= 0)

    # Check that some heights are non-zero (there should be vegetation somewhere)
    @test any(Array(parent(raw_height)) .> 0)

    # Check that maximum height is reasonable (should be < 50m based on CLM data)
    @test maximum(Array(parent(raw_height))) < 50.0

    # Try capping the height
    capped_height = ClimaLand.Canopy.clm_canopy_height(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
        max_height = FT(1.0),
    )
    @test maximum(Array(parent(capped_height))) <= 1.0
end
