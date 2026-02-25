using Test

# Ensure extensions are loaded for SpaceVaryingInput
using ClimaCore
using NCDatasets

import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using Dates

using ClimaLand
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

const FT = Float64

@testset "PrescribedBiomassModel with scalar height" begin
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
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    # Create LAI input
    LAI_fun = (t) -> FT(5.0)
    LAI = TimeVaryingInput(LAI_fun)

    # Test with scalar height
    scalar_height = FT(2.5)
    biomass = Canopy.PrescribedBiomassModel{FT}(
        surface_domain,
        LAI,
        toml_dict;
        height = scalar_height,
    )

    # Verify that the biomass model was created correctly
    @test biomass.height == scalar_height
    @test typeof(biomass.height) == FT
end

@testset "PrescribedBiomassModel with Field height" begin
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
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    # Create LAI input
    LAI_fun = (t) -> FT(5.0)
    LAI = TimeVaryingInput(LAI_fun)

    # Create a spatially-varying height field
    coords = ClimaCore.Fields.coordinate_field(surface_space)
    field_height = coords.lat .* 0 .+ FT(3.0)  # Constant field of 3m for testing

    # Test with Field height (main-branch pattern)
    biomass = Canopy.PrescribedBiomassModel{FT}(
        surface_domain,
        LAI,
        toml_dict;
        height = field_height,
    )

    # Verify that the biomass model was created correctly
    @test biomass.height isa ClimaCore.Fields.Field
    @test axes(biomass.height) == surface_space
    @test all(Array(parent(biomass.height)) .== FT(3.0))
end

@testset "CanopyModel integration with spatially-varying height" begin
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
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    # Set up forcing
    start_date = DateTime(2005)

    u_atmos = t -> 10  # m/s
    liquid_precip = (t) -> 0  # m
    snow_precip = (t) -> 0  # m
    T_atmos = t -> 290  # K
    q_atmos = t -> 0.001  # kg/kg
    P_atmos = t -> 1e5  # Pa
    h_atmos = FT(10)  # m - atmospheric reference height
    c_atmos = (t) -> 4.11e-4  # mol/mol

    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput(liquid_precip),
        TimeVaryingInput(snow_precip),
        TimeVaryingInput(T_atmos),
        TimeVaryingInput(u_atmos),
        TimeVaryingInput(q_atmos),
        TimeVaryingInput(P_atmos),
        start_date,
        h_atmos,
        toml_dict;
        c_co2 = TimeVaryingInput(c_atmos),
    )

    SW_d = (t) -> 300  # W/m^2
    LW_d = (t) -> 250  # W/m^2
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(SW_d),
        TimeVaryingInput(LW_d),
        start_date,
    )

    ground = ClimaLand.PrescribedGroundConditions{FT}()
    forcing = (; atmos, radiation, ground)

    # Create LAI input
    LAI_fun = (t) -> FT(5.0)
    LAI = TimeVaryingInput(LAI_fun)

    # Read and cap spatially-varying canopy height (two-step pattern)
    capped_height = Canopy.clm_canopy_height(
        surface_space;
        regridder_type = regridder_type,
        extrapolation_bc = extrapolation_bc,
        max_height = FT(8),
    )

    # Create biomass model with spatially-varying height
    biomass = Canopy.PrescribedBiomassModel{FT}(
        surface_domain,
        LAI,
        toml_dict;
        height = capped_height,
    )

    # Create canopy model with custom biomass
    canopy = Canopy.CanopyModel{FT}(
        surface_domain,
        forcing,
        LAI,
        toml_dict;
        biomass = biomass,
    )

    # Verify that the canopy model was created correctly
    @test canopy.biomass.height isa ClimaCore.Fields.Field
    @test axes(canopy.biomass.height) == surface_space
    z_0min = canopy.boundary_conditions.turbulent_flux_parameterization.z_0min
    @test canopy.boundary_conditions.turbulent_flux_parameterization.displ ==
          toml_dict["canopy_d_coeff"] .* capped_height
    @test canopy.boundary_conditions.turbulent_flux_parameterization.z_0m ==
          toml_dict["canopy_z_0m_coeff"] .* capped_height .+ z_0min
    @test canopy.boundary_conditions.turbulent_flux_parameterization.z_0b ==
          toml_dict["canopy_z_0b_coeff"] .* capped_height .+ z_0min

end
