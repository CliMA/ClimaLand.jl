using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
using Dates
using ClimaParams
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Default constructors, FT = $FT" begin
        domain = ClimaLand.Domains.global_domain(FT)
        atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT)
        forcing = (; atmos, radiation)
        toml_dict = ClimaLand.Parameters.create_toml_dict(
            FT,
            joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml"),
        )
        prognostic_land_components = (:canopy, :soil, :soilco2)

        # Soil model
        soil = EnergyHydrology{FT}(
            domain,
            forcing,
            toml_dict;
            prognostic_land_components,
            additional_sources = (ClimaLand.RootExtraction{FT}(),),
        )

        # SoilCO2 model
        co2_prognostic_soil =
            Soil.Biogeochemistry.PrognosticMet(soil.parameters)
        soil_organic_carbon = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput((t) -> 5),
        )
        soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(
            co2_prognostic_soil,
            soil_organic_carbon,
            atmos,
        )
        soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(domain, soilco2_drivers)

        # Canopy model
        canopy_domain = ClimaLand.Domains.obtain_surface_domain(domain)
        LAI = TimeVaryingInput((t) -> FT(1.0))
        ground = ClimaLand.PrognosticGroundConditions{FT}()
        canopy_forcing = (; atmos, radiation, ground)
        canopy = Canopy.CanopyModel{FT}(
            canopy_domain,
            canopy_forcing,
            LAI,
            toml_dict;
            prognostic_land_components,
        )

        model = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain)
        # The constructor has many asserts that check the model
        # components, so we don't need to check them again here.
        Y, p, cds = initialize(model)
        # check that albedos have been added to cache
        @test haskey(p.soil, :PAR_albedo)
        @test haskey(p.soil, :NIR_albedo)
        # initialize cache, then check that albedos are set to the correct values
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, 0.0)
        canopy_bc = model.canopy.boundary_conditions
        α_soil_PAR = Canopy.ground_albedo_PAR(
            Val(canopy_bc.prognostic_land_components),
            canopy_bc.ground,
            Y,
            p,
            0.0,
        )
        @test p.soil.PAR_albedo == α_soil_PAR
        α_soil_NIR = Canopy.ground_albedo_NIR(
            Val(canopy_bc.prognostic_land_components),
            canopy_bc.ground,
            Y,
            p,
            0.0,
        )
        @test p.soil.NIR_albedo == α_soil_NIR
    end
end
