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
    toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
    @testset "Default constructors, FT = $FT" begin
        domain = ClimaLand.Domains.global_domain(FT)
        atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
        forcing = (; atmos, radiation)
        toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
        LAI = TimeVaryingInput((t) -> FT(1.0))
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
    @testset "Hydraulic options, FT = $FT" begin
        domain = ClimaLand.Domains.Column(;
            zlim = (FT(-1), FT(0)),
            nelements = 10,
            longlat = (FT(-110), FT(34)),
        )
        toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)
        surface_space = domain.space.surface
        start_date = DateTime(2008)
        era5_ncdata_path =
            ClimaLand.Artifacts.era5_land_forcing_data2008_path(; lowres = true)
        forcing = ClimaLand.prescribed_forcing_era5(
            era5_ncdata_path,
            surface_space,
            start_date,
            earth_param_set,
            FT;
        )
        ground = ClimaLand.PrognosticGroundConditions{FT}()
        LAI = ClimaLand.Canopy.prescribed_lai_modis(
            surface_space,
            start_date,
            start_date,
        )
        surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)

        land = SoilCanopyModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain;
            canopy = ClimaLand.CanopyModel{FT}(
                surface_domain,
                (; forcing.atmos, forcing.radiation, ground),
                LAI,
                toml_dict,
                soil_moisture_stress = ClimaLand.Canopy.PiecewiseMoistureStressModel{
                    FT,
                }(
                    domain,
                    toml_dict,
                ),
                prognostic_land_components = (:canopy, :soil, :soilco2),
            ),
        )
        Y, p, _ = initialize(land)

        # Create the model with a steady state hydraulics
        land_ss = SoilCanopyModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain;
            canopy = ClimaLand.CanopyModel{FT}(
                surface_domain,
                (; forcing.atmos, forcing.radiation, ground),
                LAI,
                toml_dict,
                hydraulics = ClimaLand.Canopy.PlantHydraulics.SteadyStateModel{
                    FT,
                }(),
                soil_moisture_stress = ClimaLand.Canopy.PiecewiseMoistureStressModel{
                    FT,
                }(
                    domain,
                    toml_dict,
                ),
                prognostic_land_components = (:canopy, :soil, :soilco2),
            ),
        )
        Y_ss, p_ss, _ = initialize(land_ss)
        set_initial_cache! = make_set_initial_cache(land)
        set_initial_cache_ss! = make_set_initial_cache(land_ss)
        # Soil IC
        ϑ_l0 = land.soil.parameters.ν ./ 2
        θ_i0 = land.soil.parameters.ν ./ 5
        T = FT(270.0)
        ρc_s = @. ClimaLand.Soil.volumetric_heat_capacity(
            ϑ_l0,
            θ_i0,
            land.soil.parameters.ρc_ds,
            earth_param_set,
        )
        ρe_int0 = @. ClimaLand.Soil.volumetric_internal_energy(
            θ_i0,
            ρc_s,
            T,
            earth_param_set,
        )
        Y.soil.ϑ_l .= ϑ_l0
        Y_ss.soil.ϑ_l .= ϑ_l0
        Y.soil.θ_i .= θ_i0
        Y_ss.soil.θ_i .= θ_i0
        Y.soil.ρe_int = ρe_int0
        Y_ss.soil.ρe_int = ρe_int0

        # Canopy IC
        ϑ0 = land.canopy.hydraulics.parameters.ν / 2
        CTemp0 = FT(290.5)
        Y.canopy.hydraulics.ϑ_l.:1 .= ϑ0
        Y.canopy.energy.T .= CTemp0
        Y_ss.canopy.energy.T .= CTemp0

        set_initial_cache!(p, Y, 0.0)
        set_initial_cache_ss!(p_ss, Y_ss, 0.0)

        @test p_ss.canopy.turbulent_fluxes.transpiration ==
              p.canopy.turbulent_fluxes.transpiration
        # Check that root extraction seen by soil matches transpiration
        sfc_field = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaCore.Operators.column_integral_definite!(
            sfc_field,
            p_ss.root_extraction,
        )
        @test sfc_field == p_ss.canopy.turbulent_fluxes.transpiration
    end
end
