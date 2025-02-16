using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
using Dates
import ClimaLand.Parameters as LP
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy.PlantHydraulics
for FT in (Float32, Float64)
    @testset "PrognosticSoil, FT = $FT" begin
        # Discretize radiation into 2 bands
        spectral_discretization = ClimaLand.TwoBandSpectralDiscretization{FT}()
        # Only care about PAR and NIR albedo values
        albedo = (0.2, 0.4)
        # setup SoilCanopyModel with dummy params
        earth_param_set = LP.LandParameters(FT)
        radius = FT(6378.1e3)
        depth = FT(50)
        nelements = (50, 10)
        domain = ClimaLand.Domains.SphericalShell(;
            radius = radius,
            depth = depth,
            nelements = nelements,
            npolynomial = 1,
            dz_tuple = FT.((10.0, 0.1)),
        )
        hcm = vanGenuchten{FT}(; α = FT(0), n = FT(0))
        # Radiation
        start_date = DateTime(2005)
        atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT)
        top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(atmos, radiation)
        zero_water_flux = WaterFluxBC((p, t) -> 0.0)
        zero_heat_flux = HeatFluxBC((p, t) -> 0.0)
        boundary_fluxes = (;
            top = top_bc,
            bottom = WaterHeatBC(;
                water = zero_water_flux,
                heat = zero_heat_flux,
            ),
        )
        soil_params = ClimaLand.Soil.EnergyHydrologyParameters(
            FT;
            ν = FT(0.1),
            ν_ss_om = FT(0.01),
            ν_ss_quartz = FT(0.01),
            ν_ss_gravel = FT(0.01),
            hydrology_cm = hcm,
            K_sat = FT(0),
            S_s = FT(0),
            θ_r = FT(0),
            albedo,
            emissivity = FT(0),
            z_0m = FT(0),
            z_0b = FT(0),
        )
        soil_args = (domain = domain, parameters = soil_params)
        soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
        soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
        soilco2_sources = ()
        soilco2_boundary_conditions =
            (; top = soilco2_top_bc, bottom = soilco2_bot_bc)
        soilco2_ps = SoilCO2ModelParameters(FT; D_ref = FT(0.0))
        soilco2_args = (;
            boundary_conditions = soilco2_boundary_conditions,
            sources = soilco2_sources,
            domain = domain,
            parameters = soilco2_ps,
        )
        canopy_component_types = (;
            autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
            radiative_transfer = Canopy.TwoStreamModel{FT},
            photosynthesis = Canopy.FarquharModel{FT},
            conductance = Canopy.MedlynConductanceModel{FT},
            hydraulics = Canopy.PlantHydraulicsModel{FT},
            energy = Canopy.BigLeafEnergyModel{FT},
        )
        autotrophic_respiration_args =
            (; parameters = Canopy.AutotrophicRespirationParameters(FT))
        radiative_transfer_args = (;
            parameters = Canopy.TwoStreamParameters(
                FT;
                spectral_discretization = spectral_discretization,
                Ω = FT(0),
                ρ_leaf = (FT(0), FT(0)),
                τ_leaf = (FT(0), FT(0)),
            )
        )
        conductance_args =
            (; parameters = Canopy.MedlynConductanceParameters(FT; g1 = FT(0)))
        photosynthesis_args = (;
            parameters = Canopy.FarquharParameters(FT, FT(0); Vcmax25 = FT(0))
        )
        ai_parameterization = PlantHydraulics.PrescribedSiteAreaIndex{FT}(
            TimeVaryingInput(identity),
            FT(0),
            FT(0),
        )

        retention_model = PlantHydraulics.LinearRetentionCurve{FT}(FT(0))
        conductivity_model = PlantHydraulics.Weibull{FT}(FT(0), FT(0), FT(0))
        plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
            ai_parameterization = ai_parameterization,
            ν = FT(0),
            S_s = FT(0),
            rooting_depth = FT(0),
            conductivity_model = conductivity_model,
            retention_model = retention_model,
        )
        plant_hydraulics_args = (
            parameters = plant_hydraulics_ps,
            n_stem = 0,
            n_leaf = 1,
            compartment_midpoints = [FT(0.0)],
            compartment_surfaces = [FT(0.0), FT(0.0)],
        )
        energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(FT(0)),)
        canopy_component_args = (;
            autotrophic_respiration = autotrophic_respiration_args,
            radiative_transfer = radiative_transfer_args,
            photosynthesis = photosynthesis_args,
            conductance = conductance_args,
            hydraulics = plant_hydraulics_args,
            energy = energy_args,
        )
        shared_params =
            Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
                FT(0),
                FT(0),
                earth_param_set,
            )
        canopy_model_args = (;
            parameters = shared_params,
            domain = ClimaLand.obtain_surface_domain(domain),
        )
        runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
            f_over = FT(0),
            f_max = fill(FT(0), domain.space.surface),
            R_sb = FT(0),
        )
        Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(
            TimeVaryingInput(identity),
        )
        land_input = (
            atmos = atmos,
            radiation = radiation,
            runoff = runoff_model,
            soil_organic_carbon = Csom,
        )
        land_input = (
            atmos = atmos,
            radiation = radiation,
            runoff = runoff_model,
            soil_organic_carbon = Csom,
        )
        model = SoilCanopyModel{FT}(;
            soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT},
            soilco2_args = soilco2_args,
            land_args = land_input,
            soil_model_type = Soil.EnergyHydrology{FT},
            soil_args = soil_args,
            canopy_component_types = canopy_component_types,
            canopy_component_args = canopy_component_args,
            canopy_model_args = canopy_model_args,
        )

        # initialize model
        Y, p, coords = initialize(model)
        # check that albedos have been added to cache
        @test haskey(p.soil, :albedo)
        # initialize cache, then check that albedos are set to the correct values
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, 0.0)
        canopy_bc = model.canopy.boundary_conditions
        @test all(
            Array(
                parent(
                    Canopy.ground_albedo(
                        Val(canopy_bc.prognostic_land_components),
                        canopy_bc.ground,
                        Y,
                        p,
                        0.0,
                    ),
                ),
            ) .== ground_albedo,
        )
    end
end
