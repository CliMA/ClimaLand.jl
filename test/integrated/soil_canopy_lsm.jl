using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
using Dates
using ClimaParams
import SurfaceFluxes
import ClimaLand.Parameters as LP
using ClimaUtilities
for FT in (Float32, Float64)
    toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
    @testset "Default constructors, FT = $FT" begin
        domain = ClimaLand.Domains.global_domain(FT)
        atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT; toml_dict) # h_atmos = 1 relative to the ground
        forcing = (; atmos, radiation)
        toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
        LAI = TimeVaryingInput((t) -> FT(1.0))
        Δt = FT(450)
        # Choose a canopy height such that h_atmos > displ
        biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(;
            LAI,
            SAI = FT(0),
            RAI = FT(1),
            rooting_depth = FT(1),
            height = FT(0.3),
        )
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            ClimaLand.Domains.obtain_surface_domain(domain),
            (;
                atmos = forcing.atmos,
                radiation = forcing.radiation,
                ground = ClimaLand.PrognosticGroundConditions{FT}(),
            ),
            LAI,
            toml_dict;
            prognostic_land_components = (:canopy, :soil, :soilco2),
            biomass,
        )
        model = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain; canopy)
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
    @testset "Root extraction with Steady State and Piecewise SMS, FT = $FT" begin
        domain = ClimaLand.Domains.Column(;
            zlim = (FT(-2), FT(0)),
            nelements = 10,
            longlat = (FT(-92.2), FT(38.7)),
        )
        LAI = TimeVaryingInput((t) -> 1.0)
        start_date = DateTime(2005)
        stop_date = start_date+Day(1)
        forcing = prescribed_forcing_era5(
            start_date,
            stop_date,
            domain.space.surface,
            toml_dict,
            FT;
            use_lowres_forcing = true,
        )
        model = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain)
        Y, p, cds = initialize(model)
        set_ic! =
            ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
                model,
            )
        t0 = ClimaUtilities.TimeManager.ITime(0, Second(1), start_date)
        set_ic!(Y, p, t0, model)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t0)
        z = domain.fields.z
        norm = p.scratch1 # A surface scratch field
        (; θ_high, θ_low, c) = model.canopy.soil_moisture_stress
        βm = @. Canopy.compute_piecewise_moisture_stress(
            θ_high,
            θ_low,
            c,
            p.soil.θ_l,
        )
        root_dryness_weight =
            @. Canopy.root_distribution(z, model.canopy.biomass.rooting_depth)*βm
        ClimaCore.Operators.column_integral_definite!(norm, root_dryness_weight)
        @test all(
            parent(
                @. p.root_extraction -
                   root_dryness_weight * p.canopy.turbulent_fluxes.vapor_flux /
                   (norm+sqrt(eps(FT)))
            ) .≈ 0,
        )

        root_weight =
            @. Canopy.root_distribution(z, model.canopy.biomass.rooting_depth)
        ClimaCore.Operators.column_integral_definite!(norm, root_weight)
        tmp = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaCore.Operators.column_integral_definite!(
            tmp,
            βm .* root_weight ./ norm,
        )
        @test all(parent(@. p.canopy.soil_moisture_stress.βm - tmp .≈ 0))
    end
end

@testset "Soil below a canopy is forced above the ground" begin
    FT = Float64
    toml_dict = ClimaLand.Parameters.create_toml_dict(FT)
    domain = ClimaLand.Domains.Column(;
        zlim = (FT(-2), FT(0)),
        nelements = 10,
        longlat = (FT(-92.2), FT(38.7)),
    )
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    start_date = DateTime(2008, 7, 1)
    # Forcing at a height above the apparent sink of the canopy, as for
    # reanalysis forcing
    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput((t) -> 0.0),
        TimeVaryingInput((t) -> 0.0),
        TimeVaryingInput((t) -> 298.0),
        TimeVaryingInput((t) -> 3.0),
        TimeVaryingInput((t) -> 0.012),
        TimeVaryingInput((t) -> 101325.0),
        start_date,
        FT(10),
        toml_dict;
        reference = SurfaceFluxes.ReferenceAboveApparentSink(),
    )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput((t) -> 600.0),
        TimeVaryingInput((t) -> 380.0),
        start_date;
        cosθs = (t, start_date) -> FT(0.8),
        frac_diff = TimeVaryingInput((t) -> 0.3),
        toml_dict,
    )
    forcing = (; atmos, radiation)
    LAI = TimeVaryingInput((t) -> FT(2))
    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI = FT(0.5),
        RAI = FT(1),
        rooting_depth = FT(1),
        height = FT(10),
    )
    canopy = Canopy.CanopyModel{FT}(
        surface_domain,
        (;
            atmos,
            radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components = (:canopy, :soil, :soilco2),
        biomass,
    )
    model = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain; canopy)
    Y, p, _ = initialize(model)
    earth_param_set = model.soil.parameters.earth_param_set
    Y.soil.ϑ_l .= FT(0.25)
    Y.soil.θ_i .= FT(0)
    ρc_s = Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        model.soil.parameters.ρc_ds,
        earth_param_set,
    )
    Y.soil.ρe_int .= Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        FT(295),
        earth_param_set,
    )
    Y.soilco2.CO2 .= FT(0.000412)
    Y.soilco2.O2 .= FT(0.21)
    Y.soilco2.SOC .= FT(5)
    Y.canopy.hydraulics.ϑ_l .= model.canopy.hydraulics.parameters.ν
    Y.canopy.energy.T .= FT(298)
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, 0.0)
    ClimaLand.Simulations.set_canopy_component_initial_conditions!(
        Y,
        p,
        model.canopy.photosynthesis,
        model.canopy,
    )
    set_initial_cache!(p, Y, 0.0)

    # The soil exchanges with the canopy air at the apparent sink of the
    # canopy, d + z_0m above the ground; the soil has no displacement height
    sfp = model.canopy.boundary_conditions.turbulent_flux_parameterization
    @test ClimaLand.ground_forcing(model, Y, p).reference ==
          SurfaceFluxes.ReferenceAboveSurface()
    @test all(parent(p.soil.turbulent_fluxes.Δz_eff) .≈ sfp.displ + sfp.z_0m)
    @test all(isfinite, parent(p.soil.turbulent_fluxes.shf))
    @test all(isfinite, parent(p.soil.turbulent_fluxes.lhf))
end
