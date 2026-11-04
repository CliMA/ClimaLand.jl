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
    @testset "Steady state plant hydraulics with prognostic soil, FT = $FT" begin
        domain = ClimaLand.Domains.Column(;
            zlim = (FT(-2), FT(0)),
            nelements = 10,
            longlat = (FT(-92.2), FT(38.7)),
        )
        LAI = TimeVaryingInput((t) -> 1.0)
        # Near local solar noon in July, so that the canopy transpires
        start_date = DateTime(2005, 7, 1, 18)
        stop_date = start_date + Day(1)
        forcing = prescribed_forcing_era5(
            start_date,
            stop_date,
            domain.space.surface,
            toml_dict,
            FT;
            use_lowres_forcing = true,
        )
        t0 = ClimaUtilities.TimeManager.ITime(0, Second(1), start_date)
        z = domain.fields.z
        function initialize_state(model; ϑ_l = nothing)
            Y, p, _ = initialize(model)
            set_ic! =
                ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
                    model,
                )
            set_ic!(Y, p, t0, model)
            isnothing(ϑ_l) || (Y.soil.ϑ_l .= ϑ_l)
            make_set_initial_cache(model)(p, Y, t0)
            return Y, p
        end
        steady_state_canopy(soil_moisture_stress) = Canopy.CanopyModel{FT}(
            ClimaLand.Domains.obtain_surface_domain(domain),
            (;
                atmos = forcing.atmos,
                radiation = forcing.radiation,
                ground = ClimaLand.PrognosticGroundConditions{FT}(),
            ),
            LAI,
            toml_dict;
            prognostic_land_components = (:canopy, :soil, :soilco2),
            soil_moisture_stress,
            hydraulics = Canopy.SteadyStateModel{FT}(),
        )
        transpiration(p) = Array(parent(p.canopy.turbulent_fluxes.vapor_flux))
        fa_roots(p) = Array(parent(p.canopy.hydraulics.fa_roots))

        model = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain)
        @test model.canopy.hydraulics isa Canopy.SteadyStateModel
        (; θ_high, θ_low, c, ψ_wilt) = model.canopy.soil_moisture_stress
        Y, p = initialize_state(model)
        ϑ_wet = copy(Y.soil.ϑ_l)
        T = transpiration(p)
        @test all(T .> 1e-9)
        # Root extraction balances transpiration, so the plant water is constant
        @test fa_roots(p) ≈ T rtol = 100 * eps(FT)
        dY = similar(Y)
        ClimaLand.make_compute_exp_tendency(
            model.canopy.hydraulics,
            model.canopy,
        )(
            dY,
            Y,
            p,
            t0,
        )
        @test all(abs.(parent(dY.canopy.hydraulics.ϑ_l)) .<= 100 * eps(FT) .* T)

        # The root zone stress factor is the root-weighted layer stress factor
        βm = @. Canopy.compute_piecewise_moisture_stress(
            θ_high,
            θ_low,
            c,
            p.soil.θ_l,
        )
        norm = p.scratch1 # A surface scratch field
        root_weight =
            @. Canopy.root_distribution(z, model.canopy.biomass.rooting_depth)
        ClimaCore.Operators.column_integral_definite!(norm, root_weight)
        tmp = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaCore.Operators.column_integral_definite!(
            tmp,
            βm .* root_weight ./ norm,
        )
        @test all(parent(@. p.canopy.soil_moisture_stress.βm - tmp .≈ 0))

        # Layers at the wilting point give no water to the canopy, and the
        # leaf water potential is a weighted potential of the wetter layers
        dry = Array(parent(z)) .> -0.4
        Y, p = initialize_state(model; ϑ_l = @. ifelse(z > -0.4, θ_low, ϑ_wet))
        root_extraction = Array(parent(p.root_extraction))
        @test all(root_extraction[dry] .== 0)
        @test all(root_extraction[.!dry] .> 0)
        @test fa_roots(p) ≈ transpiration(p) rtol = 100 * eps(FT)
        ψ_wet = Array(parent(p.soil.ψ))[.!dry]
        ψ_leaf = only(Array(parent(p.canopy.hydraulics.ψ)))
        @test minimum(ψ_wet) * (1 + sqrt(eps(FT))) <=
              ψ_leaf <=
              maximum(ψ_wet) * (1 - sqrt(eps(FT)))

        # At the wilting point everywhere, no water is extracted and the leaf
        # is at the wilting potential
        Y, p = initialize_state(model; ϑ_l = θ_low)
        @test all(parent(p.root_extraction) .== 0)
        @test all(fa_roots(p) .== 0)
        @test all(parent(p.canopy.hydraulics.ψ) .== ψ_wilt)

        # A non-default wilting potential, with spatially varying soil
        # parameters, gives consistent thresholds
        soil_params = (;
            ν = model.soil.parameters.ν,
            θ_r = model.soil.parameters.θ_r,
            hydrology_cm = model.soil.parameters.hydrology_cm,
        )
        sms_wet_wilt = Canopy.PiecewiseMoistureStressModel{FT}(
            domain,
            toml_dict;
            ψ_wilt = FT(-100),
            soil_params,
        )
        @test all(parent(sms_wet_wilt.θ_low) .> parent(θ_low))
        @test SoilCanopyModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain;
            soil = model.soil,
            canopy = steady_state_canopy(sms_wet_wilt),
        ) isa SoilCanopyModel
        # Thresholds not computed from the soil parameters are rejected
        sms_no_cm = Canopy.PiecewiseMoistureStressModel{FT}(
            domain,
            toml_dict;
            soil_params = (; soil_params.ν, soil_params.θ_r),
        )
        @test_throws AssertionError SoilCanopyModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain;
            soil = model.soil,
            canopy = steady_state_canopy(sms_no_cm),
        )

        # Moisture stress independent of soil moisture: layers at the residual
        # water content give no water to the canopy
        model_no_stress = SoilCanopyModel{FT}(
            forcing,
            LAI,
            toml_dict,
            domain;
            soil = model.soil,
            canopy = steady_state_canopy(Canopy.NoMoistureStressModel{FT}()),
        )
        θ_r = model.soil.parameters.θ_r
        Y, p = initialize_state(
            model_no_stress;
            ϑ_l = @. ifelse(z > -0.4, θ_r, ϑ_wet)
        )
        root_extraction = Array(parent(p.root_extraction))
        @test all(root_extraction[dry] .== 0)
        @test all(root_extraction[.!dry] .> 0)
        @test fa_roots(p) ≈ transpiration(p) rtol = 100 * eps(FT)
        ψ_leaf = only(Array(parent(p.canopy.hydraulics.ψ)))
        @test isfinite(ψ_leaf)
        @test ψ_leaf <= maximum(Array(parent(p.soil.ψ)))
    end

    @testset "Plant water initial conditions from the soil, FT = $FT" begin
        domain = ClimaLand.Domains.Column(;
            zlim = (FT(-2), FT(0)),
            nelements = 10,
            longlat = (FT(-92.2), FT(38.7)),
        )
        forcing = ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
        LAI = TimeVaryingInput((t) -> FT(1.0))
        model = SoilCanopyModel{FT}(
            (; atmos = forcing[1], radiation = forcing[2]),
            LAI,
            toml_dict,
            domain,
        )
        ψ_soil = ClimaCore.Fields.zeros(domain.space.subsurface) .- FT(50)
        ψ_roots =
            ClimaLand.Simulations.root_weighted_soil_potential(ψ_soil, model)
        @test all(parent(ψ_roots) .≈ FT(-50))

        hydraulics = Canopy.PlantHydraulicsModel{FT}(
            ClimaLand.Domains.obtain_surface_domain(domain),
            toml_dict,
        )
        (; ν, retention_model) = hydraulics.parameters
        ϑ_l = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaLand.Simulations.set_plant_water_from_potential!(
            ϑ_l,
            hydraulics,
            ψ_roots,
        )
        @test all(parent(ϑ_l) .≈ ν * (1 + retention_model.a * FT(-50)))
        # Below the potential at zero saturation, the water content is floored
        ψ_very_dry = ClimaCore.Fields.zeros(domain.space.surface) .- FT(1e4)
        ClimaLand.Simulations.set_plant_water_from_potential!(
            ϑ_l,
            hydraulics,
            ψ_very_dry,
        )
        @test all(parent(ϑ_l) .== eps(FT))
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
