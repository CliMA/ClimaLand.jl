using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Thermodynamics
using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP
using Dates


for FT in (Float32, Float64)
    @testset "Surface fluxes and radiation for soil, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        earth_param_set = LP.LandParameters(toml_dict)

        soil_domains = [
            ClimaLand.Domains.Column(;
                zlim = FT.((-100.0, 0.0)),
                nelements = 10,
            ),
            ClimaLand.Domains.HybridBox(;
                xlim = FT.((-1.0, 0.0)),
                ylim = FT.((-1.0, 0.0)),
                zlim = FT.((-100.0, 0.0)),
                nelements = (2, 2, 10),
                periodic = (true, true),
            ),
        ]
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        vg_m = FT(1) - FT(1) / vg_n
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        S_c = hcm.S_c

        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        emissivity = FT(0.99)
        z_0m = FT(0.001)
        z_0b = z_0m
        # Radiation
        start_date = DateTime(2005)
        SW_d = (t) -> 500
        LW_d = (t) -> 5.67e-8 * 280.0^4.0
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(SW_d),
            TimeVaryingInput(LW_d),
            start_date,
        )
        # Atmos
        precip = (t) -> 1e-8
        precip_snow = (t) -> 0
        T_atmos = (t) -> 285
        u_atmos = (t) -> 3
        q_atmos = (t) -> 0.005
        h_atmos = FT(3)
        P_atmos = (t) -> 101325
        atmos = PrescribedAtmosphere(
            TimeVaryingInput(precip),
            TimeVaryingInput(precip_snow),
            TimeVaryingInput(T_atmos),
            TimeVaryingInput(u_atmos),
            TimeVaryingInput(q_atmos),
            TimeVaryingInput(P_atmos),
            start_date,
            h_atmos,
            toml_dict,
        )
        @test atmos.gustiness == FT(1)
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

        for domain in soil_domains
            NIR_albedo_dry = fill(FT(0.4), domain.space.surface)
            PAR_albedo_dry = fill(FT(0.2), domain.space.surface)
            NIR_albedo_wet = fill(FT(0.3), domain.space.surface)
            PAR_albedo_wet = fill(FT(0.1), domain.space.surface)
            albedo = ClimaLand.Soil.CLMTwoBandSoilAlbedo{FT}(;
                NIR_albedo_dry,
                NIR_albedo_wet,
                PAR_albedo_dry,
                PAR_albedo_wet,
            )
            params = ClimaLand.Soil.EnergyHydrologyParameters(
                toml_dict;
                ν,
                ν_ss_om,
                ν_ss_quartz,
                ν_ss_gravel,
                hydrology_cm = hcm,
                K_sat,
                S_s,
                θ_r,
                albedo,
                emissivity,
                z_0m,
                z_0b,
            )
            model = Soil.EnergyHydrology{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_fluxes,
                sources = (),
            )
            @test ClimaComms.context(model) == ClimaComms.context()
            @test ClimaComms.device(model) == ClimaComms.device()
            drivers = ClimaLand.get_drivers(model)
            @test drivers == (atmos, radiation)
            Y, p, coords = initialize(model)
            Δz_top = model.domain.fields.Δz_top
            @test propertynames(p.drivers) == (
                :P_liq,
                :P_snow,
                :T,
                :P,
                :u,
                :q,
                :c_co2,
                :SW_d,
                :LW_d,
                :cosθs,
                :frac_diff,
            )
            @test propertynames(p.soil.turbulent_fluxes) ==
                  (:lhf, :shf, :vapor_flux_liq, :vapor_flux_ice)
            @test propertynames(p.soil) == (
                :total_water,
                :total_energy,
                :K,
                :ψ,
                :θ_l,
                :T,
                :κ,
                :Tf_depressed,
                :bidiag_matrix_scratch,
                :full_bidiag_matrix_scratch,
                :turbulent_fluxes,
                :R_n,
                :top_bc,
                :top_bc_wvec,
                :sfc_scratch,
                :PAR_albedo,
                :NIR_albedo,
                :sub_sfc_scratch,
                :infiltration,
                :bottom_bc,
                :bottom_bc_wvec,
            )
            function init_soil!(Y, z, params)
                ν = params.ν
                FT = eltype(ν)
                Y.soil.ϑ_l .= ν / 2
                Y.soil.θ_i .= 0
                T = FT(280)
                ρc_s = Soil.volumetric_heat_capacity(
                    ν / 2,
                    FT(0),
                    params.ρc_ds,
                    params.earth_param_set,
                )
                Y.soil.ρe_int =
                    Soil.volumetric_internal_energy.(
                        FT(0),
                        ρc_s,
                        T,
                        params.earth_param_set,
                    )
            end

            t = Float64(0)
            init_soil!(Y, coords.subsurface.z, model.parameters)
            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, t)
            space = axes(p.drivers.P_liq)
            @test p.drivers.P_liq == zeros(space) .+ FT(1e-8)
            @test p.drivers.P_snow == zeros(space) .+ FT(0)
            @test p.drivers.T == zeros(space) .+ FT(285)
            @test p.drivers.u == zeros(space) .+ FT(3)
            @test p.drivers.q == zeros(space) .+ FT(0.005)
            @test p.drivers.P == zeros(space) .+ FT(101325)
            @test p.drivers.LW_d == zeros(space) .+ FT(5.67e-8 * 280.0^4.0)
            @test p.drivers.SW_d == zeros(space) .+ FT(500)
            face_space =
                ClimaCore.Spaces.face_space(model.domain.space.subsurface)
            N = ClimaCore.Spaces.nlevels(face_space)
            surface_space = model.domain.space.surface
            z_sfc = ClimaCore.Fields.Field(
                ClimaCore.Fields.field_values(
                    ClimaCore.Fields.level(
                        ClimaCore.Fields.coordinate_field(face_space).z,
                        ClimaCore.Utilities.PlusHalf(N - 1),
                    ),
                ),
                surface_space,
            )
            T_sfc = ClimaCore.Fields.zeros(surface_space) .+ FT(280.0)
            @test ClimaLand.surface_emissivity(model, Y, p) == emissivity
            PAR_albedo = p.soil.PAR_albedo
            NIR_albedo = p.soil.NIR_albedo
            @test Base.Broadcast.materialize(
                ClimaLand.surface_albedo(model, Y, p),
            ) == PAR_albedo ./ 2 .+ NIR_albedo ./ 2
            @test ClimaLand.component_temperature(model, Y, p) == T_sfc

            conditions = copy(p.soil.turbulent_fluxes)
            ClimaLand.turbulent_fluxes!(
                conditions,
                model.boundary_conditions.top.atmos,
                model,
                Y,
                p,
                t,
            )
            R_n_copy = copy(p.soil.R_n)
            ClimaLand.net_radiation!(
                R_n_copy,
                model.boundary_conditions.top.radiation,
                model,
                Y,
                p,
                t,
            )
            @test R_n_copy == p.soil.R_n
            @test conditions == p.soil.turbulent_fluxes

            ClimaLand.Soil.soil_boundary_fluxes!(
                top_bc,
                ClimaLand.TopBoundary(),
                model,
                nothing,
                Y,
                p,
                t,
            )
            computed_water_flux = p.soil.top_bc.water
            computed_energy_flux = p.soil.top_bc.heat

            expected_water_flux = @. FT(precip(t)) .+ conditions.vapor_flux_liq
            @test computed_water_flux == expected_water_flux
            expected_energy_flux = @. R_n_copy +
               conditions.lhf +
               conditions.shf +
               FT(precip(t)) * Soil.volumetric_internal_energy_liq(
                   FT(T_atmos(t)),
                   earth_param_set,
               )
            @test computed_energy_flux == expected_energy_flux

            # Test soil resistances for liquid water
            θ_sfc = range(θ_r + eps(FT), ν, 5)
            S_sfc = @. ClimaLand.Soil.effective_saturation(ν, θ_sfc, θ_r)
            (; d_ds, evap_p, evap_α, hydrology_cm, ν, θ_r) = params
            S_c = hydrology_cm.S_c * evap_α
            dsl = @. ClimaLand.Soil.dry_soil_layer_thickness(
                S_sfc,
                S_c,
                d_ds,
                evap_p,
            )
            @test extrema(dsl ./ d_ds)[1] == 0
            @test extrema(dsl ./ d_ds)[2] ==
                  ((S_c - S_sfc[1]) / S_sfc[1])^evap_p
            _D_vapor = FT(LP.D_vapor(earth_param_set))
            S_sfc = @. ClimaLand.Soil.effective_saturation(ν, θ_sfc, θ_r)
            gsoil = @. ClimaLand.Soil.soil_conductance(
                S_sfc,
                hydrology_cm.S_c,
                d_ds,
                evap_p,
                evap_α,
                _D_vapor,
            )
            @test gsoil[1] < gsoil[2]
        end
    end
end

# Test CompositionBasedSoilAlbedo
for FT in (Float32, Float64)
    @testset "CompositionBasedSoilAlbedo, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)

        domain = ClimaLand.Domains.Column(;
            zlim = FT.((-100.0, 0.0)),
            nelements = 10,
        )

        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100)
        S_s = FT(1e-3)
        vg_n = FT(2.0)
        vg_α = FT(2.6)
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)

        # Test with different soil compositions
        # Case 1: Sandy desert (high vg_n, low OM) - expect high albedo
        ν_ss_om_desert = FT(0.0)
        ν_ss_quartz_desert = FT(0.7)
        ν_ss_gravel_desert = FT(0.3)

        # Case 2: Organic soil (low vg_n, high OM)
        ν_ss_om_organic = FT(0.1)
        ν_ss_quartz_organic = FT(0.2)
        ν_ss_gravel_organic = FT(0.1)

        emissivity = FT(0.99)
        z_0m = FT(0.001)
        z_0b = z_0m

        # Create composition-based albedo with default coefficients
        comp_albedo = ClimaLand.Soil.CompositionBasedSoilAlbedo{FT}()

        # Test that coefficients have physically correct signs:
        # - Organic matter darkens soil (negative coefficient)
        # - Higher vg_n means sandier = brighter (positive coefficient)
        # - Coarse fragments (gravel/rocks) are typically bright (positive coefficient)
        @test comp_albedo.c_om_PAR < 0  # OM darkens
        @test comp_albedo.c_vgn_PAR > 0  # sandier = brighter
        @test comp_albedo.c_cf_PAR > 0  # rocks are bright
        @test comp_albedo.c_om_NIR < 0
        @test comp_albedo.c_vgn_NIR > 0
        @test comp_albedo.c_cf_NIR > 0

        # Test that albedo bounds are physically reasonable
        @test comp_albedo.α_min >= FT(0)
        @test comp_albedo.α_max <= FT(1)
        @test comp_albedo.α_min < comp_albedo.α_max

        # Test that moisture parameters are physically reasonable
        @test FT(0) < comp_albedo.f_wet <= FT(1)  # moisture darkening factor
        @test comp_albedo.β > FT(0)  # positive exponent

        # Test with custom coefficients (verify they are applied)
        custom_albedo = ClimaLand.Soil.CompositionBasedSoilAlbedo{FT}(;
            η₀_PAR = FT(-4.0),
            c_vgn_PAR = FT(2.0),
            f_wet = FT(0.4),
        )
        @test custom_albedo.η₀_PAR == FT(-4.0)
        @test custom_albedo.c_vgn_PAR == FT(2.0)
        @test custom_albedo.f_wet == FT(0.4)

        # Test dry_albedo_from_composition: logistic output is always bounded
        α_min = comp_albedo.α_min
        α_max = comp_albedo.α_max

        # Even with extreme inputs, output must be in [α_min, α_max]
        extreme_high = ClimaLand.Soil.dry_albedo_from_composition(
            FT(10.0),   # very high intercept
            FT(0.0),
            FT(100.0),  # extreme positive vg_n effect
            FT(100.0),
            FT(0.0),    # no OM
            FT(3.0),    # high vg_n
            FT(1.0),    # high CF
            α_min,
            α_max,
        )
        @test extreme_high <= α_max
        @test extreme_high >= α_min

        extreme_low = ClimaLand.Soil.dry_albedo_from_composition(
            FT(-10.0),   # very low intercept
            FT(-100.0),  # extreme OM darkening
            FT(0.0),
            FT(0.0),
            FT(1.0),     # high OM
            FT(1.0),     # low vg_n
            FT(0.0),     # no CF
            α_min,
            α_max,
        )
        @test extreme_low <= α_max
        @test extreme_low >= α_min

        # Test physical relationships using actual default coefficients
        # Higher vg_n (sandier texture) → brighter albedo
        albedo_low_n = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.02),  # low OM
            FT(1.5),   # lower vg_n (clay-like)
            FT(0.2),
            α_min,
            α_max,
        )
        albedo_high_n = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.02),  # same low OM
            FT(2.5),   # higher vg_n (sandy)
            FT(0.2),   # same CF
            α_min,
            α_max,
        )
        @test albedo_high_n > albedo_low_n  # sandier soil should be brighter

        # More coarse fragments → brighter albedo
        albedo_low_cf = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.02),
            FT(2.0),
            FT(0.1),  # low gravel
            α_min,
            α_max,
        )
        albedo_high_cf = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.02),
            FT(2.0),
            FT(0.5),  # high gravel
            α_min,
            α_max,
        )
        @test albedo_high_cf > albedo_low_cf  # rocky terrain should be brighter

        # More organic matter → darker albedo
        albedo_low_om = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.01),  # low OM
            FT(2.0),
            FT(0.2),
            α_min,
            α_max,
        )
        albedo_high_om = ClimaLand.Soil.dry_albedo_from_composition(
            comp_albedo.η₀_PAR,
            comp_albedo.c_om_PAR,
            comp_albedo.c_vgn_PAR,
            comp_albedo.c_cf_PAR,
            FT(0.15),  # high OM
            FT(2.0),
            FT(0.2),
            α_min,
            α_max,
        )
        @test albedo_high_om < albedo_low_om  # organic soil should be darker

        # Test albedo_with_nonlinear_moisture helper: moisture always darkens
        α_dry = FT(0.3)
        f_wet = FT(0.5)
        β = FT(0.5)

        # Dry soil (S_e = 0) should have full dry albedo
        α_at_dry = ClimaLand.Soil.albedo_with_nonlinear_moisture(α_dry, FT(0), f_wet, β)
        @test α_at_dry ≈ α_dry

        # Wet soil should be darker than dry
        α_at_wet = ClimaLand.Soil.albedo_with_nonlinear_moisture(α_dry, FT(1), f_wet, β)
        @test α_at_wet < α_dry
        @test α_at_wet ≈ α_dry * (FT(1) - f_wet)  # at full saturation

        # Intermediate moisture should be between dry and wet
        α_at_mid = ClimaLand.Soil.albedo_with_nonlinear_moisture(α_dry, FT(0.5), f_wet, β)
        @test α_at_mid < α_dry
        @test α_at_mid > α_at_wet

        # Test full model integration with sandy desert soil
        # Radiation
        start_date = DateTime(2005)
        SW_d = (t) -> 500
        LW_d = (t) -> 5.67e-8 * 280.0^4.0
        radiation = PrescribedRadiativeFluxes(
            FT,
            TimeVaryingInput(SW_d),
            TimeVaryingInput(LW_d),
            start_date,
        )

        # Atmos
        precip = (t) -> 1e-8
        precip_snow = (t) -> 0
        T_atmos = (t) -> 285
        u_atmos = (t) -> 3
        q_atmos = (t) -> 0.005
        h_atmos = FT(3)
        P_atmos = (t) -> 101325
        atmos = PrescribedAtmosphere(
            TimeVaryingInput(precip),
            TimeVaryingInput(precip_snow),
            TimeVaryingInput(T_atmos),
            TimeVaryingInput(u_atmos),
            TimeVaryingInput(q_atmos),
            TimeVaryingInput(P_atmos),
            start_date,
            h_atmos,
            toml_dict,
        )

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

        params = ClimaLand.Soil.EnergyHydrologyParameters(
            toml_dict;
            ν,
            ν_ss_om = ν_ss_om_desert,
            ν_ss_quartz = ν_ss_quartz_desert,
            ν_ss_gravel = ν_ss_gravel_desert,
            hydrology_cm = hcm,
            K_sat,
            S_s,
            θ_r,
            albedo = comp_albedo,
            emissivity,
            z_0m,
            z_0b,
        )

        model = Soil.EnergyHydrology{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_fluxes,
            sources = (),
        )

        Y, p, coords = initialize(model)

        # Initialize with dry conditions
        Y.soil.ϑ_l .= θ_r + FT(0.01)
        Y.soil.θ_i .= 0
        T = FT(280)
        ρc_s = Soil.volumetric_heat_capacity(
            θ_r + FT(0.01),
            FT(0),
            params.ρc_ds,
            params.earth_param_set,
        )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                FT(0),
                ρc_s,
                T,
                params.earth_param_set,
            )

        t = Float64(0)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t)

        # Check that PAR and NIR are different (composition model has different coefficients)
        par_vals = Array(parent(p.soil.PAR_albedo))
        nir_vals = Array(parent(p.soil.NIR_albedo))
        @test !all(par_vals .≈ nir_vals)

        # Check NIR > PAR (typical for mineral soils - NIR reflects more than visible)
        @test all(nir_vals .> par_vals)

        # For dry soil, albedo should be within physical bounds
        @test all(par_vals .> comp_albedo.α_min)
        @test all(par_vals .< comp_albedo.α_max)
        @test all(nir_vals .> comp_albedo.α_min)
        @test all(nir_vals .< comp_albedo.α_max)

        # Test with wet conditions
        Y.soil.ϑ_l .= ν - FT(0.01)
        set_initial_cache!(p, Y, t)

        par_vals_wet = Array(parent(p.soil.PAR_albedo))
        nir_vals_wet = Array(parent(p.soil.NIR_albedo))

        # Wet albedo should be lower than dry (moisture darkens soil)
        @test all(par_vals_wet .< par_vals)
        @test all(nir_vals_wet .< nir_vals)

        # Check surface albedo function
        α_sfc =
            Base.Broadcast.materialize(ClimaLand.surface_albedo(model, Y, p))
        expected_mean = (p.soil.PAR_albedo .+ p.soil.NIR_albedo) ./ 2
        @test Array(parent(α_sfc)) ≈ Array(parent(expected_mean))
    end
end
