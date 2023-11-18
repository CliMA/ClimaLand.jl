using Test
using ClimaCore
import CLIMAParameters as CP
using Thermodynamics
using ClimaLSM
using ClimaLSM.Soil
import ClimaLSM
import ClimaLSM.Parameters as LSMP
using Insolation
using Dates
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

for FT in (Float32, Float64)
    @testset "Surface fluxes and radiation for soil, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)

        soil_domains = [
            ClimaLSM.Domains.Column(;
                zlim = FT.((-100.0, 0.0)),
                nelements = 10,
            ),
            ClimaLSM.Domains.HybridBox(;
                xlim = FT.((-1.0, 0.0)),
                ylim = FT.((-1.0, 0.0)),
                zlim = FT.((-100.0, 0.0)),
                nelements = (2, 2, 10),
                npolynomial = 1,
                periodic = (true, true),
            ),
        ]
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        vg_m = FT(1) - FT(1) / vg_n
        hcm = vanGenuchten(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        S_c = hcm.S_c
        @test Soil.dry_soil_layer_thickness(FT(1), S_c, FT(1)) == FT(0)
        @test Soil.dry_soil_layer_thickness(FT(0), S_c, FT(1)) == FT(1)


        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        κ_minerals = FT(2.5)
        κ_om = FT(0.25)
        κ_quartz = FT(8.0)
        κ_air = FT(0.025)
        κ_ice = FT(2.21)
        κ_liq = FT(0.57)
        ρp = FT(2.66 / 1e3 * 1e6)
        ρc_ds = FT(2e6 * (1.0 - ν))
        κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
        κ_dry_soil = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
        κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
        κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
        emissivity = FT(0.99)
        PAR_albedo = FT(0.2)
        NIR_albedo = FT(0.4)
        z_0m = FT(0.001)
        z_0b = z_0m
        # Radiation
        ref_time = DateTime(2005)
        SW_d = (t) -> 500
        LW_d = (t) -> 5.67e-8 * 280.0^4.0
        radiation = PrescribedRadiativeFluxes(
            FT,
            SW_d,
            LW_d,
            ref_time;
            orbital_data = Insolation.OrbitalData(),
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
            precip,
            precip_snow,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            ref_time,
            h_atmos,
        )
        @test atmos.gustiness == FT(1)
        top_bc = ClimaLSM.Soil.AtmosDrivenFluxBC(atmos, radiation)
        zero_flux = FluxBC((p, t) -> 0.0)
        boundary_fluxes =
            (; top = top_bc, bottom = (water = zero_flux, heat = zero_flux))
        params = ClimaLSM.Soil.EnergyHydrologyParameters{FT}(;
            κ_dry = κ_dry_soil,
            κ_sat_frozen = κ_sat_frozen,
            κ_sat_unfrozen = κ_sat_unfrozen,
            ρc_ds = ρc_ds,
            ν = ν,
            ν_ss_om = ν_ss_om,
            ν_ss_quartz = ν_ss_quartz,
            ν_ss_gravel = ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
            PAR_albedo = PAR_albedo,
            NIR_albedo = NIR_albedo,
            emissivity = emissivity,
            z_0m = z_0m,
            z_0b = z_0b,
            earth_param_set = earth_param_set,
        )

        for domain in soil_domains
            model = Soil.EnergyHydrology{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_fluxes,
                sources = (),
            )

            Y, p, coords = initialize(model)
            function init_soil!(Y, z, params)
                ν = params.ν
                FT = eltype(ν)
                Y.soil.ϑ_l .= ν / 2
                Y.soil.θ_i .= 0
                T = FT(280)
                ρc_s = Soil.volumetric_heat_capacity(ν / 2, FT(0), params)
                Y.soil.ρe_int =
                    Soil.volumetric_internal_energy.(
                        FT(0),
                        ρc_s,
                        T,
                        Ref(params),
                    )
            end

            t = Float64(0)
            init_soil!(Y, coords.subsurface.z, model.parameters)
            set_initial_aux_state! = make_set_initial_aux_state(model)
            set_initial_aux_state!(p, Y, t)


            face_space = ClimaLSM.Domains.obtain_face_space(
                model.domain.space.subsurface,
            )
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
            @test ClimaLSM.surface_emissivity(model, Y, p) == emissivity
            @test ClimaLSM.surface_evaporative_scaling(model, Y, p) == FT(1)
            @test ClimaLSM.surface_height(model, Y, p) == z_sfc
            @test ClimaLSM.surface_albedo(model, Y, p) ==
                  PAR_albedo / 2 + NIR_albedo / 2
            @test ClimaLSM.surface_temperature(model, Y, p, t) == T_sfc

            thermo_params =
                LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
            ts_in = construct_atmos_ts(atmos, t, thermo_params)
            ρ_sfc = compute_ρ_sfc.(Ref(thermo_params), Ref(ts_in), T_sfc)
            @test ClimaLSM.surface_air_density(
                model.boundary_conditions.top.atmos,
                model,
                Y,
                p,
                t,
                T_sfc,
            ) == ρ_sfc

            q_sat =
                Thermodynamics.q_vap_saturation_generic.(
                    Ref(thermo_params),
                    T_sfc,
                    ρ_sfc,
                    Ref(Thermodynamics.Liquid()),
                )
            g = LSMP.grav(model.parameters.earth_param_set)
            M_w = LSMP.molar_mass_water(model.parameters.earth_param_set)
            R = LSMP.gas_constant(model.parameters.earth_param_set)
            ψ_sfc =
                parent(p.soil.ψ)[end] .+ ClimaCore.Fields.zeros(surface_space)
            q_sfc = @. (q_sat * exp(g * ψ_sfc * M_w / (R * T_sfc)))
            @test ClimaLSM.surface_specific_humidity(
                model,
                Y,
                p,
                T_sfc,
                ρ_sfc,
            ) == q_sfc

            conditions = ClimaLSM.surface_fluxes(
                model.boundary_conditions.top.atmos,
                model,
                Y,
                p,
                t,
            )
            R_n = ClimaLSM.net_radiation(
                model.boundary_conditions.top.radiation,
                model,
                Y,
                p,
                t,
            )

            (computed_water_flux, computed_energy_flux) =
                ClimaLSM.Soil.soil_boundary_fluxes(
                    top_bc,
                    ClimaLSM.TopBoundary(),
                    model,
                    nothing,
                    Y,
                    p,
                    t,
                )

            (; ν, θ_r, d_ds) = model.parameters
            _D_vapor = FT(LSMP.D_vapor(model.parameters.earth_param_set))
            S_l_sfc = ClimaLSM.Domains.top_center_to_surface(
                Soil.effective_saturation.(ν, Y.soil.ϑ_l, θ_r),
            )
            τ_a = ClimaLSM.Domains.top_center_to_surface(
                @. (ν - p.soil.θ_l - Y.soil.θ_i)^(FT(5 / 2)) / ν
            )
            dsl = Soil.dry_soil_layer_thickness.(S_l_sfc, S_c, d_ds)
            r_soil = @. dsl / (_D_vapor * τ_a) # [s\m]
            r_ae = conditions.r_ae
            expected_water_flux = @. FT(atmos.liquid_precip(t)) .+
               conditions.vapor_flux * r_ae / (r_soil + r_ae)
            @test computed_water_flux == expected_water_flux
            expected_energy_flux = @. R_n +
               conditions.lhf * r_ae / (r_soil + r_ae) +
               conditions.shf
            @test computed_energy_flux == expected_energy_flux
        end
    end
end
