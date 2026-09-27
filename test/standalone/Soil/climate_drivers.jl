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
                  (:lhf, :shf, :vapor_flux_liq, :vapor_flux_ice, :T_sfc)
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
                :q_sfc,
                :PAR_albedo,
                :NIR_albedo,
                :W_gap,
                :r_undercanopy,
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
                Y.soil.ρe_int = Soil.volumetric_internal_energy.(
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
            @test ClimaLand.surface_emissivity(model, Y, p) == emissivity
            PAR_albedo = p.soil.PAR_albedo
            NIR_albedo = p.soil.NIR_albedo
            @test Base.Broadcast.materialize(
                ClimaLand.surface_albedo(model, Y, p),
            ) == PAR_albedo ./ 2 .+ NIR_albedo ./ 2
            T_sfc = p.soil.turbulent_fluxes.T_sfc
            @test parent(ClimaLand.component_temperature(model, Y, p)) ==
                  parent(T_sfc)
            @test all(parent(T_sfc) .> FT(280.0))

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
            # The fluxes stored by the skin solve match a separate solve at the
            # skin temperature to within the tolerance of the Monin-Obukhov
            # solve
            stored = p.soil.turbulent_fluxes
            @test parent(conditions.T_sfc) == parent(stored.T_sfc)
            for name in (:lhf, :shf)
                F_fresh = parent(getproperty(conditions, name))
                F_stored = parent(getproperty(stored, name))
                @test all(
                    abs.(F_fresh .- F_stored) .<
                    FT(0.01) .* (abs.(F_fresh) .+ FT(1)),
                )
            end
            @test all(
                abs.(
                    parent(conditions.vapor_flux_liq) .-
                    parent(stored.vapor_flux_liq),
                ) .<
                FT(0.01) .*
                (abs.(parent(conditions.vapor_flux_liq)) .+ eps(FT)),
            )

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

            expected_water_flux = @. FT(precip(t)) .+ stored.vapor_flux_liq
            @test computed_water_flux == expected_water_flux
            expected_energy_flux = @. R_n_copy +
               stored.lhf +
               stored.shf +
               FT(precip(t)) * Soil.volumetric_internal_energy_liq(
                   FT(T_atmos(t)),
                   earth_param_set,
               )
            @test computed_energy_flux == expected_energy_flux

            # The skin temperature closes the surface energy balance: the
            # atmospheric fluxes match conduction between skin and top layer
            # to within the tolerance of the Monin-Obukhov solve.
            T_top = ClimaLand.Domains.top_center_to_surface(p.soil.T)
            κ_top = ClimaLand.Domains.top_center_to_surface(p.soil.κ)
            r_sfc = @. Δz_top / κ_top
            F_atmos = @. R_n_copy + stored.lhf + stored.shf
            F_scale = @. abs(R_n_copy) + abs(stored.lhf) + abs(stored.shf)
            G_skin = @. (T_top - T_sfc) / r_sfc
            @test all(parent(T_sfc) .!= parent(T_top))
            @test all(
                abs.(parent(F_atmos) .- parent(G_skin)) .<
                FT(0.01) .* parent(F_scale),
            )

            # Over a frozen top cell under strong sunshine and warm air, the
            # skin would rise above the melting point; it is capped at the
            # depressed freezing temperature, and the atmospheric fluxes at
            # that temperature exceed the skin-top conduction (the excess
            # melts ice in the top cell)
            Y_frozen = copy(Y)
            p_frozen = deepcopy(p)
            θ_i_frozen = FT(0.2)
            θ_l_frozen = θ_r + FT(0.05)
            T_frozen = FT(265)
            ρc_frozen = Soil.volumetric_heat_capacity(
                θ_l_frozen,
                θ_i_frozen,
                model.parameters.ρc_ds,
                earth_param_set,
            )
            Y_frozen.soil.ϑ_l .= θ_l_frozen
            Y_frozen.soil.θ_i .= θ_i_frozen
            Y_frozen.soil.ρe_int .= Soil.volumetric_internal_energy.(
                θ_i_frozen,
                ρc_frozen,
                T_frozen,
                earth_param_set,
            )
            set_initial_cache!(p_frozen, Y_frozen, t)
            Tf_sfc = ClimaLand.Domains.top_center_to_surface(
                p_frozen.soil.Tf_depressed,
            )
            T_top_frozen =
                ClimaLand.Domains.top_center_to_surface(p_frozen.soil.T)
            @test all(parent(T_top_frozen) .< parent(Tf_sfc))
            @test all(
                isapprox.(
                    parent(p_frozen.soil.turbulent_fluxes.T_sfc),
                    parent(Tf_sfc);
                    rtol = 10eps(FT),
                ),
            )
            κ_top_frozen =
                ClimaLand.Domains.top_center_to_surface(p_frozen.soil.κ)
            F_frozen = @. p_frozen.soil.R_n +
               p_frozen.soil.turbulent_fluxes.lhf +
               p_frozen.soil.turbulent_fluxes.shf
            G_frozen =
                @. (T_top_frozen - p_frozen.soil.turbulent_fluxes.T_sfc) *
                   κ_top_frozen / Δz_top
            # Net downward atmospheric flux exceeds the conduction from the
            # skin into the soil
            @test all(parent(F_frozen) .< parent(G_frozen))

            # Trace ice in a top cell above the melting point does not cap the
            # skin (the cap would remove the feedback of the surface fluxes on
            # the soil temperature)
            Y_trace = copy(Y)
            p_trace = deepcopy(p)
            Y_trace.soil.θ_i .= FT(1e-6)
            set_initial_cache!(p_trace, Y_trace, t)
            Tf_trace = ClimaLand.Domains.top_center_to_surface(
                p_trace.soil.Tf_depressed,
            )
            @test all(
                parent(
                    ClimaLand.Domains.top_center_to_surface(p_trace.soil.T),
                ) .> parent(Tf_trace),
            )
            @test all(
                parent(p_trace.soil.turbulent_fluxes.T_sfc) .>
                parent(Tf_trace) .+ 1,
            )

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
            @test extrema(dsl ./ d_ds)[2] ≈ ((S_c - S_sfc[1]) / S_c)^evap_p
            # The dry surface layer is bounded by its maximum thickness
            @test all(dsl .<= d_ds)
            _D_vapor = FT(LP.D_vapor(earth_param_set))
            S_sfc = @. ClimaLand.Soil.effective_saturation(ν, θ_sfc, θ_r)
            gsoil = @. ClimaLand.Soil.soil_conductance(
                S_sfc,
                hydrology_cm.S_c,
                d_ds,
                evap_p,
                evap_α,
                _D_vapor,
                ν,
                θ_r,
                FT(0),
            )
            @test gsoil[1] < gsoil[2]
            @test ClimaLand.Soil.soil_conductance(
                FT(0),
                hydrology_cm.S_c,
                d_ds,
                evap_p,
                evap_α,
                _D_vapor,
                ν,
                θ_r,
                FT(0),
            ) == 0
            # The dry layer is air-dry: its tortuosity does not depend on the
            # moisture of the soil below it
            τ_a = ClimaLand.Soil.soil_tortuosity(ν, θ_r)
            @test τ_a ≈ (ν - θ_r)^FT(2.5) / ν
            # Ice fills pore space that vapor would otherwise diffuse through
            @test ClimaLand.Soil.soil_tortuosity(ν, θ_r, FT(0.1)) < τ_a
            # The conductance is the inverse dry-layer resistance, times the
            # availability of mobile liquid water
            dsl1 = ClimaLand.Soil.dry_soil_layer_thickness(
                S_sfc[1],
                S_c,
                d_ds,
                evap_p,
            )
            f_avail1 = S_sfc[1]^2 / (S_sfc[1]^2 + FT(0.01)^2)
            @test gsoil[1] ≈ f_avail1 * _D_vapor * τ_a / dsl1
            # At saturation there is no dry layer: the conductance is unbounded
            @test gsoil[end] > FT(1e3)
            @test issorted(gsoil)
        end
    end
end
