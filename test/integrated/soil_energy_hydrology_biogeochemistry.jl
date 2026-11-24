using Test
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams
using ClimaLand
using ClimaLand.Domains: Column, HybridBox
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using Dates

import ClimaParams
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Soil respiration test set, FT = $FT" begin
        toml_dict = ClimaLand.Parameters.create_toml_dict(FT)

        zmax = FT(0)
        zmin = FT(-1)
        nelems = 20
        domain = Column(;
            zlim = (zmin, zmax),
            nelements = nelems,
            longlat = FT.((-118.0, 45.0)),
        )

        (atmos, radiation) = prescribed_analytic_forcing(FT; toml_dict)
        forcing = (; atmos, radiation)
        Δt = FT(450)
        model = LandSoilBiogeochemistry{FT}(forcing, toml_dict, domain)
        @test ClimaComms.context(model) == ClimaComms.context()
        @test ClimaComms.device(model) == ClimaComms.device()

        @test model.soilco2.drivers.met.ν == model.soil.parameters.ν
        @test model.soilco2.drivers.met isa ClimaLand.PrognosticMet
        drivers = ClimaLand.get_drivers(model)
        @test drivers[1] == atmos
        @test drivers[2] == radiation
        Y, p, coords = initialize(model)
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
        function init_soil!(Y, z, params)
            ν = params.ν
            FT = eltype(Y.soil.ϑ_l)
            Y.soil.ϑ_l .= FT(0.33)
            Y.soil.θ_i .= FT(0.0)
            T = FT(279.85)
            ρc_s = @. FT(
                Soil.volumetric_heat_capacity(
                    FT(0.33),
                    FT(0.0),
                    params.ρc_ds,
                    params.earth_param_set,
                ),
            )
            Y.soil.ρe_int .= Soil.volumetric_internal_energy.(
                FT(0.0),
                ρc_s,
                T,
                params.earth_param_set,
            )
        end

        function init_co2!(Y, C_0)
            Y.soilco2.CO2 .= C_0
            Y.soilco2.O2 .= FT(0.23)
            Y.soilco2.SOC .= FT(5.0)    # Initialize SOC
        end

        z = coords.subsurface.z
        init_soil!(Y, z, model.soil.parameters)
        C = FT(4)
        init_co2!(Y, C)
        t0 = FT(0.0)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t0)

        @test p.soil.T ≈ Soil.Biogeochemistry.soil_temperature(
            model.soilco2.drivers.met,
            p,
            Y,
            t0,
            z,
        )
        @test all(parent(Y.soilco2.SOC) .== FT(5.0))
        @test p.soil.θ_l ≈ Soil.Biogeochemistry.soil_moisture(
            model.soilco2.drivers.met,
            p,
            Y,
            t0,
            z,
        )
    end
    @testset "Soil CO2 and O2 initial conditions, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        domain = Column(;
            zlim = (FT(-1), FT(0)),
            nelements = 10,
            longlat = FT.((-118.0, 45.0)),
        )
        (atmos, radiation) = prescribed_analytic_forcing(FT; toml_dict)
        model =
            LandSoilBiogeochemistry{FT}((; atmos, radiation), toml_dict, domain)
        Y, p, coords = initialize(model)
        z = coords.subsurface.z
        (; ν, θ_r, ρc_ds, earth_param_set) = model.soil.parameters
        # Pores full of ice above 0.2 m, half full of water down to 0.5 m, and full
        # of water below, so dissolved storage matters
        frozen = @. z > FT(-0.2)
        @. Y.soil.ϑ_l = ifelse(
            frozen,
            θ_r + FT(0.02),
            ifelse(z < FT(-0.5), ν, θ_r + (ν - θ_r) / 2),
        )
        @. Y.soil.θ_i = ifelse(frozen, ν - Y.soil.ϑ_l, FT(0))
        ρc_s = @. Soil.volumetric_heat_capacity(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            ρc_ds,
            earth_param_set,
        )
        @. Y.soil.ρe_int = Soil.volumetric_internal_energy(
            Y.soil.θ_i,
            ρc_s,
            ifelse(frozen, FT(268), FT(288)),
            earth_param_set,
        )
        Y.soilco2 .= FT(0)
        t0 = FT(0)
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, t0)
        ClimaLand.Simulations.set_soilco2_initial_conditions!(Y, p, model)
        set_initial_cache!(p, Y, t0)

        # The soil air holds the prescribed O2 profile, whatever the water content
        params = model.soilco2.parameters
        O2_f = @. o2_fraction_from_concentration(
            Y.soilco2.O2 / p.soilco2.θ_eff_o2,
            p.soilco2.T,
            p.drivers.P,
            params,
        )
        O2_f_expected = @. FT(0.21) * exp(log(FT(0.21 / 0.15)) * z)
        @test all(
            isapprox.(
                Array(parent(O2_f)),
                Array(parent(O2_f_expected));
                rtol = 1e-3,
            ),
        )
        # and the atmospheric CO2 mole fraction
        R = LP.gas_constant(earth_param_set)
        CO2_f = @. Y.soilco2.CO2 / p.soilco2.θ_eff * R * p.soilco2.T /
           (p.drivers.P * params.M_C)
        c_co2 = Array(parent(p.drivers.c_co2))[1]
        @test all(x -> isapprox(x, c_co2; rtol = 1e-3), Array(parent(CO2_f)))
    end
    @testset "PrognosticMet, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        zmax = FT(0)
        zmin = FT(-1)
        nelems = 10
        xmin = ymin = zmin
        xmax = ymax = zmax
        col = Column(; zlim = (zmin, zmax), nelements = nelems)
        box = HybridBox(;
            xlim = (xmin, xmax),
            ylim = (ymin, ymax),
            zlim = (zmin, zmax),
            nelements = (nelems, nelems, nelems),
        )

        ν = FT(0.556)
        K_sat = FT(0.0443 / 3600 / 100)
        S_s = FT(1e-3)
        vg_n = FT(2.0)
        vg_α = FT(2.6)
        hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)


        soil_ps_col = Soil.EnergyHydrologyParameters(
            toml_dict;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            K_sat,
            hydrology_cm,
            S_s,
            θ_r,
        )

        zero_field = ClimaCore.Fields.zeros(box.space.subsurface)
        vg_fields_to_hcm_field(α::FT, n::FT) where {FT} =
            ClimaLand.Soil.vanGenuchten{FT}(;
                @NamedTuple{α::FT, n::FT}((α, n))...,
            )
        hydrology_cm_field =
            vg_fields_to_hcm_field.(vg_α .+ zero_field, vg_n .+ zero_field)
        ν_field = ν .+ zero_field
        ν_ss_om_field = ν_ss_om .+ zero_field
        ν_ss_quartz_field = ν_ss_quartz .+ zero_field
        ν_ss_gravel_field = ν_ss_gravel .+ zero_field
        K_sat_field = K_sat .+ zero_field
        S_s_field = S_s .+ zero_field
        θ_r_field = θ_r .+ zero_field

        soil_ps_box = Soil.EnergyHydrologyParameters(
            toml_dict;
            ν = ν_field,
            ν_ss_om = ν_ss_om_field,
            ν_ss_quartz = ν_ss_quartz_field,
            ν_ss_gravel = ν_ss_gravel_field,
            K_sat = K_sat_field,
            hydrology_cm = hydrology_cm_field,
            S_s = S_s_field,
            θ_r = θ_r_field,
        )

        met_col = ClimaLand.PrognosticMet(soil_ps_col)
        met_box = ClimaLand.PrognosticMet(soil_ps_box)
        @test met_box.ν == soil_ps_box.ν
        @test met_col.ν == soil_ps_col.ν
        # θ_a100 is the air content left when the matric potential is -1 m
        S_100_col = Soil.effective_saturation(
            soil_ps_col.ν,
            soil_ps_col.ν - met_col.θ_a100,
            soil_ps_col.θ_r,
        )
        @test Soil.matric_potential(soil_ps_col.hydrology_cm, S_100_col) ≈
              FT(-1)
        ψ_100_box = @. Soil.matric_potential(
            soil_ps_box.hydrology_cm,
            Soil.effective_saturation(
                soil_ps_box.ν,
                soil_ps_box.ν - met_box.θ_a100,
                soil_ps_box.θ_r,
            ),
        )
        @test all(parent(ψ_100_box) .≈ FT(-1))
        @test met_box.b ==
              @. Soil.approximate_ψ_S_slope(soil_ps_box.hydrology_cm)
        @test met_col.b == Soil.approximate_ψ_S_slope(soil_ps_col.hydrology_cm)
        vg_m = 1 - 1 / vg_n
        @test Soil.approximate_ψ_S_slope(soil_ps_col.hydrology_cm) ==
              (1 + vg_m) / (vg_n * vg_m^2)
        @test Soil.approximate_ψ_S_slope(
            BrooksCorey{FT}(; c = FT(1.0), ψb = FT(1.0)),
        ) == FT(1)
    end
end
