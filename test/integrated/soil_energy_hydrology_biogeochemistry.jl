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

        (atmos, radiation) = prescribed_analytic_forcing(FT)
        forcing = (; atmos, radiation)
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
            :thermal_state,
            :SW_d,
            :LW_d,
            :cosθs,
            :frac_diff,
            :soc,
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
            Y.soil.ρe_int .=
                Soil.volumetric_internal_energy.(
                    FT(0.0),
                    ρc_s,
                    T,
                    params.earth_param_set,
                )
        end

        function init_co2!(Y, C_0)
            Y.soilco2.C .= C_0
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
        @test all(parent(p.drivers.soc) .== FT(5.0))
        @test p.soil.θ_l ≈ Soil.Biogeochemistry.soil_moisture(
            model.soilco2.drivers.met,
            p,
            Y,
            t0,
            z,
        )
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
            FT;
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
            FT;
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
        @test met_box.θ_a100 == @. Soil.inverse_matric_potential(
            soil_ps_box.hydrology_cm,
            -FT(1),
        ) * (soil_ps_box.ν - soil_ps_box.θ_r) + soil_ps_box.θ_r
        @test met_col.θ_a100 ==
              Soil.inverse_matric_potential(soil_ps_col.hydrology_cm, -FT(1)) *
              (soil_ps_col.ν - soil_ps_col.θ_r) + soil_ps_col.θ_r
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
