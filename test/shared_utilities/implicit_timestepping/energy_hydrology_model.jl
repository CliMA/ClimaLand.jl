using Test
import ClimaComms
ClimaComms.@import_required_backends
using LinearAlgebra
import ClimaCore: MatrixFields
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using Dates

import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "Full Soil Jacobian entries, Flux BC, FT = $FT" begin

        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(1.43)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.124)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)
        params = Soil.EnergyHydrologyParameters(
            FT;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat,
            S_s,
            θ_r,
        )
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        sources = ()

        zero_water_flux_bc = WaterFluxBC((p, t) -> 0.0)
        zero_heat_flux_bc = HeatFluxBC((p, t) -> 0.0)
        boundary_fluxes = (;
            top = WaterHeatBC(;
                water = zero_water_flux_bc,
                heat = zero_heat_flux_bc,
            ),
            bottom = WaterHeatBC(;
                water = zero_water_flux_bc,
                heat = zero_heat_flux_bc,
            ),
        )
        soil = Soil.EnergyHydrology{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_fluxes,
            sources = sources,
        )

        Y, p, coords = initialize(soil)
        Y.soil.ϑ_l .= FT(0.24)
        Y.soil.θ_i .= FT(0)
        T = FT.(280.0)
        ρc_s = @. Soil.volumetric_heat_capacity(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            params.ρc_ds,
            params.earth_param_set,
        )
        @. Y.soil.ρe_int = Soil.volumetric_internal_energy(
            Y.soil.θ_i,
            ρc_s,
            T,
            params.earth_param_set,
        )
        # We do not set the initial aux state here because we want to test that it is updated correctly in the jacobian correctly.
        jacobian = ClimaLand.FieldMatrixWithSolver(Y)
        jac_tendency! = make_jacobian(soil)
        dtγ = FT(1.0)
        jac_tendency!(jacobian, Y, p, dtγ, FT(0.0))


        K_ic =
            impedance_factor(FT(0) / (FT(0.24) + FT(0) - θ_r), params.Ω) *
            viscosity_factor(FT(280), params.γ, params.γT_ref) *
            hydraulic_conductivity(
                hcm,
                K_sat,
                effective_saturation(ν, FT(0.24), θ_r),
            )
        dz = FT(0.01)
        dψdϑ_ic = dψdϑ(hcm, FT(0.24), ν, θ_r, S_s)
        jac_ϑ_l = jacobian.matrix[
            MatrixFields.@name(soil.ϑ_l),
            MatrixFields.@name(soil.ϑ_l)
        ]
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ϑ_l.entries.:2))[1] .≈
              dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ϑ_l.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        @test Array(parent(jac_ϑ_l.entries.:2))[end] .≈
              dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I

        jac_ρe = jacobian.matrix[
            MatrixFields.@name(soil.ρe_int),
            MatrixFields.@name(soil.ρe_int)
        ]
        κ_ic = thermal_conductivity(
            soil.parameters.κ_dry,
            kersten_number(
                FT(0),
                relative_saturation(FT(0.24), FT(0), ν),
                params.α,
                params.β,
                ν_ss_om,
                ν_ss_quartz,
                ν_ss_gravel,
            ),
            κ_sat(FT(0.24), FT(0), params.κ_sat_unfrozen, params.κ_sat_frozen),
        )
        dz = FT(0.01)
        dTdρ_ic =
            1 / volumetric_heat_capacity(
                FT(0.24),
                FT(0),
                params.ρc_ds,
                params.earth_param_set,
            )
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ρe.entries.:2))[1] .≈
              dtγ * (-κ_ic / dz^2 * dTdρ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ρe.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * κ_ic / dz^2 * dTdρ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        @test Array(parent(jac_ρe.entries.:2))[end] .≈
              dtγ * (-κ_ic / dz^2 * dTdρ_ic) - I

        # Off diagonal blocks, ∂T_ρe_int/∂ϑ_l
        jac_ρeϑ = jacobian.matrix[
            MatrixFields.@name(soil.ρe_int),
            MatrixFields.@name(soil.ϑ_l)
        ]

        ρe_liq_ic = ClimaLand.Soil.volumetric_internal_energy_liq(
            T,
            params.earth_param_set,
        )
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ρeϑ.entries.:2))[1] .≈
              dtγ * (-ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ρeϑ.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        @test Array(parent(jac_ρeϑ.entries.:2))[end] .≈
              dtγ * (-ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I

        @test ~(:dfluxBCdY ∈ propertynames(p.soil))
    end

    @testset "Jacobian boundary terms, Atmos Driven, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)
        zmax = FT(0)
        zmin = FT(-1.5)
        nelems = 150
        domain =
            ClimaLand.Domains.Column(; zlim = (zmin, zmax), nelements = nelems)
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
            earth_param_set,
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
        albedo_dry = fill(FT.((0.2, 0.4)), domain.space.surface)
        albedo_wet = fill(FT.((0.1, 0.3)), domain.space.surface)
        params = ClimaLand.Soil.EnergyHydrologyParameters(
            FT;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat,
            S_s,
            θ_r,
            albedo_dry,
            albedo_wet,
            emissivity,
            z_0m,
            z_0b,
        )
        soil = Soil.EnergyHydrology{FT}(;
            parameters = params,
            domain = domain,
            boundary_conditions = boundary_fluxes,
            sources = (),
        )
        Y, p, coords = initialize(soil)
        ϑ_l0 = FT(0.24)
        θ_i0 = FT(0)
        Y.soil.ϑ_l .= ϑ_l0
        Y.soil.θ_i .= θ_i0
        T = FT(280)
        ρc_s = Soil.volumetric_heat_capacity(
            ϑ_l0,
            θ_i0,
            params.ρc_ds,
            params.earth_param_set,
        )
        Y.soil.ρe_int =
            Soil.volumetric_internal_energy.(
                θ_i0,
                ρc_s,
                T,
                params.earth_param_set,
            )
        t = 0.0
        update_drivers! = make_update_drivers((atmos, radiation))
        update_drivers!(p, t)
        # We do not set the initial aux state here because we want to test that it is updated correctly in the jacobian correctly.
        jacobian = ClimaLand.FieldMatrixWithSolver(Y)
        jac_tendency! = make_jacobian(soil)
        dtγ = FT(1.0)
        jac_tendency!(jacobian, Y, p, dtγ, t)

        K_ic =
            impedance_factor(θ_i0 / (ϑ_l0 + θ_i0 - θ_r), params.Ω) *
            viscosity_factor(T, params.γ, params.γT_ref) *
            hydraulic_conductivity(
                hcm,
                K_sat,
                effective_saturation(ν, ϑ_l0, θ_r),
            )
        dz = FT(0.01)
        dψdϑ_ic = dψdϑ(hcm, ϑ_l0, ν, θ_r, S_s)
        jac_ϑ_l = jacobian.matrix[
            MatrixFields.@name(soil.ϑ_l),
            MatrixFields.@name(soil.ϑ_l)
        ]
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ϑ_l.entries.:2))[1] .≈
              dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ϑ_l.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * K_ic / dz^2 * dψdϑ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        @test Array(parent(jac_ϑ_l.entries.:2))[end] .≈
              dtγ * (-K_ic / dz^2 * dψdϑ_ic) - I

        jac_ρe = jacobian.matrix[
            MatrixFields.@name(soil.ρe_int),
            MatrixFields.@name(soil.ρe_int)
        ]
        κ_ic = thermal_conductivity(
            soil.parameters.κ_dry,
            kersten_number(
                θ_i0,
                relative_saturation(ϑ_l0, θ_i0, ν),
                params.α,
                params.β,
                ν_ss_om,
                ν_ss_quartz,
                ν_ss_gravel,
            ),
            κ_sat(ϑ_l0, θ_i0, params.κ_sat_unfrozen, params.κ_sat_frozen),
        )
        dz = FT(0.01)
        dTdρ_ic =
            1 / volumetric_heat_capacity(
                ϑ_l0,
                θ_i0,
                params.ρc_ds,
                params.earth_param_set,
            )
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ρe.entries.:2))[1] .≈
              dtγ * (-κ_ic / dz^2 * dTdρ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ρe.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * κ_ic / dz^2 * dTdρ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        ∂F∂T = 0
        dfluxBCdY_heat = ∂F∂T * dTdρ_ic
        @test Array(parent(jac_ρe.entries.:2))[end] .≈
              dtγ * (-κ_ic / dz^2 * dTdρ_ic - dfluxBCdY_heat / dz) - I

        # Off diagonal blocks, ∂T_ρe_int/∂ϑ_l
        jac_ρeϑ = jacobian.matrix[
            MatrixFields.@name(soil.ρe_int),
            MatrixFields.@name(soil.ϑ_l)
        ]

        ρe_liq_ic = ClimaLand.Soil.volumetric_internal_energy_liq(
            T,
            params.earth_param_set,
        )
        # Check the main diagonal, entry corresponding to bottom of column
        @test Array(parent(jac_ρeϑ.entries.:2))[1] .≈
              dtγ * (-ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I
        # Check the main diagonal, entries corresponding to interior of domain
        @test all(
            Array(parent(jac_ρeϑ.entries.:2))[2:(end - 1)] .≈
            dtγ * (-2 * ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I,
        )
        # Check the main diagonal, entry corresponding to top of column
        @test Array(parent(jac_ρeϑ.entries.:2))[end] .≈
              dtγ * (-ρe_liq_ic * K_ic / dz^2 * dψdϑ_ic) - I

        @test :dfluxBCdY ∈ propertynames(p.soil)
    end

end
