using Test
import ClimaComms
ClimaComms.@import_required_backends
using LinearAlgebra
import ClimaCore: MatrixFields
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil

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
        jacobian = ImplicitEquationJacobian(Y)
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
    end
end
