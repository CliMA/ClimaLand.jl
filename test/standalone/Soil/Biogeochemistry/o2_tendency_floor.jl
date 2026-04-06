# Regression test for the O2_f tendency bracket. With CFL-violating explicit
# diffusion (hot/dry column, fine near-surface layers, large Δt) and a steep
# initial O2_f profile, the unlimited tendency drives Y_new below 0 *and*
# above the atmospheric O2 fraction (2-cell oscillation). The bracket
# (`SoilCO2Model.Δt > 0`) caps dY so 0 ≤ Y + Δt*dY ≤ O2_f_atm.

using Test
import ClimaComms
ClimaComms.@import_required_backends
using Dates
using ClimaCore
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil.Biogeochemistry
import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    @testset "O2_f tendency bracket prevents [0, O2_f_atm] excursions, FT = $FT" begin
        toml_dict = LP.create_toml_dict(FT)
        parameters = SoilCO2ModelParameters(toml_dict)

        # Stretched column with a very thin top layer to force CFL violation
        # under the chosen Δt (~100 s). dz_top must be small for diffusion
        # to dominate within one step.
        nelems = 20
        domain = Column(;
            zlim = (FT(-2.0), FT(0.0)),
            nelements = nelems,
            dz_tuple = (FT(0.5), FT(0.005)),
        )

        # Hot, dry conditions → large D_o2, small θ_eff_o2 → CFL bound is tiny
        T_soil = (z, t) -> eltype(z)(320)  # 47 °C
        θ_l = (z, t) -> eltype(z)(0.05)   # very dry
        ν = FT(0.5)
        θ_r = FT(0.0)
        hcm = ClimaLand.Soil.vanGenuchten{FT}(; α = FT(0.1), n = FT(2))
        prescribed_met = PrescribedMet{FT}(T_soil, θ_l, ν, θ_r, hcm)

        atmos = ClimaLand.PrescribedAtmosphere(
            TimeVaryingInput((t) -> 1.0),
            TimeVaryingInput((t) -> 1.0),
            TimeVaryingInput((t) -> 1.0),
            TimeVaryingInput((t) -> 1.0),
            TimeVaryingInput((t) -> 1.0),
            TimeVaryingInput((t) -> 100000.0),
            Dates.now(),
            FT(30),
            toml_dict;
            c_co2 = TimeVaryingInput((t) -> 4.2e-4),
        )
        drivers = SoilDrivers(prescribed_met, atmos)
        sources = (MicrobeProduction{FT}(),)

        Δt_floor = FT(100.0)
        model_floor = SoilCO2Model{FT}(
            domain,
            drivers,
            toml_dict,
            Δt_floor;
            parameters,
            sources,
        )
        # Δt = 0 disables the O2_f tendency floor.
        model_nofloor = SoilCO2Model{FT}(
            domain,
            drivers,
            toml_dict,
            FT(0);
            parameters,
            sources,
        )

        # Steep exponential O2_f IC: ~0.21 at z=0, ~0 at z=-1m. This is the
        # initial condition used in multi_column_soilco2.jl that triggers
        # the bug.
        function set_o2_ic!(Y)
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.O2_f)).z
            k = FT(log(21))
            @. Y.soilco2.O2_f = FT(0.21) * exp(k * z)
            Y.soilco2.CO2 .= FT(1e-3)
            Y.soilco2.SOC .= FT(5.0)
        end

        # Step manually with explicit Euler. The CFL-driven 2-cell oscillation
        # appears within a few steps; 20 steps is enough to catch the bug
        # without long runs in CI.
        nsteps = 20
        O2_f_atm = FT(parameters.O2_f_atm)
        function run_steps(model, Δt_step)
            Y, p, _ = initialize(model)
            set_o2_ic!(Y)
            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, FT(0))
            dY = similar(Y)
            exp_tendency! = make_exp_tendency(model)
            min_o2 = FT(Inf)
            max_o2 = -FT(Inf)
            for _ in 1:nsteps
                # update_aux! mutates Y (clamps O2_f); call before each tendency
                set_initial_cache!(p, Y, FT(0))
                exp_tendency!(dY, Y, p, FT(0))
                @. Y.soilco2.O2_f += Δt_step * dY.soilco2.O2_f
                @. Y.soilco2.CO2 += Δt_step * dY.soilco2.CO2
                min_o2 = min(min_o2, minimum(parent(Y.soilco2.O2_f)))
                max_o2 = max(max_o2, maximum(parent(Y.soilco2.O2_f)))
            end
            return min_o2, max_o2
        end

        min_o2_floor, max_o2_floor = run_steps(model_floor, Δt_floor)
        min_o2_nofloor, max_o2_nofloor = run_steps(model_nofloor, Δt_floor)

        # With the bracket on, O2_f stays inside [0, O2_f_atm] up to round-off.
        # Without it, the CFL-driven oscillation pushes O2_f well below zero
        # AND well above O2_f_atm at adjacent cells.
        tol = FT(10) * eps(FT)
        @test min_o2_floor >= -tol
        @test max_o2_floor <= O2_f_atm + tol
        @test min_o2_nofloor < -FT(0.01)
        @test max_o2_nofloor > O2_f_atm + FT(0.01)
    end
end
