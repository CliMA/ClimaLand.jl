import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP

@testset "Test specific ClimaAtmos outcomes" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
    thermo_params = param_set.thermo_params
    u_in = (FT(-19.07545), FT(16.88031))
    u_sfc = (FT(0.0), FT(0.0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    Δz = FT(15.000001)
    ΔDSEᵥ = FT(598.96875)
    ts_in =
        TD.PhaseEquil{FT}(1.2595116, 99902.82, 12337.749, 0.0044478197, 275.624)
    ts_sfc = TD.PhaseEquil{FT}(
        1.2544012,
        99335.55,
        11996.086,
        0.0044396375,
        275.1768,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)

    u_in = (FT(-0.168524f0), FT(-0.000566946f0))
    u_sfc = (FT(0.0), FT(0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    Δz = FT(15.000001)
    ΔDSEᵥ = FT(662.09375)
    ts_in = TD.PhaseEquil{FT}(
        1.2605726f0,
        100331.47f0,
        7956.4053f0,
        0.002202735f0,
        276.95068f0,
    )
    ts_sfc = TD.PhaseEquil{FT}(
        1.2499729f0,
        99303.92f0,
        13258.002f0,
        0.0047157165f0,
        276.01752f0,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)

    u_in = (FT(-14.154735), FT(-5.1905923))
    u_sfc = (FT(0.0), FT(0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    Δz = FT(15.000001)
    ΔDSEᵥ = FT(21.34375)
    ts_in = TD.PhaseEquil{FT}(
        1.1730341f0,
        98689.72f0,
        43302.703f0,
        0.012817842f0,
        290.8733f0,
    )
    ts_sfc = TD.PhaseEquil{FT}(
        1.1740736f0,
        98819.375f0,
        43671.336f0,
        0.012941063f0,
        290.97592f0,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)

    u_in = (FT(-13.526638), FT(-8.794365))
    u_sfc = (FT(0.0), FT(0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    ΔDSEᵥ = FT(24.46875)
    ts_in = TD.PhaseEquil{FT}(
        1.1698402f0,
        98647.89f0,
        44855.285f0,
        0.013289474f0,
        291.46088f0,
    )
    ts_sfc = TD.PhaseEquil{FT}(
        1.1708081f0,
        98770.55f0,
        45266.523f0,
        0.013432522f0,
        291.5569f0,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)

    u_in = (FT(-41.34482), FT(-23.609104))
    u_sfc = (FT(0.0), FT(0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    Δz = FT(15.000001)
    ΔDSEᵥ = FT(0.25)
    ts_in = TD.PhaseEquil{FT}(
        1.2182463f0,
        96874.9f0,
        13805.914f0,
        0.0048752176f0,
        276.25174f0,
    )
    ts_sfc = TD.PhaseEquil{FT}(
        1.2197124f0,
        97042.68f0,
        14087.365f0,
        0.004953378f0,
        276.38446f0,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)

    u_in = (FT(-0.75088084), FT(-0.09317328))
    u_sfc = (FT(0.0), FT(0))
    z0m = FT(1.0e-5)
    z0b = FT(1.0e-5)
    Δz = FT(15.000001)
    ΔDSEᵥ = FT(1024.0)
    ts_in = TD.PhaseEquil{FT}(
        1.2317619f0,
        99965.086f0,
        12921.355f0,
        0.0026684932f0,
        282.31366f0,
    )
    ts_sfc = TD.PhaseEquil{FT}(
        1.214932f0,
        98294.87f0,
        21252.703f0,
        0.006637053f0,
        280.76575f0,
    )
    state_sfc = SF.StateValues(FT(0), u_sfc, ts_sfc)
    state_in = SF.StateValues(Δz, u_in, ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    result = SF.surface_conditions(param_set, sc)
    @test result.L_MO > FT(0)
end
