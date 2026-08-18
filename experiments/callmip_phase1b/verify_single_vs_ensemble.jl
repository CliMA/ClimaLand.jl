# Phase 1 correctness harness for the CalLMIP multi-site ColumnEnsemble setup.
#
# Extends experiments/integrated/era5/column_ensemble_comparison.jl (duplicated
# columns, one shared forcing) to the CalLMIP case of genuinely different
# forcing per column:
#   1. DK-Sor on a single `Column` vs DK-Sor as column 1 of a 2-site
#      `ColumnEnsemble` (DK-Sor, US-NR1): no cross-column contamination and
#      geometry effects stay at float level (rtol 1e-9, era5 precedent).
#   2. Row-order flip (US-NR1, DK-Sor): the DK-Sor column must be EXACTLY
#      identical regardless of its position in the ensemble.
#   3. Forcing fidelity: evaluating each driver at raw met-file sample times
#      must reproduce the file values per column (row order + interpolation).
#
# Run:  julia --project=.buildkite experiments/callmip_phase1b/verify_single_vs_ensemble.jl
# GPU:  CLIMACOMMS_DEVICE=CUDA julia --project=.buildkite ...

import ClimaComms
ClimaComms.@import_required_backends

using Test
using Dates
import ClimaCore

include(joinpath(@__DIR__, "..", "integrated", "fluxnet", "callmip_forcing.jl"))
include(joinpath(@__DIR__, "..", "integrated", "era5", "comparison_utils.jl"))

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
toml_dict = LP.create_toml_dict(FT)

# UTC offsets from experiments/callmip_phase1b/DATA_MANIFEST.md (verified
# against SWdown sunrise/sunset midpoints; do not derive from longitude)
dksor = read_callmip_met(
    ClimaLand.Artifacts.callmip_phase1_forcing_path("DK-Sor"),
    1,
)
usnr1 = read_callmip_met(
    ClimaLand.Artifacts.callmip_phase1_forcing_path("US-NR1"),
    -7,
)

start_date = DateTime(2010, 7, 1)
stop_date = start_date + Day(30)
Δt = 450.0

@info "Sim A: DK-Sor on a single Column" start_date stop_date Δt
sim_A = build_callmip_simulation(
    [dksor],
    start_date,
    stop_date,
    Δt,
    toml_dict,
    FT;
    domain_type = :column,
)
@time solve!(sim_A)

@info "Sim B: ColumnEnsemble [DK-Sor, US-NR1]"
sim_B = build_callmip_simulation(
    [dksor, usnr1],
    start_date,
    stop_date,
    Δt,
    toml_dict,
    FT,
)
@time solve!(sim_B)

@info "Sim C: ColumnEnsemble [US-NR1, DK-Sor] (flipped order)"
sim_C = build_callmip_simulation(
    [usnr1, dksor],
    start_date,
    stop_date,
    Δt,
    toml_dict,
    FT,
)
@time solve!(sim_C)

@testset "CalLMIP single vs ensemble (different forcing per column)" begin
    for sim in (sim_A, sim_B, sim_C)
        @test sim._integrator.t == sim._integrator.sol.prob.tspan[2]
    end
    Y_A, p_A = sim_A._integrator.u, sim_A._integrator.p
    Y_B, p_B = sim_B._integrator.u, sim_B._integrator.p
    Y_C, p_C = sim_C._integrator.u, sim_C._integrator.p

    @testset "DK-Sor: single Column vs ensemble column 1" begin
        rtol = 1e-9
        diffs_Y = field_diffs(Y_A, Y_B; col2 = 1, name = "Y")
        diffs_p = field_diffs(p_A, p_B; col2 = 1, name = "p")
        report_diffs(diffs_Y; label = "Y: single vs ensemble col 1 (DK-Sor)")
        report_diffs(diffs_p; label = "p: single vs ensemble col 1 (DK-Sor)")
        for diffs in (diffs_Y, diffs_p)
            @test isempty(filter(((_, d),) -> !(d.err <= rtol), diffs))
        end
    end

    @testset "DK-Sor: ensemble col 1 (B) vs ensemble col 2 (C, flipped)" begin
        diffs_Y = field_diffs(Y_B, Y_C; col1 = 1, col2 = 2, name = "Y")
        diffs_p = field_diffs(p_B, p_C; col1 = 1, col2 = 2, name = "p")
        report_diffs(diffs_Y; label = "Y: DK-Sor col 1 of B vs col 2 of C")
        report_diffs(diffs_p; label = "p: DK-Sor col 1 of B vs col 2 of C")
        for diffs in (diffs_Y, diffs_p)
            @test isempty(filter(((_, d),) -> !(d.err <= 0.0), diffs))
        end
    end

    @testset "US-NR1: ensemble col 2 (B) vs ensemble col 1 (C, flipped)" begin
        diffs_Y = field_diffs(Y_B, Y_C; col1 = 2, col2 = 1, name = "Y")
        diffs_p = field_diffs(p_B, p_C; col1 = 2, col2 = 1, name = "p")
        for diffs in (diffs_Y, diffs_p)
            @test isempty(filter(((_, d),) -> !(d.err <= 0.0), diffs))
        end
    end
end

@testset "Forcing fidelity at raw met sample times" begin
    sites = [dksor, usnr1]
    domain = ColumnEnsemble(;
        zlim = FT.((-2, 0)),
        nelements = 10,
        longlat = [FT.((s.long, s.lat)) for s in sites],
    )
    surface_space = domain.space.surface
    forcing = prescribed_forcing_callmip(
        sites,
        surface_space,
        start_date,
        toml_dict,
        FT,
    )
    utc_dates, aligned = align_sites(sites)
    seconds = Float64.(Dates.value.(utc_dates .- start_date)) ./ 1000
    dest = ClimaCore.Fields.zeros(surface_space)

    # Raw driver series in file units vs the TVIs built from them; evaluating
    # at a sample time must reproduce the sample (linear interp at a node)
    checks = (
        ("Tair", forcing.atmos.T, s -> s.Tair),
        ("Wind", forcing.atmos.u, s -> s.Wind),
        ("Qair", forcing.atmos.q, s -> s.Qair),
        ("Psurf", forcing.atmos.P, s -> s.Psurf),
        ("SWdown", forcing.radiation.SW_d, s -> s.SWdown),
        ("LWdown", forcing.radiation.LW_d, s -> s.LWdown),
        ("CO2air", forcing.atmos.c_co2, s -> s.CO2air .* 1e-6),
    )
    n = length(utc_dates)
    ks = unique(round.(Int, range(1, n, length = 9)))
    for (name, tvi, raw) in checks
        for k in ks
            evaluate!(dest, tvi, seconds[k])
            vals = Array(ClimaCore.Fields.field2array(dest))
            expected = [raw(s)[k] for s in aligned]
            @test all(isapprox.(vals, expected; rtol = 1e-12)) ||
                  error("$name mismatch at k=$k: $vals vs $expected")
        end
    end
end
@info "verify_single_vs_ensemble.jl: all testsets finished"
