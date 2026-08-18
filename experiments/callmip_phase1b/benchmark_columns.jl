# Phase 3: throughput (simulated years per wall day) vs number of ensemble
# columns, on the active ClimaComms device. Real 21-site forcing, cycled to
# reach larger N (duplicate columns are physically independent, so scaling is
# representative). 1-day warmup run first (JIT); each N then times a fresh
# 10-day simulation.
#
# Run:  julia --project=.buildkite experiments/callmip_phase1b/benchmark_columns.jl
# GPU:  CLIMACOMMS_DEVICE=CUDA julia --project=.buildkite ...

import ClimaComms
ClimaComms.@import_required_backends

using Dates

include(joinpath(@__DIR__, "forcing_phase1b.jl"))

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
toml_dict = LP.create_toml_dict(FT)
device = ClimaComms.device()

start_date = DateTime(2010, 6, 1)
warm_stop = start_date + Day(1)
stop_date = start_date + Day(10)
Δt = 450.0
window = (start_date - Day(1), stop_date + Day(1))

base = load_calibration_sites(; window)
Ns = [1, 2, 5, 10, 21, 30, 100, 300]

@info "Column-scaling benchmark" device Δt timed_days = 10
sim_warm = build_callmip_simulation(
    base[1:2],
    start_date,
    warm_stop,
    Δt,
    toml_dict,
    FT,
)
solve!(sim_warm) # JIT warmup
results = NamedTuple[]
for N in Ns
    sites_N = [base[mod1(i, length(base))] for i in 1:N]
    sim = build_callmip_simulation(
        sites_N,
        start_date,
        stop_date,
        Δt,
        toml_dict,
        FT,
    )
    GC.gc() # steady allocation baseline before timing
    stats = @timed solve!(sim)
    sypd = (10 / 365.25) / (stats.time / 86400)
    col_sypd = N * sypd
    push!(results, (; N, wall = stats.time, sypd, col_sypd))
    @info "N=$N" wall_s = round(stats.time; digits = 2) SYPD =
        round(sypd; digits = 1) column_SYPD = round(col_sypd; digits = 1)
end
println("\ndevice, N, wall_s_10days, SYPD, column_SYPD")
for r in results
    println(
        "$(nameof(typeof(device))), $(r.N), $(round(r.wall; digits = 2)), ",
        "$(round(r.sypd; digits = 1)), $(round(r.col_sypd; digits = 1))",
    )
end
