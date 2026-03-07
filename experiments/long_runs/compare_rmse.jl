# Compare RMSEs between two global simulation runs
#
# Usage:
#   julia --project=. experiments/long_runs/compare_rmse.jl <baseline_dir> <changes_dir>
#
# Example:
#   julia --project=. experiments/long_runs/compare_rmse.jl \
#       main_uncal_gpu/global_diagnostics/output_active \
#       changes_uncal_gpu/global_diagnostics/output_active

using ClimaAnalysis
using Dates
using Printf

# --- Parse arguments ---
if length(ARGS) < 2
    error("""
    Usage: julia --project=. experiments/long_runs/compare_rmse.jl <baseline_dir> <changes_dir>

    Example:
      julia --project=. experiments/long_runs/compare_rmse.jl \\
          main_uncal_gpu/global_diagnostics/output_active \\
          changes_uncal_gpu/global_diagnostics/output_active
    """)
end

baseline_dir = ARGS[1]
changes_dir = ARGS[2]

isdir(baseline_dir) || error("Baseline directory not found: $baseline_dir")
isdir(changes_dir) || error("Changes directory not found: $changes_dir")

@info "Comparing RMSEs" baseline = baseline_dir changes = changes_dir

# --- Load simulation directories ---
baseline_simdir = ClimaAnalysis.SimDir(baseline_dir)
changes_simdir = ClimaAnalysis.SimDir(changes_dir)

baseline_vars = Set(ClimaAnalysis.available_vars(baseline_simdir))
changes_vars = Set(ClimaAnalysis.available_vars(changes_simdir))
common_vars = intersect(baseline_vars, changes_vars)

# Variables of interest for energy flux evaluation
energy_vars = ["shf", "lhf", "lwu", "swu", "et", "gpp"]
# Include ghf if present in both
if "ghf" in common_vars
    push!(energy_vars, "ghf")
end
eval_vars = filter(v -> v in common_vars, energy_vars)

@info "Variables available in both runs" vars = sort(collect(common_vars))
@info "Evaluating" vars = eval_vars

# --- Load ERA5 observations (reuse leaderboard data_sources logic) ---
# We access ERA5 data via the same mechanism as the leaderboard
include(joinpath(@__DIR__, "..", "ilamb", "ilamb_conversion.jl"))

# --- Compute RMSEs ---
"""
    compute_annual_mean_rmse(sim_var, obs_var; spin_up_months=12)

Compute the global RMSE averaged over the last 12 months of the simulation,
after discarding `spin_up_months` months of spin-up.
"""
function compute_annual_mean_rmse(sim_var, obs_var; spin_up_months = 12)
    # Remove spin-up
    spinup_cutoff = spin_up_months * 30 * 86400.0
    if ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff
        sim_var =
            ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff)
    end

    # Align time windows
    obs_times = ClimaAnalysis.times(obs_var)
    sim_times = ClimaAnalysis.times(sim_var)
    min_time = maximum(first.((obs_times, sim_times)))
    max_time = minimum(last.((obs_times, sim_times)))

    sim_var = ClimaAnalysis.window(sim_var, "time", left = min_time, right = max_time)
    obs_var = ClimaAnalysis.window(obs_var, "time", left = min_time, right = max_time)

    obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
    obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

    # Compute monthly RMSEs and average
    times = ClimaAnalysis.times(sim_var)
    rmses = [
        ClimaAnalysis.global_rmse(
            ClimaAnalysis.slice(sim_var, time = t),
            ClimaAnalysis.slice(obs_var, time = t),
        ) for t in times
    ]
    return mean(rmses), rmses
end

using Statistics: mean

# --- Main comparison ---
println("\n" * "="^80)
println("  RMSE COMPARISON: baseline (main_uncal) vs changes (changes_uncal)")
println("="^80)
@printf("%-8s  %12s  %12s  %12s  %8s\n",
    "Var", "Baseline", "Changes", "Δ RMSE", "Change%")
println("-"^60)

results = Dict{String, NamedTuple}()

for short_name in eval_vars
    try
        # Load simulation data
        baseline_var = get(baseline_simdir; short_name)
        changes_var = get(changes_simdir; short_name)

        # We need observation data — use ERA5 leaderboard paths
        # The leaderboard's get_obs_var_dict handles this, but we're outside the ext.
        # For now, compute sim-vs-sim difference and note that obs comparison
        # requires running inside the ext context.

        # Compute sim-only global means for quick comparison
        baseline_times = ClimaAnalysis.times(baseline_var)
        changes_times = ClimaAnalysis.times(changes_var)

        # Use last 12 months after 12 months spin-up
        spinup_cutoff = 12 * 30 * 86400.0

        if baseline_times[end] >= spinup_cutoff
            baseline_var = ClimaAnalysis.window(baseline_var, "time", left = spinup_cutoff)
        end
        if changes_times[end] >= spinup_cutoff
            changes_var = ClimaAnalysis.window(changes_var, "time", left = spinup_cutoff)
        end

        # Align times
        bt = ClimaAnalysis.times(baseline_var)
        ct = ClimaAnalysis.times(changes_var)
        min_t = maximum(first.((bt, ct)))
        max_t = minimum(last.((bt, ct)))
        baseline_var = ClimaAnalysis.window(baseline_var, "time", left = min_t, right = max_t)
        changes_var = ClimaAnalysis.window(changes_var, "time", left = min_t, right = max_t)

        # Compute per-month global RMSE between baseline and changes
        times = ClimaAnalysis.times(baseline_var)
        rmse_diff = [
            ClimaAnalysis.global_rmse(
                ClimaAnalysis.slice(baseline_var, time = t),
                ClimaAnalysis.slice(changes_var, time = t),
            ) for t in times
        ]

        # Compute global means
        baseline_means = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(baseline_var, time = t),
                ),
            ).data[] for t in times
        ]
        changes_means = [
            ClimaAnalysis.weighted_average_lonlat(
                ClimaAnalysis.apply_oceanmask(
                    ClimaAnalysis.slice(changes_var, time = t),
                ),
            ).data[] for t in times
        ]

        avg_baseline = mean(baseline_means)
        avg_changes = mean(changes_means)
        avg_rmse_diff = mean(rmse_diff)

        @printf("%-8s  %12.3f  %12.3f  %12.3f  %+7.1f%%\n",
            short_name, avg_baseline, avg_changes, avg_rmse_diff,
            100.0 * (avg_changes - avg_baseline) / abs(avg_baseline + 1e-10))

        results[short_name] = (
            baseline_mean = avg_baseline,
            changes_mean = avg_changes,
            rmse_between = avg_rmse_diff,
        )
    catch e
        @warn "Could not process $short_name" exception = e
    end
end

println("-"^60)
println("\nNote: The numbers above are global lonlat-weighted means (W/m² or native units)")
println("      and RMSE between the two runs (not vs observations).")
println("      For full obs-based RMSE comparison, see the leaderboard plots in each run's output dir.")
println("\n      Leaderboard plots location:")
println("        Baseline: $(basename(baseline_dir))/../ERA5_global_rmse_and_bias_graphs.png")
println("        Changes:  $(basename(changes_dir))/../ERA5_global_rmse_and_bias_graphs.png")
