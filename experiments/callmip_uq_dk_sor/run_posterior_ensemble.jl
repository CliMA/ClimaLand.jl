"""
Forward-model ensemble driven by the posterior distribution from
`emulate_sample.jl`.

Draws `N_SAMPLES` parameter vectors from the constrained posterior, writes
ClimaCalibrate-compatible parameter TOMLs, then runs the ClimaLand forward
model for each sample using the same `model_interface.jl` as the calibration.

The resulting daily NEE / Qle / Qh diagnostics can be used to construct
empirical uncertainty bands on land-surface fluxes for CalLMIP Phase 1a.

Workflow
--------
1. Load the posterior from `emulate_sample.jl` output
2. Draw N_SAMPLES parameter vectors from the constrained posterior
3. Stage parameter TOMLs in the ClimaCalibrate directory structure
   (iteration_000/member_NNN/parameters.toml)
4. Run the forward model for each sample (SLURM-parallel via ClimaCalibrate)
5. Save a summary of the posterior parameter samples for analysis

Usage
-----
    julia --project=.buildkite \\
          experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl

Required input (adjust POSTERIOR_DIR to the actual path produced by
emulate_sample.jl, e.g. experiments/callmip_uq_dk_sor/output_posterior_uq/):
    posterior_its1to7.jld2  -- constrained posterior samples

Output directory:
    experiments/callmip_uq_dk_sor/output_posterior_ensemble/
        iteration_000/member_NNN/   -- one directory per posterior sample
            parameters.toml         -- sampled parameter values
            daily_diagnostics.jld2  -- NEE, Qle, Qh timeseries
        posterior_parameter_samples.jld2  -- all sampled parameter vectors
"""

using Distributed
import Random
import JLD2
import ClimaCalibrate
import ClimaLand
import EnsembleKalmanProcesses.ParameterDistributions as PD
using LinearAlgebra, Dates

# ── Configuration ─────────────────────────────────────────────────────────────
const SITE_ID    = "DK-Sor"
const DT         = Float64(450)     # model time step (seconds), matches calibration
const N_SAMPLES  = 50               # number of posterior parameter draws to run

const climaland_dir  = pkgdir(ClimaLand)
const cal_dir        = joinpath(climaland_dir, "experiments", "calibrate_dk_sor")  # calibration artifacts
const exp_dir        = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")  # UQ outputs

# Directory produced by emulate_sample.jl — adjust if you changed `case`:
const POSTERIOR_DIR  = joinpath(exp_dir, "output_posterior_uq")

# Output directory for this ensemble (treated as "iteration 0" by ClimaCalibrate):
const OUTPUT_DIR     = joinpath(exp_dir, "output_posterior_ensemble")
const OBS_FILEPATH   = joinpath(cal_dir, "observations.jld2")

# ── Locate and load the posterior ─────────────────────────────────────────────
function find_posterior_file(dir)
    files = filter(f -> startswith(f, "posterior_its") && endswith(f, ".jld2"),
                   readdir(dir))
    isempty(files) &&
        error("No posterior_its*.jld2 found in $dir — run emulate_sample.jl first")
    path = joinpath(dir, last(sort(files)))   # latest file
    @info "Loading posterior from $path"
    return path
end

posterior_data = JLD2.load(find_posterior_file(POSTERIOR_DIR))
constrained_posterior  = posterior_data["constrained_posterior"]   # (n_params × n_mcmc)
constrained_ekp_optimal = posterior_data["constrained_ekp_optimal"]
param_names            = posterior_data["param_names"]

n_params    = size(constrained_posterior, 1)
n_available = size(constrained_posterior, 2)

@info "Posterior loaded: $n_params parameters, $n_available MCMC samples" *
      "\n  Drawing $N_SAMPLES samples for forward ensemble"

length(param_names) == n_params ||
    error("param_names length ($(length(param_names))) ≠ posterior rows ($n_params)")

# ── Draw parameter samples ────────────────────────────────────────────────────
rng = Random.MersenneTwister(2025_03_23)
idx = rand(rng, 1:n_available, N_SAMPLES)
posterior_samples = constrained_posterior[:, idx]   # (n_params × N_SAMPLES)

println("\nEKI optimum (calibration best-estimate):")
for (name, v) in zip(param_names, constrained_ekp_optimal)
    println("  $(rpad(name, 38)) $(round(v; sigdigits = 4))")
end

# ── Stage parameter TOMLs in ClimaCalibrate directory structure ───────────────
# ClimaCalibrate.parameter_path(output_dir, iteration, member) resolves to:
#   output_dir/iteration_NNN/member_MMM/parameters.toml
# We stage all posterior draws as iteration=0, member=1…N_SAMPLES.

isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

function write_parameter_toml(path, names, values)
    open(path, "w") do io
        for (name, val) in zip(names, values)
            # soilCO2_pre_exponential_factor uses "Land" scope;
            # all other parameters use "getindex" scope
            used_in = name == "soilCO2_pre_exponential_factor" ||
                      name == "michaelis_constant"              ||
                      name == "O2_michaelis_constant" ? "[\"Land\"]" : "[\"getindex\"]"
            println(io, "[\"$name\"]")
            println(io, "value = $val")
            println(io, "type = \"float\"")
            println(io, "used_in = $used_in")
            println(io)
        end
    end
end

for m in 1:N_SAMPLES
    member_dir = ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, 0, m)
    isdir(member_dir) || mkpath(member_dir)
    write_parameter_toml(
        ClimaCalibrate.parameter_path(OUTPUT_DIR, 0, m),
        param_names,
        posterior_samples[:, m],
    )
end

# Also write the EKI optimal as member 0 (useful as a reference run)
ekp_optimal_dir  = joinpath(
    ClimaCalibrate.path_to_iteration(OUTPUT_DIR, 0), "member_000",
)
isdir(ekp_optimal_dir) || mkpath(ekp_optimal_dir)
write_parameter_toml(
    joinpath(ekp_optimal_dir, "parameters.toml"),
    param_names,
    constrained_ekp_optimal,
)

# Save the drawn samples for downstream analysis
JLD2.jldsave(
    joinpath(OUTPUT_DIR, "posterior_parameter_samples.jld2");
    param_names             = param_names,
    posterior_samples       = posterior_samples,
    constrained_ekp_optimal = constrained_ekp_optimal,
    sample_indices          = idx,
    n_samples               = N_SAMPLES,
)
@info "Parameter TOMLs written for $N_SAMPLES members + 1 EKI-optimal reference"

# ── Add SLURM workers ─────────────────────────────────────────────────────────
curr_backend = ClimaCalibrate.get_backend()
if curr_backend == ClimaCalibrate.CaltechHPCBackend ||
   curr_backend == ClimaCalibrate.GCPBackend
    @info "SLURM detected — adding workers (ntasks should be ≥ $N_SAMPLES in SLURM script)"
    addprocs(ClimaCalibrate.SlurmManager())
elseif nworkers() == 1
    @info "No SLURM workers detected — running serially (JuliaBackend fallback)"
end

# ── Broadcast configuration to all workers ────────────────────────────────────
@everywhere using Distributed
@everywhere import ClimaLand
@everywhere const SITE_ID     = $SITE_ID
@everywhere const OUTPUT_DIR  = $OUTPUT_DIR
@everywhere const OBS_FILEPATH = $OBS_FILEPATH
@everywhere const DT          = $DT
@everywhere include(
    joinpath(pkgdir(ClimaLand), "experiments", "calibrate_dk_sor", "model_interface.jl"),
)

# ── Run posterior ensemble ────────────────────────────────────────────────────
@info "Launching $N_SAMPLES forward model runs (posterior ensemble)…"

results = pmap(1:N_SAMPLES) do m
    try
        ClimaCalibrate.forward_model(0, m)
        return (member = m, status = :ok)
    catch e
        @error "Forward model failed for member $m" exception = e
        return (member = m, status = :failed)
    end
end

n_ok     = count(r -> r.status == :ok,     results)
n_failed = count(r -> r.status == :failed, results)
@info "Posterior ensemble complete: $n_ok succeeded, $n_failed failed"
@info "Results in: $OUTPUT_DIR"

# ── Collect summary (optional — produces a convenience file for plotting) ──────
@info "Collecting per-member diagnostics…"

member_results = Dict{Int, Any}()
for m in 1:N_SAMPLES
    diag_path = joinpath(
        ClimaCalibrate.path_to_ensemble_member(OUTPUT_DIR, 0, m),
        "daily_diagnostics.jld2",
    )
    isfile(diag_path) || continue
    d = JLD2.load(diag_path)
    member_results[m] = (dates = d["dates"], nee = d["nee"],
                          qle  = d["qle"],   qh  = d["qh"])
end

JLD2.jldsave(
    joinpath(OUTPUT_DIR, "posterior_ensemble_diagnostics.jld2");
    member_results          = member_results,
    posterior_samples       = posterior_samples,
    constrained_ekp_optimal = constrained_ekp_optimal,
    param_names             = param_names,
    n_samples               = N_SAMPLES,
)
@info "Ensemble diagnostics saved → $(joinpath(OUTPUT_DIR, "posterior_ensemble_diagnostics.jld2"))"
@info "Use `analyze_posterior_ensemble.jl` to compute flux uncertainty bands."
