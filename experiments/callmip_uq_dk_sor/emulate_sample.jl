"""
CalibrateEmulateSample (CES) v1.0.0 uncertainty quantification for the DK-Sor
single-site calibration.

Reads the EKP object saved by `run_calibration.jl`, trains a Gaussian-Process
(GP) emulator on the UKI (parameter ↔ model-output) pairs, then runs pCN-MH
MCMC to obtain an approximate posterior over the 14 calibrated parameters.

Workflow
--------
1. Load the latest EKP from the calibration output directory
2. Recreate the prior (must mirror `run_calibration.jl` exactly)
3. Load the observation noise covariance from `observations.jld2`
4. Split EKI iterations into training / test sets
5. Train GP emulator with input + output dimension reduction
6. Validate emulator on the held-out test set
7. Run MCMC (pCN-MH) to approximate the posterior
8. Save emulator, chain, and posterior; plot prior vs posterior marginals

For CalLMIP Phase 1a: run this script after `run_calibration.jl` completes,
then use `run_posterior_ensemble.jl` to propagate posterior uncertainty forward
through the land model.

Usage
-----
    julia --project=experiments/callmip_uq_dk_sor \\
          experiments/callmip_uq_dk_sor/emulate_sample.jl

Required input files (relative to `pkgdir(ClimaLand)`):
    experiments/calibrate_dk_sor/output/   -- calibration output directory
        iteration_NNN/eki_file.jld2        -- EKP objects from calibration
    experiments/calibrate_dk_sor/observations.jld2  -- obs + noise covariance

Tunable settings (top of script):
    case           -- label for the output directory
    N_train        -- number of EKI iterations used for GP training
    skip           -- stride over iterations (1 = use all)
    retain_var_in  -- fraction of input variance retained after PCA
    retain_var_out -- fraction of output variance retained after PCA
    chain_length   -- number of MCMC steps
    init_stepsize  -- initial step size for step-size optimisation
"""

using JLD2, TOML
using Plots
using Distributions, LinearAlgebra, Statistics, Random

using CalibrateEmulateSample
using CalibrateEmulateSample.Utilities
using CalibrateEmulateSample.ParameterDistributions
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.MarkovChainMonteCarlo
using CalibrateEmulateSample.DataContainers
import CalibrateEmulateSample.EnsembleKalmanProcesses as EKP
import CalibrateEmulateSample.EnsembleKalmanProcesses.ParameterDistributions as PD

import ClimaLand

rng = Random.MersenneTwister(60732)

# ── Tunable settings ──────────────────────────────────────────────────────────
case           = "posterior_uq"   # output subdirectory label
N_train        = 8                # EKI iterations used to train the emulator (10-iter calibration → max=8)
skip           = 1                # stride over iterations (1 = all)
retain_var_in  = 0.99             # fraction of input variance kept after PCA
retain_var_out = 0.95             # fraction of output variance kept after PCA
chain_length   = 100_000          # MCMC chain length
init_stepsize  = 0.1              # initial step size for optimisation (pCN β)

# ── Paths ─────────────────────────────────────────────────────────────────────
const climaland_dir  = abspath(joinpath(@__DIR__, "..", ".."))
const cal_dir        = joinpath(climaland_dir, "experiments", "calibrate_dk_sor")  # calibration artifacts
const exp_dir        = joinpath(climaland_dir, "experiments", "callmip_uq_dk_sor")  # UQ outputs
const cal_output_dir = joinpath(cal_dir, "output")
const obs_filepath   = joinpath(cal_dir, "observations.jld2")

out_dir = joinpath(exp_dir, "output_$(case)")
isdir(out_dir) || mkdir(out_dir)

# Archive a copy of this script for reproducibility
cp(abspath(@__FILE__), joinpath(out_dir, basename(@__FILE__)); force = true)

@info "=== CalibrateEmulateSample UQ: DK-Sor ===" *
      "\n  case   : $case" *
      "\n  outputs: $out_dir"

# ── Priors ────────────────────────────────────────────────────────────────────
# Loaded from experiments/calibrate_dk_sor/priors.jl — the same file used by
# run_calibration.jl.  This eliminates the risk of the two scripts diverging.
#
include(joinpath(cal_dir, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]

@info "Prior built for $(length(param_names)) parameters:" *
      "\n  $(join(param_names, "\n  "))"

# ── Load calibration EKP ──────────────────────────────────────────────────────
# The calibration saves eki_file.jld2 at each iteration under:
#   cal_output_dir/iteration_NNN/eki_file.jld2
# We search for the highest completed iteration.

function find_latest_ekp_path(cal_out_dir)
    iter = -1
    while isfile(
        joinpath(cal_out_dir, "iteration_$(lpad(iter + 1, 3, '0'))", "eki_file.jld2"),
    )
        iter += 1
    end
    iter == -1 &&
        error("No eki_file.jld2 found in $cal_out_dir. Run run_calibration.jl first.")
    path = joinpath(cal_out_dir, "iteration_$(lpad(iter, 3, '0'))", "eki_file.jld2")
    @info "Loading EKP from iteration $iter: $path"
    return path
end

ekp_path = find_latest_ekp_path(cal_output_dir)
ekp      = JLD2.load_object(ekp_path)

N_iterations_max = length(EKP.get_u(ekp))
@info "EKP loaded: $N_iterations_max stored parameter ensembles " *
      "($(N_iterations_max - 1) completed UKI iterations, " *
      "ensemble size = $(EKP.get_N_ens(ekp)))"

# ── Load observation noise covariance ─────────────────────────────────────────
obs_data  = JLD2.load(obs_filepath)
noise_cov = obs_data["noise_cov"]   # Diagonal{Float64,…} saved by generate_observations.jl

Γ = Matrix(noise_cov)

@info "Observation vector length: $(length(obs_data["y_obs"])) " *
      "(3 × $(length(obs_data["obs_dates"])) valid days)"

# ── Training / test split ─────────────────────────────────────────────────────
# get_training_points(ekp, i) returns the pair (u_{i-1}, G_i).
# We cannot use the final stored u (index N_iterations_max - 1) because no
# corresponding G_final was computed after the last EKI update.
train_iterations = 1:skip:min(skip * N_train, N_iterations_max - 2)
test_iterations  = setdiff(1:(N_iterations_max - 1), train_iterations)

@info "Training iterations : $(collect(train_iterations))"
@info "Test iterations     : $(collect(test_iterations))"

train_pairs = Utilities.get_training_points(ekp, train_iterations)
test_pairs  = Utilities.get_training_points(ekp, test_iterations)

# ── Prune NaN columns ─────────────────────────────────────────────────────────
# Prune sigma-points with NaN outputs or extreme outlier values.
# Extreme sigma-points (forward-model blow-ups) dominate the SVD basis and
# corrupt the decorrelated GP emulator even if they are not NaN.
# threshold=1000 cleanly separates the observed blow-ups (max ≈ 4k–140k in
# iters 1–3) from well-behaved members (max ≤ 175 in iters 4–9).
function prune_bad_columns(pairs::PairedDataContainer; threshold::Real = 1_000)
    out         = get_outputs(pairs)
    col_max     = maximum(abs.(out); dims = 1)[:]
    bad_nan     = findall(sum(isnan.(out); dims = 1)[:] .> 0)
    bad_extreme = findall(col_max .> threshold)
    bad         = union(bad_nan, bad_extreme)
    isempty(bad_nan)     || @warn "  Pruning $(length(bad_nan)) NaN columns"
    isempty(bad_extreme) || @warn "  Pruning $(length(bad_extreme)) extreme columns " *
                                  "(max |G| > $threshold; max found = $(round(maximum(col_max[bad_extreme]), sigdigits=4)))"
    keep = setdiff(1:size(out, 2), bad)
    return PairedDataContainer(get_inputs(pairs)[:, keep], out[:, keep])
end

train_pairs = prune_bad_columns(train_pairs)
test_pairs  = prune_bad_columns(test_pairs)

@info "After pruning: $(size(get_inputs(train_pairs), 2)) train columns, " *
      "$(size(get_inputs(test_pairs), 2)) test columns"

# ── Build and train GP emulator ───────────────────────────────────────────────
gppackage = GPJL()
mlt       = GaussianProcess(gppackage; noise_learn = false)

emulator_path = joinpath(out_dir, "emulator_its$(first(train_iterations))to$(last(train_iterations)).jld2")

if isfile(emulator_path)
    @info "Loading cached emulator from $emulator_path (skipping GP training)"
    emulator = JLD2.load(emulator_path, "emulator")
else
    @info "Building GP emulator (output dim before reduction: " *
          "$(size(get_outputs(train_pairs), 1)))…"

    # CES v1.0.0 encoder_schedule API (replaces deprecated obs_noise_cov /
    # normalize_inputs / decorrelate / retained_svd_frac kwargs from v0.7).
    # decorrelate_sample_cov  : PCA on EKP input sample covariance, retaining
    #                           retain_var_in fraction of variance.
    # decorrelate_structure_mat : SVD on the output structure matrix (Γ),
    #                              retaining retain_var_out fraction.
    encoder_schedule = [
        (decorrelate_sample_cov(retain_var = retain_var_in),    "in"),
        (decorrelate_structure_mat(retain_var = retain_var_out), "out"),
    ]
    emulator = Emulator(
        mlt,
        train_pairs;
        encoder_schedule = deepcopy(encoder_schedule),
        encoder_kwargs   = (; obs_noise_cov = Γ),
        verbose          = true,
    )
    optimize_hyperparameters!(emulator)

    JLD2.save(emulator_path, "emulator", emulator)
    @info "Emulator saved → $emulator_path"
end

# ── Emulator validation ───────────────────────────────────────────────────────
# Per Ollie's review: go beyond MSE — visualize predicted vs actual to catch
# physical/structural issues that aggregate metrics can hide.
y_mean_train, y_var_train = Emulators.predict(emulator, get_inputs(train_pairs); transform_to_real = true)
y_mean_test,  y_var_test  = Emulators.predict(emulator, get_inputs(test_pairs);  transform_to_real = true)

y_true_train = get_outputs(train_pairs)
y_true_test  = get_outputs(test_pairs)

mse_train = mean((y_true_train .- y_mean_train) .^ 2)
mse_test  = mean((y_true_test  .- y_mean_test)  .^ 2)

# Per-output-dimension RMSE (in the decorrelated/reduced space)
rmse_per_dim_train = sqrt.(mean((y_true_train .- y_mean_train) .^ 2; dims = 2))
rmse_per_dim_test  = sqrt.(mean((y_true_test  .- y_mean_test)  .^ 2; dims = 2))

println("\n╔═══════════════════════════════════════╗")
println("║        Emulator validation            ║")
println("╠═══════════════════════════════════════╣")
println("║  GP MSE (train): $(lpad(round(mse_train, sigdigits=4), 15))      ║")
println("║  GP MSE (test) : $(lpad(round(mse_test,  sigdigits=4), 15))      ║")
println("╚═══════════════════════════════════════╝")
println("NOTE: A test MSE much larger than the train MSE indicates overfitting.")

# ── Validation scatter plots ─────────────────────────────────────────────────
# Plot emulator predicted vs actual (in the reduced/decorrelated output space)
# for train and test sets.  Visual inspection can reveal systematic biases,
# heteroskedastic errors, or outlier GP fits.

let n_out_dims = size(y_mean_train, 1)
    # Summarise all output dimensions into one scatter per train/test
    p_scatter = scatter(
        vec(y_true_train), vec(y_mean_train);
        label       = "train ($(size(y_mean_train, 2)) pts × $n_out_dims dims)",
        alpha       = 0.3,
        markersize  = 2,
        markerstrokewidth = 0,
        xlabel      = "True (decorrelated output)",
        ylabel      = "GP predicted (decorrelated output)",
        title       = "Emulator: predicted vs actual\n(latent output space)",
    )
    scatter!(
        vec(y_true_test), vec(y_mean_test);
        label       = "test ($(size(y_mean_test, 2)) pts × $n_out_dims dims)",
        alpha       = 0.5,
        markersize  = 3,
        markerstrokewidth = 0,
    )
    # 1:1 reference line
    lo = min(minimum(y_true_train), minimum(y_true_test))
    hi = max(maximum(y_true_train), maximum(y_true_test))
    plot!(p_scatter, [lo, hi], [lo, hi]; lc = :black, ls = :dash, label = "1:1")

    scatter_path = joinpath(out_dir,
        "emulator_validation_scatter_its$(first(train_iterations))to$(last(train_iterations)).png")
    savefig(p_scatter, scatter_path)
    @info "Emulator validation scatter plot → $scatter_path"

    # Per output-dimension RMSE bar chart (shows which latent dims are poorly fit)
    p_rmse = bar(
        1:n_out_dims,
        vec(rmse_per_dim_train);
        label  = "train",
        alpha  = 0.7,
        xlabel = "Output dimension (decorrelated)",
        ylabel = "RMSE",
        title  = "Per-dimension RMSE\n(latent output space)",
    )
    bar!(
        1:n_out_dims, vec(rmse_per_dim_test);
        label = "test", alpha = 0.5,
    )
    rmse_path = joinpath(out_dir,
        "emulator_validation_rmse_its$(first(train_iterations))to$(last(train_iterations)).png")
    savefig(p_rmse, rmse_path)
    @info "Emulator per-dimension RMSE plot → $rmse_path"
end

# ── MCMC posterior sampling ───────────────────────────────────────────────────
init_sample = EKP.get_u_mean(ekp, maximum(train_iterations))

# CES v1.0.0: MCMCWrapper expects an ObservationSeries (not a raw vector).
observation_series = EKP.get_observation_series(ekp)

mcmc = MCMCWrapper(
    pCNMHSampling(),
    observation_series,
    prior,
    emulator;
    init_params = init_sample,
)

@info "Optimising MCMC step size (N=2000 test steps)…"
new_step = try
    optimize_stepsize(rng, mcmc; init_stepsize, N = 2000, discard_initial = 0)
catch e
    fallback = init_stepsize / 10
    @warn "optimize_stepsize did not converge ($(e)). Falling back to $fallback."
    fallback
end
@info "Step size selected: $new_step  — sampling chain of length $chain_length…"

chain = MarkovChainMonteCarlo.sample(
    rng,
    mcmc,
    chain_length;
    stepsize        = new_step,
    discard_initial = 2_000,
)
display(chain)

posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)

mcmc_path = joinpath(out_dir, "mcmc_and_chain_its$(first(train_iterations))to$(last(train_iterations)).jld2")
JLD2.save(mcmc_path, "mcmc", mcmc, "chain", chain)
@info "MCMC chain saved → $mcmc_path"

# ── Back-transform to constrained (physical) parameter space ──────────────────
posterior_samples = reduce(vcat, [get_distribution(posterior)[name] for name in get_name(posterior)])
constrained_posterior   = PD.transform_unconstrained_to_constrained(
    prior, posterior_samples,
)
constrained_ekp_optimal = PD.transform_unconstrained_to_constrained(prior, init_sample)

posterior_path = joinpath(
    out_dir,
    "posterior_its$(first(train_iterations))to$(last(train_iterations)).jld2",
)
JLD2.save(
    posterior_path,
    "posterior",               posterior,
    "constrained_posterior",   constrained_posterior,
    "constrained_ekp_optimal", constrained_ekp_optimal,
    "param_names",             param_names,
    "train_iterations",        collect(train_iterations),
)
@info "Posterior saved → $posterior_path"

# ── Summary statistics ────────────────────────────────────────────────────────
posterior_mean = vec(mean(constrained_posterior; dims = 2))
posterior_std  = vec(std(constrained_posterior; dims = 2))

println("\n╔══════════════════════════════════════════════════════════════════╗")
println("║         Posterior summary (constrained / physical space)         ║")
println("╠══════════════════════════════╦═══════════════╦══════════════════╣")
println("║ parameter                    ║  EKI optimal  ║  post μ ± σ      ║")
println("╠══════════════════════════════╬═══════════════╬══════════════════╣")
for (i, name) in enumerate(param_names)
    println(
        "║ ",
        rpad(name, 29),
        "║ ",
        lpad(round(constrained_ekp_optimal[i]; sigdigits = 4), 13),
        "  ║ ",
        rpad("$(round(posterior_mean[i]; sigdigits=4)) ± $(round(posterior_std[i]; sigdigits=3))", 16),
        "║",
    )
end
println("╚══════════════════════════════╩═══════════════╩══════════════════╝")

# ── Plot: prior vs posterior marginals ────────────────────────────────────────
p = plot(prior; fill = :lightgray, label = "prior")
plot!(posterior; fill = :darkblue, alpha = 0.5, label = "posterior")
for (i, sp) in enumerate(p.subplots)
    vline!(sp, [constrained_ekp_optimal[i]];
           lc = :red, ls = :dash, lw = 2, label = "EKI optimum")
    title!(sp, param_names[i]; titlefontsize = 7)
end
plot!(p; size = (1800, 1200), margin = 4Plots.mm)

fig_path = joinpath(
    out_dir,
    "posterior_marginals_its$(first(train_iterations))to$(last(train_iterations)).png",
)
savefig(p, fig_path)
@info "Posterior marginals plot → $fig_path"

@info "=== CES UQ pipeline complete ==="
@info "Next step: run `run_posterior_ensemble.jl` to propagate the posterior " *
      "through the land model for CalLMIP Phase 1a."
