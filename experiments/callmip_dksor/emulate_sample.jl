"""
Calibrate-Emulate-Sample (CES) uncertainty quantification for DK-Sor CalLMIP.

Trains a Gaussian Process (GP) emulator on the EKI (parameter, G) pairs from
run_calibration.jl, then runs MCMC to obtain an approximate posterior. This
is ~10,000× cheaper than running full ensembles.

Workflow
--------
1. Load EKP checkpoint from output_calibration/ekp_checkpoint.jld2
2. Extract (u, G) training pairs from all EKI iterations
3. Prune NaN / extreme-blowup members from training data
4. Train GP emulator with noise/prior covariance from EKP object
5. Validate emulator (train vs test MSE, scatter plot)
6. Run pCN-MH MCMC (100k steps) with optimised step size
7. Save posterior samples and posterior mean parameters

Outputs (in output_calibration/)
---------------------------------
  emulator.jld2              — trained GP emulator
  mcmc_chain.jld2            — MCMC chain object
  posterior_samples.jld2     — constrained posterior samples (N_params × N_samples)
  posterior_mean.jld2        — posterior mean vector (N_params,) + param names
  posterior_marginals.png    — prior vs posterior marginals plot

Known issue (PR #1693): GP emulator can fail near prior mean; posterior coverage
may be < 50%. This script reports coverage diagnostics. If coverage < 0.5, a
warning is printed and the EKI final mean is used as the fallback posterior.

Usage
-----
  julia --threads=8 --project=.buildkite \
        experiments/callmip_dksor/emulate_sample.jl
"""

using JLD2
using Distributions, LinearAlgebra, Statistics, Random
using CairoMakie

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
const CHAIN_LENGTH   = 100_000
const INIT_STEPSIZE  = 0.1

# ── Paths ─────────────────────────────────────────────────────────────────────
const CLIMALAND_DIR = pkgdir(ClimaLand)
const OUTDIR        = joinpath(CLIMALAND_DIR,
    "experiments/callmip_dksor/output_calibration")
mkpath(OUTDIR)

# ── Load priors (must match run_calibration.jl exactly) ──────────────────────
include(joinpath(@__DIR__, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]
@info "Prior: $(length(param_names)) parameters"

# ── Load EKP checkpoint ───────────────────────────────────────────────────────
ekp_path = joinpath(OUTDIR, "ekp_checkpoint.jld2")
isfile(ekp_path) || error("EKP checkpoint not found: $ekp_path\n" *
                           "Run run_calibration.jl first.")
ekp = JLD2.load(ekp_path, "ekp")

N_iters = length(EKP.get_u(ekp)) - 1   # completed iterations
N_ens   = EKP.get_N_ens(ekp)
@info "EKP loaded: $N_iters completed iterations, ensemble size $N_ens"

# ── Build training/test pairs from EKI (u, G) history ─────────────────────────
# get_training_points(ekp, i) returns the pair (u_{i-1}, G_i).
# Hold out the last iteration as a test set.
train_iters = 1:max(1, N_iters - 1)
test_iters  = [N_iters]

train_pairs = Utilities.get_training_points(ekp, train_iters)
test_pairs  = Utilities.get_training_points(ekp, test_iters)

# ── Prune NaN / extreme outlier columns ───────────────────────────────────────
function prune_bad_columns(pairs::PairedDataContainer; threshold::Real = 1_000)
    out     = get_outputs(pairs)
    bad_nan = findall(sum(isnan.(out); dims = 1)[:] .> 0)
    bad_ext = findall(maximum(abs.(out); dims = 1)[:] .> threshold)
    bad     = union(bad_nan, bad_ext)
    isempty(bad_nan) || @warn "  Pruning $(length(bad_nan)) NaN columns"
    isempty(bad_ext) || @warn "  Pruning $(length(bad_ext)) extreme columns (|G| > $threshold)"
    keep = setdiff(1:size(out, 2), bad)
    PairedDataContainer(get_inputs(pairs)[:, keep], out[:, keep])
end

train_pairs = prune_bad_columns(train_pairs)
test_pairs  = prune_bad_columns(test_pairs)

@info "After pruning: $(size(get_inputs(train_pairs), 2)) train / " *
      "$(size(get_inputs(test_pairs), 2)) test columns"

# ── Build and train GP emulator ───────────────────────────────────────────────
# encoder_kwargs_from(ObservationSeries) extracts obs_noise_cov from the first
# observation (with a warning if covariances differ across minibatches).
# This is the correct API for use with ObservationSeries minibatching.
encoder_kwargs = Utilities.encoder_kwargs_from(EKP.get_observation_series(ekp))

emulator_path = joinpath(OUTDIR, "emulator.jld2")
emulator = if isfile(emulator_path)
    @info "Loading cached emulator from $emulator_path"
    JLD2.load(emulator_path, "emulator")
else
    @info "Training GP emulator (output dim $(size(get_outputs(train_pairs), 1)))…"
    mlt = GaussianProcess(GPJL(); noise_learn = false)
    em  = Emulator(mlt, train_pairs; encoder_kwargs, verbose = true)
    optimize_hyperparameters!(em)
    JLD2.save(emulator_path, "emulator", em)
    @info "Emulator saved → $emulator_path"
    em
end

# ── Emulator validation ───────────────────────────────────────────────────────
y_pred_train, _ = Emulators.predict(emulator, get_inputs(train_pairs); transform_to_real = true)
y_pred_test,  _ = Emulators.predict(emulator, get_inputs(test_pairs);  transform_to_real = true)
y_true_train    = get_outputs(train_pairs)
y_true_test     = get_outputs(test_pairs)

mse_train = mean((y_true_train .- y_pred_train) .^ 2)
mse_test  = mean((y_true_test  .- y_pred_test)  .^ 2)
@info "Emulator validation:"
@info "  MSE train = $(round(mse_train; sigdigits=3))"
@info "  MSE test  = $(round(mse_test;  sigdigits=3))"
(mse_test > 10 * mse_train) &&
    @warn "Test MSE >> train MSE — emulator may be overfitting"

# Validation scatter plot
let
    fig = Figure(size = (600, 600))
    ax  = Axis(fig[1, 1];
               xlabel = "True (decorrelated output)",
               ylabel = "GP predicted",
               title  = "Emulator validation: predicted vs actual")
    scatter!(ax, vec(y_true_train), vec(y_pred_train);
             color = (:blue, 0.3), markersize = 3, label = "train")
    scatter!(ax, vec(y_true_test), vec(y_pred_test);
             color = (:red, 0.5), markersize = 4, label = "test")
    lo = min(minimum(y_true_train), minimum(y_true_test))
    hi = max(maximum(y_true_train), maximum(y_true_test))
    lines!(ax, [lo, hi], [lo, hi]; color = :black, linestyle = :dash, label = "1:1")
    axislegend(ax)
    save(joinpath(OUTDIR, "emulator_validation_scatter.png"), fig)
end
@info "Emulator scatter plot → $(joinpath(OUTDIR, "emulator_validation_scatter.png"))"

# ── MCMC posterior sampling ───────────────────────────────────────────────────
# Initialise at the EKI posterior mean (in unconstrained space)
init_sample = EKP.get_u_mean(ekp, N_iters)

observation_series = EKP.get_observation_series(ekp)

mcmc = MCMCWrapper(
    pCNMHSampling(),
    observation_series,
    prior,
    emulator;
    init_params = init_sample,
)

@info "Optimising MCMC step size (2000 test steps)…"
new_step = try
    optimize_stepsize(rng, mcmc; init_stepsize = INIT_STEPSIZE, N = 2000,
                      discard_initial = 0)
catch e
    fallback = INIT_STEPSIZE / 10
    @warn "optimize_stepsize failed ($e). Falling back to $fallback."
    fallback
end
@info "Step size: $new_step. Sampling $CHAIN_LENGTH steps…"

chain = MarkovChainMonteCarlo.sample(
    rng, mcmc, CHAIN_LENGTH;
    stepsize        = new_step,
    discard_initial = 2_000,
)
display(chain)

posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)
JLD2.save(joinpath(OUTDIR, "mcmc_chain.jld2"), "mcmc", mcmc, "chain", chain)
@info "MCMC chain saved"

# ── Back-transform to constrained (physical) parameter space ──────────────────
posterior_samples_raw = reduce(vcat,
    [get_distribution(posterior)[name] for name in get_name(posterior)])
constrained_posterior = PD.transform_unconstrained_to_constrained(
    prior, posterior_samples_raw)
eki_mean_constrained  = PD.transform_unconstrained_to_constrained(prior, init_sample)

JLD2.save(joinpath(OUTDIR, "posterior_samples.jld2");
    constrained_posterior, param_names)
@info "Posterior samples saved"

# ── Coverage diagnostic ───────────────────────────────────────────────────────
# Check that ~90% of the prior intervals contain the posterior mean.
# If coverage < 50%, the GP emulator likely failed (see PR #1693).
coverage_flags = Bool[]
for (i, name) in enumerate(param_names)
    p5, p95 = quantile(constrained_posterior[i, :], [0.05, 0.95])
    push!(coverage_flags, p5 <= eki_mean_constrained[i] <= p95)
end
coverage = mean(coverage_flags)
@info "Posterior coverage (fraction within prior 90% CI): $(round(coverage; digits=2))"

USE_EKI_FALLBACK = coverage < 0.5
if USE_EKI_FALLBACK
    @warn "Coverage $coverage < 0.5 — GP emulator likely failed. " *
          "Using EKI final mean as posterior fallback."
    posterior_mean_vec = eki_mean_constrained
else
    posterior_mean_vec = vec(mean(constrained_posterior; dims = 2))
end

# ── Save posterior mean ───────────────────────────────────────────────────────
JLD2.save(joinpath(OUTDIR, "posterior_mean.jld2");
    posterior_mean    = posterior_mean_vec,
    param_names,
    eki_mean          = eki_mean_constrained,
    coverage,
    used_eki_fallback = USE_EKI_FALLBACK,
)
@info "Posterior mean saved → $(joinpath(OUTDIR, "posterior_mean.jld2"))"

@info "Posterior mean:"
for (i, name) in enumerate(param_names)
    @info "  $(rpad(name, 35)) = $(round(posterior_mean_vec[i]; sigdigits=4))"
end

# ── Plot prior vs posterior marginals ─────────────────────────────────────────
let
    n  = length(param_names)
    nc = 5
    nr = ceil(Int, n / nc)
    fig = Figure(size = (nc * 300, nr * 200))
    for (i, name) in enumerate(param_names)
        r, c = divrem(i - 1, nc)
        ax   = Axis(fig[r + 1, c + 1]; title = name, titlesize = 9)

        post_samp = vec(constrained_posterior[i, :])
        hist!(ax, post_samp; bins = 40, color = (:steelblue, 0.7),
              normalization = :pdf, label = "posterior")

        # Prior density (approximate via 5000 prior samples)
        prior_samp = PD.transform_unconstrained_to_constrained(
            prior, EKP.construct_initial_ensemble(
                Random.MersenneTwister(0), prior, 5000)[i:i, :])
        density!(ax, vec(prior_samp); color = (:gray, 0.3), linestyle = :dash,
                 label = "prior")
        vlines!(ax, [eki_mean_constrained[i]]; color = :red, linestyle = :dash,
                linewidth = 2, label = "EKI mean")
    end
    save(joinpath(OUTDIR, "posterior_marginals.png"), fig)
    @info "Posterior marginals → $(joinpath(OUTDIR, "posterior_marginals.png"))"
end

@info "=== CES UQ pipeline complete ==="
