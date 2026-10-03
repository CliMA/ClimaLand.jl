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

Known issue: the GP emulator can fail near the prior mean; posterior coverage
may be < 50%. This script reports coverage diagnostics. If coverage < 0.5, a
warning is printed and the EKI final mean is used as the fallback posterior.

Usage
-----
  julia --threads=8 --project=.buildkite \
        experiments/callmip_dksor/emulate_sample.jl
"""

using JLD2
using LinearAlgebra, Statistics, Random
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
# Respect CALLMIP_OUTDIR (matches run_calibration.jl), so CES reads the same
# directory the calibration wrote to (e.g. an E2E/test outdir). Falls back to the
# script-relative default for the standard full run.
const OUTDIR        = get(ENV, "CALLMIP_OUTDIR", joinpath(@__DIR__, "output_calibration"))
mkpath(OUTDIR)

# ── Load priors (must match run_calibration.jl exactly) ──────────────────────
include(joinpath(@__DIR__, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]
@info "Prior: $(length(param_names)) parameters"

# ── Source the training (u, G) pairs + the conditioning observation ───────────
# HYBRID path (preferred): run_calibration.jl ran a minibatched calibration, then a
# separate full-obs "CES training pass" whose (u_unconstrained, G_full) pairs have a
# CONSTANT dimension. We load those directly (get_training_points can't stack the
# variable-dim minibatched G). FALLBACK: full-batch ekp history via get_training_points.
ces_train_path = joinpath(OUTDIR, "ces_training.jld2")
ekp_path       = joinpath(OUTDIR, "ekp_checkpoint.jld2")

if isfile(ces_train_path)
    @info "Hybrid CES: loading constant-dim full-obs pairs from $ces_train_path"
    d = JLD2.load(ces_train_path)
    u_all, G_all = d["u_train"], d["G_train"]
    combined_obs = d["combined_obs"]               # full 396-dim EKP.Observation
    init_sample  = d["u_mean_final"]               # MCMC init (unconstrained post. mean)
    obs_series_mcmc = EKP.ObservationSeries(combined_obs)

    # Drop NaN/Inf columns, then split a held-out test set (last ~⅕, ≥1 col).
    good   = findall(j -> all(isfinite, G_all[:, j]), 1:size(G_all, 2))
    @info "  $(length(good))/$(size(G_all,2)) finite training columns"
    u_all, G_all = u_all[:, good], G_all[:, good]
    n      = size(G_all, 2)
    n_test = clamp(div(n, 5), 1, n - 1)
    test_c = (n - n_test + 1):n
    train_c = 1:(n - n_test)
    train_pairs = PairedDataContainer(u_all[:, train_c], G_all[:, train_c])
    test_pairs  = PairedDataContainer(u_all[:, test_c],  G_all[:, test_c])
    # Output encoder: data-driven decorrelation (decorrelate_sample_cov) with
    # retain_var < 1 to TRUNCATE to the leading-variance dims. This removes the flat
    # GP "stripes" — low-variance output dims a thin training set can't fit and that
    # otherwise leave params under-constrained (bimodal marginals). Input: standardize.
    # The realistic obs noise still drives the MCMC likelihood (obs_series_mcmc), so
    # the coverage fix is preserved.
    RETAIN = parse(Float64, get(ENV, "CALLMIP_CES_RETAIN", "0.99"))
    @info "  Encoder: decorrelate_sample_cov, output retain_var=$RETAIN (truncated)"
    emu_encoder_args = (; encoder_schedule = [
        (Utilities.decorrelate_sample_cov(), "in"),
        (Utilities.decorrelate_sample_cov(; retain_var = RETAIN), "out"),
    ])
    USE_HYBRID = true
else
    isfile(ekp_path) || error("Neither $ces_train_path nor $ekp_path found.\n" *
                              "Run run_calibration.jl first.")
    ekp = JLD2.load(ekp_path, "ekp")
    N_iters = length(EKP.get_u(ekp)) - 1   # completed iterations
    N_ens   = EKP.get_N_ens(ekp)
    @info "EKP loaded: $N_iters completed iterations, ensemble size $N_ens"
    train_iters = 1:max(1, N_iters - 1)
    test_iters  = [N_iters]
    train_pairs = Utilities.get_training_points(ekp, train_iters)
    test_pairs  = Utilities.get_training_points(ekp, test_iters)
    emu_encoder_args = (; encoder_kwargs = Utilities.encoder_kwargs_from(EKP.get_observation_series(ekp)))
    obs_series_mcmc = EKP.get_observation_series(ekp)
    init_sample     = EKP.get_u_mean(ekp, N_iters)
    USE_HYBRID = false
end

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
# emu_encoder_args set above (hybrid: truncating decorrelate_sample_cov schedule;
# fallback: encoder_kwargs_from the ekp observation series).
emulator_path = joinpath(OUTDIR, "emulator.jld2")
emulator = if isfile(emulator_path)
    @info "Loading cached emulator from $emulator_path"
    JLD2.load(emulator_path, "emulator")
else
    @info "Training GP emulator (output dim $(size(get_outputs(train_pairs), 1)))…"
    mlt = GaussianProcess(GPJL(); noise_learn = false)
    em  = Emulator(mlt, train_pairs; emu_encoder_args..., verbose = true)
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

# ── Emulator predictive-variance ('y_var') diagnostic ────────────────────────
# Compare the GP predictive std at the test points against the emulator's
# prediction error. A well-calibrated emulator has error ≈ its own predicted std
# (ratio ~1). Ratio ≫1 ⇒ emulator is over-confident (its variance is too small,
# the usual cause of MCMC under-dispersion / low coverage); ≪1 ⇒ over-cautious.
try
    _, y_var_test = Emulators.predict(emulator, get_inputs(test_pairs);
                                      transform_to_real = true)
    # y_var_test may be an (out_dim × N) matrix or a vector of covariance
    # matrices — collect the per-point variances robustly in either case.
    vars = y_var_test isa AbstractVector ?
        reduce(vcat, [vec(diag(v isa AbstractMatrix ? v : reshape([v], 1, 1))) for v in y_var_test]) :
        vec(y_var_test)
    pred_std = sqrt.(max.(Float64.(vars), 0.0))
    abs_err  = abs.(vec(y_true_test) .- vec(y_pred_test))
    mean_std = mean(filter(isfinite, pred_std))
    mean_err = mean(filter(isfinite, abs_err))
    @info "Emulator predictive-variance (y_var) diagnostic:"
    @info "  mean GP predictive std  = $(round(mean_std; sigdigits=3))"
    @info "  mean |prediction error| = $(round(mean_err; sigdigits=3))"
    @info "  error/std ratio = $(round(mean_err / max(mean_std, eps()); sigdigits=3)) (≈1 well-calibrated; ≫1 over-confident → low coverage)"
catch e
    @warn "y_var diagnostic skipped ($e)"
end

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
# CES is Calibrate→Emulate→Sample: the Sample step is anchored on the Calibrate (EKI)
# result. On the NEE-unidentifiable AR/HR directions the likelihood is flat, so an
# uninformative prior lets the emulator-MCMC drift unphysical (AR > GPP, NEE source). Using
# the EKI posterior (per-parameter mean ± std of the final ensemble) as the MCMC prior
# regularizes those sloppy directions and keeps the CES posterior physical and consistent
# with the calibration. Set CALLMIP_CES_EKI_PRIOR=false to use the uninformative prior.
mcmc_prior = prior
if get(ENV, "CALLMIP_CES_EKI_PRIOR", "true") == "true" && isfile(ekp_path)
    ekp_cal = JLD2.load(ekp_path, "ekp")
    ϕ = Float64.(EKP.get_ϕ_final(prior, ekp_cal))           # constrained final ensemble (11 × N)
    mcmc_prior, _ = build_dk_sor_eki_prior(vec(mean(ϕ; dims = 2)), vec(std(ϕ; dims = 2)))
    @info "CES Sample step anchored on the EKI (Calibrate) posterior — regularizes the AR/HR ridge"
end

# init_sample (unconstrained posterior mean) and obs_series_mcmc (full-obs
# conditioning) were set above for both the hybrid and fallback paths.
mcmc = MCMCWrapper(
    pCNMHSampling(),
    obs_series_mcmc,
    mcmc_prior,
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
    mcmc_prior, posterior_samples_raw)
eki_mean_constrained  = PD.transform_unconstrained_to_constrained(mcmc_prior, init_sample)

# NB: jldsave (not JLD2.save) for the keyword/semicolon form — JLD2.save(path; kw=val)
# silently writes NO file (and does not error), which previously broke the whole
# downstream UQ chain.
jldsave(joinpath(OUTDIR, "posterior_samples.jld2");
    constrained_posterior, param_names)
@info "Posterior samples saved"

# ── Coverage diagnostic ───────────────────────────────────────────────────────
# Check that ~90% of the prior intervals contain the posterior mean.
# If coverage < 50%, the GP emulator likely failed.
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
jldsave(joinpath(OUTDIR, "posterior_mean.jld2");
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
# Each panel is drawn in its OWN x-range with PEAK-NORMALISED step histograms.
# (A prior `density!` KDE produced degenerate autolimits for parameters whose
# constrained scale is far from O(1) — e.g. soilCO2_reference_rate ~1e-6,
# activation_energy ~1e5 — which left every panel but the first blank.)
let
    n  = length(param_names)
    nc = 5
    nr = ceil(Int, n / nc)
    fig = Figure(size = (nc * 320, nr * 230))
    prior_ens = PD.transform_unconstrained_to_constrained(
        prior, EKP.construct_initial_ensemble(Random.MersenneTwister(0), prior, 5000))
    peaknorm(x; nb = 40) = begin
        lo, hi = minimum(x), maximum(x); hi <= lo && (hi = lo + eps(lo) + 1e-30)
        w = (hi - lo) / nb; ctr = [lo + (b - 0.5) * w for b in 1:nb]; cnt = zeros(nb)
        for v in x; cnt[clamp(Int(fld(v - lo, w)) + 1, 1, nb)] += 1; end
        m = maximum(cnt); m == 0 && (m = 1.0); (ctr, cnt ./ m)
    end
    for (i, name) in enumerate(param_names)
        r, c = divrem(i - 1, nc)
        ax   = Axis(fig[r + 1, c + 1]; title = name, titlesize = 10,
                    yticksvisible = false, yticklabelsvisible = false)
        prx, pry = peaknorm(vec(prior_ens[i, :]))
        pox, poy = peaknorm(vec(constrained_posterior[i, :]))
        band!(ax, prx, zeros(length(prx)), pry; color = (:gray, 0.35))
        band!(ax, pox, zeros(length(pox)), poy; color = (:steelblue, 0.55))
        vlines!(ax, [eki_mean_constrained[i]]; color = :red, linestyle = :dash, linewidth = 1.5)
        vlines!(ax, [posterior_mean_vec[i]];   color = :black, linewidth = 1.5)
        xs = vcat(vec(prior_ens[i, :]), vec(constrained_posterior[i, :]))
        xlims!(ax, minimum(xs), maximum(xs)); ylims!(ax, 0, 1.1)
    end
    Legend(fig[nr + 1, 1:nc],
        [PolyElement(color = (:gray, 0.35)), PolyElement(color = (:steelblue, 0.55)),
         LineElement(color = :red, linestyle = :dash), LineElement(color = :black)],
        ["prior", "posterior", "EKI mean", "posterior mean"];
        orientation = :horizontal, framevisible = false)
    save(joinpath(OUTDIR, "posterior_marginals.png"), fig)
    @info "Posterior marginals → $(joinpath(OUTDIR, "posterior_marginals.png"))"
end

@info "=== CES UQ pipeline complete ==="
