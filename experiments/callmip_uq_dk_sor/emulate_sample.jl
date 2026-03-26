"""
CalibrateEmulateSample (CES) uncertainty quantification for the DK-Sor
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
N_train        = 7                # EKI iterations used to train the emulator
skip           = 1                # stride over iterations (1 = all)
retain_var_in  = 0.99             # fraction of input variance kept after PCA
retain_var_out = 0.95             # fraction of output variance kept after PCA
chain_length   = 100_000          # MCMC chain length
init_stepsize  = 0.1              # initial step size for optimisation

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
# These MUST mirror the `priors` vector in `run_calibration.jl` EXACTLY as it
# was when the EKP was created — 12 parameters (9 canopy + 3 DAMM soilCO2).
# NOTE: root_leaf_nitrogen_ratio and stem_leaf_nitrogen_ratio were added to
# run_calibration.jl after the stored EKP was generated and are NOT calibrated.
# They are handled via hardcoded defaults in run_callmip_simulations.jl.
priors_vec = [
    # Canopy / conductance parameters
    PD.constrained_gaussian("moisture_stress_c",              0.27,     0.15,     0.01,    5.0),
    PD.constrained_gaussian("pmodel_cstar",                   0.43,     0.15,     0.05,    2.0),
    PD.constrained_gaussian("pmodel_β",                      51.0,     20.0,      5.0,  500.0),
    PD.constrained_gaussian("leaf_Cd",                        0.07,     0.04,     0.005,   1.0),
    PD.constrained_gaussian("canopy_z_0m_coeff",              0.02,     0.01,     0.001,   0.3),
    PD.constrained_gaussian("canopy_z_0b_coeff",              0.0007,   0.0003,   1e-5,  0.005),
    PD.constrained_gaussian("canopy_d_coeff",                 0.007,    0.004,    0.001,   0.1),
    PD.constrained_gaussian("canopy_K_lw",                    0.85,     0.25,     0.1,     2.0),
    PD.constrained_gaussian("canopy_emissivity",              0.98,     0.01,     0.9,     1.0),
    # DAMM soil-CO₂
    PD.constrained_gaussian("soilCO2_pre_exponential_factor", 23835.0,  10000.0,  1000.0, 200000.0),
    PD.constrained_gaussian("michaelis_constant",             0.005,    0.003,    1e-4,    0.1),
    PD.constrained_gaussian("O2_michaelis_constant",          0.004,    0.002,    1e-4,    0.1),
]

prior       = PD.combine_distributions(priors_vec)
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
Γ         = Matrix(noise_cov)       # dense form required by the CES encoder

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
# Occasional forward-model failures produce NaN outputs; drop those columns.
function prune_nans(pairs::PairedDataContainer)
    out  = get_outputs(pairs)
    bad  = findall(sum(isnan.(out); dims = 1)[:] .> 0)
    isempty(bad) || @warn "  Pruning $(length(bad)) NaN columns"
    keep = setdiff(1:size(out, 2), bad)
    return PairedDataContainer(get_inputs(pairs)[:, keep], out[:, keep])
end

train_pairs = prune_nans(train_pairs)
test_pairs  = prune_nans(test_pairs)

@info "After NaN pruning: $(size(get_inputs(train_pairs), 2)) train columns, " *
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

    emulator = Emulator(
        mlt,
        train_pairs;
        obs_noise_cov    = Γ,
        normalize_inputs = true,
        decorrelate      = true,
        retained_svd_frac = retain_var_out,
    )
    # Bound the log-length-scales to [-10, 10] and log-variance to [-10, 10] to
    # prevent LBFGS divergence to Inf for nearly-flat output dimensions.
    n_gp_inputs = size(get_inputs(train_pairs), 1)
    kb = [fill(-10.0, n_gp_inputs + 1), fill(10.0, n_gp_inputs + 1)]
    optimize_hyperparameters!(emulator; kernbounds = kb)

    JLD2.save(emulator_path, "emulator", emulator)
    @info "Emulator saved → $emulator_path"
end

# ── Emulator validation ───────────────────────────────────────────────────────
y_mean_train, _ = Emulators.predict(emulator, get_inputs(train_pairs); transform_to_real = true)
y_mean_test,  _ = Emulators.predict(emulator, get_inputs(test_pairs);  transform_to_real = true)

mse_train = mean((get_outputs(train_pairs) .- y_mean_train) .^ 2)
mse_test  = mean((get_outputs(test_pairs)  .- y_mean_test)  .^ 2)

println("\n╔═══════════════════════════════════════╗")
println("║        Emulator validation            ║")
println("╠═══════════════════════════════════════╣")
println("║  GP MSE (train): $(lpad(round(mse_train, sigdigits=4), 15))      ║")
println("║  GP MSE (test) : $(lpad(round(mse_test,  sigdigits=4), 15))      ║")
println("╚═══════════════════════════════════════╝")
println("NOTE: A test MSE much larger than the train MSE indicates overfitting.")
println("      Proper scientific validation (coverage tests, rank histograms)")
println("      should be performed before using this posterior in publications.")

# ── MCMC posterior sampling ───────────────────────────────────────────────────
init_sample     = EKP.get_u_mean(ekp, maximum(train_iterations))
obs_sample      = EKP.get_obs(ekp)   # observation vector (AMorAV required by MCMCWrapper)

mcmc = MCMCWrapper(
    pCNMHSampling(),
    obs_sample,
    prior,
    emulator;
    init_params = init_sample,
)

@info "Optimising MCMC step size (N=2000 test steps)…"
new_step = try
    optimize_stepsize(rng, mcmc; init_stepsize, N = 2000, discard_initial = 0, max_iter = 40)
catch e
    @warn "optimize_stepsize did not converge ($(e)). Falling back to init_stepsize=$init_stepsize."
    init_stepsize
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
