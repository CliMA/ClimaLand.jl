# # Sample — ClimaLand CES UQ Pipeline (Step 3 of 3)
#
# Runs Markov chain Monte Carlo (MCMC) over the emulator trained by emulate.jl
# to obtain full posterior distributions for the calibrated ClimaLand parameters.
#
# The emulator is cheap to evaluate, so the MCMC can draw hundreds of thousands
# of samples without requiring any additional model runs. This yields marginal
# posterior distributions, credible intervals, and parameter correlations.
#
# Prerequisites
# -------------
# Run emulate.jl first. The emulator file must exist at:
#   <CALIBRATION_OUTPUT_DIR>/ces/emulators.jld2
#
# Usage (from the ClimaLand.jl root)
# -----------------------------------
#   julia --project=.buildkite experiments/calibration/sample.jl
#
# Environment variables
# ----------------------
#   CALIBRATION_OUTPUT_DIR  Must match what was used in emulate.jl
#                           (default: experiments/calibration/land_model)
#   CES_N_SAMPLES           Number of MCMC samples  (default: 100_000)
#   CES_INIT_STEPSIZE       Initial RWMH step size  (default: 0.1)
#   CES_DISCARD_INITIAL     Burn-in samples to drop (default: 2_000)

using LinearAlgebra
using Statistics
using Printf
using JLD2: JLD2
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions
using ClimaCalibrate: ClimaCalibrate
using ClimaLand: ClimaLand

using CalibrateEmulateSample: CalibrateEmulateSample
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.MarkovChainMonteCarlo

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
const OUTPUT_DIR = get(
    ENV,
    "CALIBRATION_OUTPUT_DIR",
    joinpath("experiments", "calibration", "land_model"),
)
ces_output_dir = joinpath(OUTPUT_DIR, "ces")
emulator_file = joinpath(ces_output_dir, "emulators.jld2")

isfile(emulator_file) ||
    error("Emulator file not found: $emulator_file\n" * "Run emulate.jl first.")

n_samples = parse(Int, get(ENV, "CES_N_SAMPLES", "100000"))
init_stepsize = parse(Float64, get(ENV, "CES_INIT_STEPSIZE", "0.1"))
discard_initial = parse(Int, get(ENV, "CES_DISCARD_INITIAL", "2000"))

# ------------------------------------------------------------------
# Load emulator and calibration metadata
# ------------------------------------------------------------------
@info "Loading emulator from: $emulator_file"
emulator = JLD2.load(emulator_file, "emulator")
prior = JLD2.load(emulator_file, "prior")
y_obs = JLD2.load(emulator_file, "y_obs")

# ------------------------------------------------------------------
# Load the EKP to initialise MCMC near the posterior mode
# ------------------------------------------------------------------
# Starting MCMC from the EKI posterior mean reduces burn-in significantly.
ekp = ClimaCalibrate.load_latest_ekp(OUTPUT_DIR)
isnothing(ekp) && error("Could not load EKP from $OUTPUT_DIR")
init_params = EKP.get_u_mean_final(ekp)
@info "MCMC initial point (unconstrained space): $init_params"

# ------------------------------------------------------------------
# Run MCMC
# ------------------------------------------------------------------
@info "Running MCMC: $n_samples samples (discarding first $discard_initial)..."
mcmc = MCMCWrapper(RWMHSampling(), y_obs, prior, emulator; init_params)
new_step = optimize_stepsize(mcmc; init_stepsize, N = 2000, discard_initial = 0)
@info "MCMC step size optimised: $new_step"
chain = MarkovChainMonteCarlo.sample(
    mcmc,
    n_samples;
    stepsize = new_step,
    discard_initial,
)

display(chain)

# ------------------------------------------------------------------
# Save posterior and report statistics
# ------------------------------------------------------------------
posterior_file = joinpath(ces_output_dir, "posterior.jld2")
posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)
JLD2.save_object(posterior_file, posterior)
@info "Posterior saved → $posterior_file"

# Transform from unconstrained (MCMC) space back to physical parameter space
constrained_posterior = Emulators.transform_unconstrained_to_constrained(
    prior,
    MarkovChainMonteCarlo.get_distribution(posterior),
)

param_names = ParameterDistributions.name(prior)
@info "--- Posterior Statistics (physical/constrained space) ---"
for pname in param_names
    samples = vec(constrained_posterior[pname])
    @info @sprintf(
        "  %-30s  mean = %10.4g   std = %10.4g   95%% CI = [%10.4g, %10.4g]",
        pname,
        mean(samples),
        std(samples),
        quantile(samples, 0.025),
        quantile(samples, 0.975),
    )
end

@info "CES pipeline complete. Posterior samples are in: $posterior_file"
