# Test of the full CES UQ pipeline for ClimaLand calibration
#
# This script validates the Calibrate → Emulate → Sample workflow without
# running the full ClimaLand model. It creates a synthetic calibration using
# a cheap toy forward model G(u) that mimics the real calibration structure,
# then tests the exact same ClimaCalibrate API calls made by emulate.jl and
# sample.jl.
#
# Usage (from the ClimaLand.jl root):
#   julia --project=.buildkite experiments/calibration/test/test_ces_pipeline.jl

using Test
using LinearAlgebra
using Statistics
using Printf
using Random: Random
using JLD2: JLD2
import EnsembleKalmanProcesses as EKP
using ClimaCalibrate: ClimaCalibrate
using ClimaLand: ClimaLand

using CalibrateEmulateSample: CalibrateEmulateSample
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.Utilities
import CalibrateEmulateSample.MarkovChainMonteCarlo

include(joinpath(pkgdir(ClimaLand), "experiments", "calibration", "api.jl"))

# ------------------------------------------------------------------
# Synthetic problem setup
# ------------------------------------------------------------------
# Mirrors the energy_fluxes.jl config: one parameter, emissivity_bare_soil.
# True value is 0.96 (as commented in run_calibration.jl).
# Toy forward model: G(emissivity) ≈ C * emissivity  (linearised LWU)
# This gives the same ClimaCalibrate data structures as the real model.

const TOY_TRUE_EMISSIVITY = 0.96
const TOY_C = 240.0  # ~ LWU in W/m2 at emissivity=1
const TOY_NOISE_VAR = 1.0  # observation noise variance

toy_G(emissivity) = [TOY_C * emissivity]

# Prior matches energy_fluxes.jl
prior = EKP.combine_distributions([
    EKP.constrained_gaussian("emissivity_bare_soil", 0.82, 0.12, 0.0, 2.0)
])

# Generate synthetic observation at the true parameter value
rng = Random.MersenneTwister(42)
y_obs_true = toy_G(TOY_TRUE_EMISSIVITY)
y_obs = y_obs_true .+ sqrt(TOY_NOISE_VAR) .* randn(rng, 1)
noise = TOY_NOISE_VAR * I(1)

@info "Synthetic observation: y_obs = $y_obs  (true = $y_obs_true)"

# ------------------------------------------------------------------
# Step 1 — Calibrate (synthetic EKI with ClimaCalibrate checkpoints)
# ------------------------------------------------------------------
output_dir = joinpath(tempdir(), "climaland_ces_test_$(getpid())")
mkpath(output_dir)
@info "Using synthetic output dir: $output_dir"

N_ens = 10
N_iterations = 8

# Build EKP directly (V6O6K's initialize only accepts an existing EKP)
rng_ekp = Random.MersenneTwister(42)
initial_ensemble = EKP.construct_initial_ensemble(rng_ekp, prior, N_ens)
ekp = EKP.EnsembleKalmanProcess(
    initial_ensemble, y_obs, noise, EKP.Inversion(); rng=rng_ekp
)

# Save the initial EKP to iteration 0 so load_latest_ekp can find it
# (load_latest_ekp searches from iter 0 upward)
mkpath(ClimaCalibrate.path_to_iteration(output_dir, 0))
ClimaCalibrate.save_eki_and_parameters(ekp, output_dir, 0, prior)

@testset "Calibration (synthetic EKI)" begin
    for i in 0:(N_iterations - 1)
        params_i = EKP.get_ϕ_final(prior, ekp)
        G_ens = hcat([toy_G(params_i[1, j]) for j in 1:N_ens]...)
        # update_ensemble! saves iteration i+1/eki_file.jld2
        ClimaCalibrate.update_ensemble!(ekp, G_ens, output_dir, i, prior)
    end

    final_mean = EKP.get_ϕ_mean_final(prior, ekp)[1]
    @info "EKI final mean: $final_mean  (true: $TOY_TRUE_EMISSIVITY)"
    @test abs(final_mean - TOY_TRUE_EMISSIVITY) < 0.05
end

# ------------------------------------------------------------------
# Step 2 — Emulate  (mirrors emulate.jl logic)
# ------------------------------------------------------------------
@testset "Emulate (GP emulator)" begin
    global ekp_loaded = ClimaCalibrate.load_latest_ekp(output_dir)
    @test !isnothing(ekp_loaded)

    # EKP may terminate early via DataMisfitController — use actual count
    global n_iter_loaded = length(ekp_loaded.g)
    @test n_iter_loaded >= 1

    # Extract (u, G(u)) pairs — same call as emulate.jl
    global iop = Utilities.get_training_points(ekp_loaded, n_iter_loaded)
    @test size(Utilities.get_inputs(iop), 2) == N_ens * n_iter_loaded

    # Extract y_obs and noise from EKP (same code as emulate.jl)
    obs_series = EKP.get_observation_series(ekp_loaded)
    first_obs = first(EKP.get_observations(obs_series))
    global y_obs_from_ekp = EKP.get_obs(first_obs)
    global obs_noise_cov = EKP.get_obs_noise_cov(first_obs)
    @test length(y_obs_from_ekp) == 1
    @test isapprox(y_obs_from_ekp[1], y_obs[1]; atol=1e-10)

    # Train GP emulator — same call as emulate.jl
    @info "Training GP emulator..."
    gppackage = GPJL()
    gauss_proc = GaussianProcess(gppackage; noise_learn=false)
    global emulator = Emulator(gauss_proc, iop; obs_noise_cov)
    optimize_hyperparameters!(emulator)
    @test emulator isa Emulator

    # Save for sampling step — same as emulate.jl
    ces_dir = joinpath(output_dir, "ces")
    mkpath(ces_dir)
    global emulator_file = joinpath(ces_dir, "emulators.jld2")
    JLD2.jldsave(
        emulator_file;
        emulator,
        prior,
        y_obs=y_obs_from_ekp,
        obs_noise_cov,
        n_iterations=n_iter_loaded,
    )
    @test isfile(emulator_file)
    @info "Emulator saved → $emulator_file"
end

# ------------------------------------------------------------------
# Step 3 — Sample  (mirrors sample.jl logic)
# ------------------------------------------------------------------
@testset "Sample (MCMC posterior)" begin
    # Load from file — same as sample.jl
    emulator_loaded = JLD2.load(emulator_file, "emulator")
    prior_loaded = JLD2.load(emulator_file, "prior")
    y_obs_loaded = JLD2.load(emulator_file, "y_obs")

    # Start MCMC from EKI posterior mean — same as sample.jl
    init_params = EKP.get_u_mean_final(ekp_loaded)
    @info "MCMC init (unconstrained): $init_params"

    # Run MCMC — same call as sample.jl
    @info "Running MCMC (20_000 samples)..."
    mcmc = MarkovChainMonteCarlo.MCMCWrapper(
        MarkovChainMonteCarlo.RWMHSampling(),
        y_obs_loaded,
        prior_loaded,
        emulator_loaded;
        init_params,
    )
    new_step = MarkovChainMonteCarlo.optimize_stepsize(
        mcmc; init_stepsize=0.1, N=2000, discard_initial=0
    )
    @info "MCMC step size: $new_step"
    chain = MarkovChainMonteCarlo.sample(
        mcmc, 20_000; stepsize=new_step, discard_initial=1_000
    )
    display(chain)

    # Save and extract posterior — same as sample.jl
    posterior_file = joinpath(output_dir, "ces", "posterior.jld2")
    posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)
    JLD2.save_object(posterior_file, posterior)
    @test isfile(posterior_file)

    constrained_posterior = Emulators.transform_unconstrained_to_constrained(
        prior_loaded, MarkovChainMonteCarlo.get_distribution(posterior)
    )

    samples_e = vec(constrained_posterior["emissivity_bare_soil"])
    post_mean = mean(samples_e)
    post_std = std(samples_e)
    ci_lo, ci_hi = quantile(samples_e, 0.025), quantile(samples_e, 0.975)

    @info @sprintf(
        "Posterior emissivity_bare_soil:  mean=%.4f  std=%.4f  95%% CI=[%.4f, %.4f]  true=%.4f",
        post_mean,
        post_std,
        ci_lo,
        ci_hi,
        TOY_TRUE_EMISSIVITY,
    )

    # Posterior mean should be within 2 prior std of truth
    @test abs(post_mean - TOY_TRUE_EMISSIVITY) < 0.15
    # Posterior should be narrower than prior (std < 0.12)
    @test post_std < 0.12
    # True value should be inside the 95% CI
    @test ci_lo <= TOY_TRUE_EMISSIVITY <= ci_hi
end

@info "All tests passed. Output in: $output_dir"
