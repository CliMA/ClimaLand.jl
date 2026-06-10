# # Emulate — ClimaLand CES UQ Pipeline (Step 2 of 3)
#
# Builds a Gaussian Process emulator from the EKI training data produced by
# the calibration step (run_calibration.jl / run_calibration.sh).
#
# The calibration pipeline (ClimaCalibrate) already handles Step 1: it runs the
# full ClimaLand ensemble for each EKI iteration and stores the EKP checkpoint
# to disk as `<output_dir>/iteration_XXX/eki_file.jld2`.
#
# This script:
#   1. Loads the final EKP from those checkpoints.
#   2. Extracts all (parameter, G(parameter)) pairs accumulated over iterations.
#   3. Trains a Gaussian Process emulator on those pairs.
#   4. Saves the emulator, prior, and observation metadata for the sampling step.
#
# Prerequisites
# -------------
# A completed calibration run. The output directory must contain EKP checkpoint
# files written by ClimaCalibrate. Run calibration first with:
#   bash experiments/calibration/run_calibration.sh <output_dir>
#
# Usage (from the ClimaLand.jl root)
# -----------------------------------
# Set CALIBRATION_CONFIG and CALIBRATION_OUTPUT_DIR to match what was used in
# the calibration run, then:
#
#   julia --project=.buildkite experiments/calibration/emulate.jl
#
# Environment variables
# ----------------------
#   CALIBRATION_CONFIG      Config file in experiments/calibration/configs/
#                           (default: energy_fluxes.jl)
#   CALIBRATION_OUTPUT_DIR  Directory written by run_calibration.sh
#                           (default: experiments/calibration/land_model)

using Dates
using LinearAlgebra
using Statistics
import Random
import JLD2
import EnsembleKalmanProcesses as EKP
import ClimaCalibrate
import ClimaLand

import CalibrateEmulateSample
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.Utilities

include(joinpath(pkgdir(ClimaLand), "experiments", "calibration", "api.jl"))

# ------------------------------------------------------------------
# Configuration — mirrors what run_calibration.jl uses
# ------------------------------------------------------------------
const OUTPUT_DIR = get(
    ENV,
    "CALIBRATION_OUTPUT_DIR",
    joinpath("experiments", "calibration", "land_model"),
)

config_file = get(ENV, "CALIBRATION_CONFIG", "energy_fluxes.jl")
include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "calibration",
        "configs",
        config_file,
    ),
)

# get_calibration_prior() and CALIBRATE_CONFIG are now defined by the config file.
prior = get_calibration_prior()

# ------------------------------------------------------------------
# Load the completed EKP from ClimaCalibrate checkpoints
# ------------------------------------------------------------------
@info "Loading latest EKP from: $OUTPUT_DIR"
ekp = ClimaCalibrate.load_latest_ekp(OUTPUT_DIR)
isnothing(ekp) && error(
    "No completed EKP found in $OUTPUT_DIR. " *
    "Run the calibration first with: bash experiments/calibration/run_calibration.sh",
)

n_iterations = length(ekp.g)
@info "EKP loaded: $(EKP.get_N_ens(ekp)) members × $n_iterations iterations = " *
      "$(EKP.get_N_ens(ekp) * n_iterations) training samples"

# ------------------------------------------------------------------
# Extract (parameter, G(parameter)) training pairs
# ------------------------------------------------------------------
# get_training_points returns a PairedDataContainer in the unconstrained parameter space.
input_output_pairs = Utilities.get_training_points(ekp, n_iterations)

# ------------------------------------------------------------------
# Extract a representative observation vector and noise covariance
# for use in the MCMC sampling step
# ------------------------------------------------------------------
# ClimaLand calibration uses EKP.ObservationSeries (one Observation per
# minibatch window). We use the first seasonal window as the MCMC target.
# For a multi-window posterior, one could concatenate all observations.
obs_series = EKP.get_observation_series(ekp)
first_obs = first(EKP.get_observations(obs_series))
y_obs = EKP.get_obs(first_obs)
obs_noise_cov = EKP.get_obs_noise_cov(first_obs)

@info "Observation vector length: $(length(y_obs))"
@info "Observation noise covariance type: $(typeof(obs_noise_cov))"

# ------------------------------------------------------------------
# Build the Gaussian Process emulator
# ------------------------------------------------------------------
@info "Training Gaussian Process emulator (this may take several minutes for large observation vectors)..."
gppackage = GPJL()
gauss_proc = GaussianProcess(gppackage; noise_learn = false)
emulator = Emulator(gauss_proc, input_output_pairs; obs_noise_cov)
optimize_hyperparameters!(emulator)
@info "Emulator trained."

# ------------------------------------------------------------------
# Save outputs for the sampling step
# ------------------------------------------------------------------
ces_output_dir = joinpath(OUTPUT_DIR, "ces")
mkpath(ces_output_dir)
emulator_file = joinpath(ces_output_dir, "emulators.jld2")

JLD2.jldsave(emulator_file; emulator, prior, y_obs, obs_noise_cov, n_iterations)
@info "Emulator saved → $emulator_file"
@info "Next step: julia --project=.buildkite experiments/calibration/sample.jl"
