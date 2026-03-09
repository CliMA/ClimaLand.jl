#!/bin/bash
#SBATCH --job-name=neon_sco2_calibration
#SBATCH --time=04:00:00
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=1
#SBATCH --output=neon_sco2_calibration_%j.out
#SBATCH --error=neon_sco2_calibration_%j.err

# CPU-only calibration (no GPUs needed for single-column runs)
export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

# Configure site (override defaults via env vars)
export NEON_SITE_ID="${NEON_SITE_ID:-NEON-srer}"
export NEON_SPINUP_DAYS="${NEON_SPINUP_DAYS:-20}"
export NEON_N_ITERATIONS="${NEON_N_ITERATIONS:-10}"

module load climacommon

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(; verbose=true)'
julia --project=.buildkite experiments/calibrate_neon/generate_observations.jl
julia --project=.buildkite experiments/calibrate_neon/run_calibration.jl
