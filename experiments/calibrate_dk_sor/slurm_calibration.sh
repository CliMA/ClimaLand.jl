#!/bin/bash
#SBATCH --job-name=dk_sor_calibration
#SBATCH --time=04:00:00
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --output=dk_sor_calibration_%j.out
#SBATCH --error=dk_sor_calibration_%j.err

# CPU-only calibration (no GPUs needed for single-column runs)
export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module load climacommon

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(; verbose=true)'
julia --project=.buildkite experiments/calibrate_dk_sor/generate_observations.jl
julia --project=.buildkite experiments/calibrate_dk_sor/run_calibration.jl
