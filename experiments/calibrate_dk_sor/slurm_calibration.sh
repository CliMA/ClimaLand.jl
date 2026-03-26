#!/bin/bash
#SBATCH --job-name=dk_sor_calibration
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --time=20:00:00
#SBATCH --ntasks=29
#SBATCH --cpus-per-task=1
#SBATCH --output=dk_sor_calibration_%j.out
#SBATCH --error=dk_sor_calibration_%j.err

# CPU-only calibration (no GPUs needed for single-column runs)
export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module use /groups/esm/modules
module load climacommon

julia --project=.buildkite -e 'using Pkg; Pkg.update(); Pkg.instantiate(; verbose=true)'
julia --project=.buildkite experiments/calibrate_dk_sor/generate_observations.jl
julia --project=.buildkite experiments/calibrate_dk_sor/run_calibration.jl
