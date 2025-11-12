#!/bin/bash
#SBATCH --job-name=slurm_calibration
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=5:00:00
#SBATCH --ntasks=17
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu

module load julia

export MODULEPATH="/groups/esm/modules:$MODULEPATH"

# Set environment variables for CliMA
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"

# Build and run the Julia code
module load climacommon
julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite/ experiments/calibration/generate_observations.jl
julia --project experiments/calibration/generate_smap_observations.jl
julia --project experiments/calibration/generate_joint_observations.jl
julia --project=.buildkite/ experiments/calibration/run_uspac_calibration.jl
