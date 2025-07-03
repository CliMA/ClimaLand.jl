#!/bin/bash
#SBATCH --job-name=slurm_calibration
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=12:00:00
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-task=1

# Set environment variables for CliMA
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"

# Build and run the Julia code
module load climacommon
julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite/ experiments/better_calibration/calibrate_land.jl
