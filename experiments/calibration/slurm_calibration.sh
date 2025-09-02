#!/bin/bash
#SBATCH --partition=a3
#SBATCH --job-name=slurm_calibration
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=12:00:00
#SBATCH --ntasks=11
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-task=1

# Set environment variables for CliMA
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"
export LONGER_RUN=""
# Build and run the Julia code
julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite experiments/calibration/run_calibration.jl
