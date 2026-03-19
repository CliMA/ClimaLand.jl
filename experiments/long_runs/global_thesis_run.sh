#!/bin/bash
#SBATCH --job-name=acharbon_global_20yr_longrun
#SBATCH --output=output_20yr.txt
#SBATCH --error=error_20yr.txt
#SBATCH --time=12:00:00
#SBATCH --gpus-per-task=1
#SBATCH --mail-user=acharbon@caltech.edu

export CLIMACOMMS_DEVICE="CUDA"
module purge
module load climacommon
julia --project=./ClimaLand.jl/.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=./ClimaLand.jl/.buildkite ./ClimaLand/experiments/long_runs/global_thesis_runs.jl