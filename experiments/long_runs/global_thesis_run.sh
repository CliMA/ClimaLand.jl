#!/bin/bash
#SBATCH --job-name=acharbon_global_test
#SBATCH --output=global_test_output.txt
#SBATCH --error=global_test_error.txt
#SBATCH --time=48:00:00
#SBATCH --gpus=1
#SBATCH --gpus-per-task=1
#SBATCH --mail-user=acharbon@caltech.edu

export CLIMACOMMS_DEVICE="CUDA"
module purge
module load climacommon
julia --project=./ClimaLand.jl/.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=./ClimaLand.jl/.buildkite ./ClimaLand.jl/experiments/long_runs/global_thesis_runs.jl