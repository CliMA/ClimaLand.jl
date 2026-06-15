#!/bin/bash
#SBATCH --job-name=dk_sor_posterior_ensemble
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --time=12:00:00
#SBATCH --ntasks=51
#SBATCH --cpus-per-task=1
#SBATCH --output=dk_sor_posterior_ensemble_%j.out
#SBATCH --error=dk_sor_posterior_ensemble_%j.err

# Forward-model ensemble driven by the CES posterior (CPU-only).
# Run this after slurm_emulate_sample.sh has completed successfully.
#
# ntasks = N_SAMPLES + 1  (1 driver + N_SAMPLES workers).
# If you change N_SAMPLES in run_posterior_ensemble.jl,
# update --ntasks accordingly (N_SAMPLES + 1).

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module use /groups/esm/modules
module load climacommon

# Instantiate the .buildkite environment before launching workers
julia --project=.buildkite \
      -e 'using Pkg; Pkg.instantiate()'

julia --project=.buildkite \
      experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl
