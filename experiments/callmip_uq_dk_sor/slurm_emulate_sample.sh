#!/bin/bash
#SBATCH --job-name=dk_sor_emulate_sample
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=dk_sor_emulate_sample_%j.out
#SBATCH --error=dk_sor_emulate_sample_%j.err

# CES emulate + MCMC sample step (CPU-only, single node).
# Run this after slurm_calibration.sh has completed successfully.
#
# The GP emulator and MCMC chain are CPU-bound and benefit from multiple
# threads (set via --cpus-per-task). Adjust memory if the observation
# vector is larger than ~7500 elements.

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module use /groups/esm/modules
module load climacommon

# Instantiate the UQ project environment (CalibrateEmulateSample etc.)
julia --project=experiments/callmip_uq_dk_sor \
      -e 'using Pkg; Pkg.instantiate()'

julia --threads=${SLURM_CPUS_PER_TASK} \
      --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/emulate_sample.jl
