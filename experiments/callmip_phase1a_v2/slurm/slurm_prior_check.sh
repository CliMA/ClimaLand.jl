#!/bin/bash
#SBATCH --job-name=dk_sor_prior_check
#SBATCH --output=experiments/callmip_phase1a_v2/slurm_logs/prior_check_%j.out
#SBATCH --error=experiments/callmip_phase1a_v2/slurm_logs/prior_check_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module purge
module use /groups/esm/modules
module load climacommon

REPO_ROOT="$SLURM_SUBMIT_DIR"
while [[ ! -d "$REPO_ROOT/.buildkite" && "$REPO_ROOT" != "/" ]]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
cd "$REPO_ROOT"
echo "Running from: $PWD"

EXP=experiments/callmip_phase1a_v2
mkdir -p ${EXP}/slurm_logs
julia --project=${EXP} -e 'using Pkg; Pkg.instantiate()'
julia --project=${EXP} ${EXP}/run_prior_check.jl
