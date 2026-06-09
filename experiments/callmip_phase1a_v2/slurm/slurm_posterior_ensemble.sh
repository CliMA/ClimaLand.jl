#!/bin/bash
#SBATCH --job-name=dk_sor_posterior_ens
#SBATCH --output=experiments/callmip_phase1a_v2/slurm_logs/posterior_ens_%j.out
#SBATCH --error=experiments/callmip_phase1a_v2/slurm_logs/posterior_ens_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=51
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=12:00:00

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

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
julia -p $(( SLURM_NTASKS - 1 )) --project=${EXP} ${EXP}/run_posterior_ensemble.jl
