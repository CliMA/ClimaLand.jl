#!/bin/bash
#SBATCH --job-name=dk_sor_emulate_sample
#SBATCH --output=experiments/callmip_phase1a_v2/slurm_logs/emulate_sample_%j.out
#SBATCH --error=experiments/callmip_phase1a_v2/slurm_logs/emulate_sample_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

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
julia --threads=${SLURM_CPUS_PER_TASK} --project=${EXP} ${EXP}/emulate_sample.jl
