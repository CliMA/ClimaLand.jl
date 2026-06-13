#!/bin/bash
#SBATCH --job-name=callmip_sims
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=experiments/callmip_dksor/slurm/logs/slurm_callmip_sims_%j.out
#SBATCH --error=experiments/callmip_dksor/slurm/logs/slurm_callmip_sims_%j.err

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module purge
module use /groups/esm/modules
module load climacommon

# Find repo root from submit dir
REPO_ROOT="$SLURM_SUBMIT_DIR"
while [[ ! -d "$REPO_ROOT/.buildkite" && "$REPO_ROOT" != "/" ]]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
cd "$REPO_ROOT" || { echo "ERROR: could not find repo root"; exit 1; }
echo "Repo root: $REPO_ROOT"
mkdir -p experiments/callmip_dksor/slurm/logs

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate()'

# 2 workers (prior + posterior) + 1 master = 3 tasks total
julia -p 2 --project=.buildkite \
      experiments/callmip_dksor/run_callmip_simulations.jl

echo "Done at $(date)"
