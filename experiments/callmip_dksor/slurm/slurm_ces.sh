#!/bin/bash
#SBATCH --job-name=callmip_ces
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=experiments/callmip_dksor/slurm/logs/slurm_ces_%j.out
#SBATCH --error=experiments/callmip_dksor/slurm/logs/slurm_ces_%j.err

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

julia --threads=8 --project=.buildkite \
      experiments/callmip_dksor/emulate_sample.jl

echo "Done at $(date)"
