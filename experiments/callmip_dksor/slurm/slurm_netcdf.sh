#!/bin/bash
#SBATCH --job-name=callmip_nc
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=experiments/callmip_dksor/slurm/logs/slurm_netcdf_%j.out
#SBATCH --error=experiments/callmip_dksor/slurm/logs/slurm_netcdf_%j.err

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

julia --project=.buildkite \
      experiments/callmip_dksor/write_callmip_netcdf.jl

echo "Done at $(date)"
