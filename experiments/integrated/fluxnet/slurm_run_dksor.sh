#!/bin/bash
#SBATCH --job-name=fluxnet_DK-Sor
#SBATCH --output=experiments/integrated/fluxnet/slurm_logs/DK-Sor_%j.out
#SBATCH --error=experiments/integrated/fluxnet/slurm_logs/DK-Sor_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

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

mkdir -p experiments/integrated/fluxnet/slurm_logs

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate()'
julia --project=.buildkite experiments/integrated/fluxnet/DK-Sor/run_DK-Sor.jl
