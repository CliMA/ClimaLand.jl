#!/bin/bash
#SBATCH --job-name=dk_sor_calibration
#SBATCH --output=experiments/callmip_phase1a_v2/slurm_logs/calibration_%j.out
#SBATCH --error=experiments/callmip_phase1a_v2/slurm_logs/calibration_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=34
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=20:00:00

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
julia --project=${EXP} ${EXP}/generate_observations.jl
julia --project=${EXP} ${EXP}/run_calibration.jl
