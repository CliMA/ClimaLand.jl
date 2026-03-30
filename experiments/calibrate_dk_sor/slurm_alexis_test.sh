#!/bin/bash
#SBATCH --job-name=dk_sor_alexis_test
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=dk_sor_alexis_test_%j.out
#SBATCH --error=dk_sor_alexis_test_%j.err

# Single forward run with Alexis's 12-param posterior to verify calibration quality.
# Parameters already written to output/iteration_999/member_001/parameters.toml.
# No Pkg.update() — environment is already instantiated.

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module use /groups/esm/modules
module load climacommon

cd "$SLURM_SUBMIT_DIR"
# Walk up to repo root if submitted from subdir
REPO_ROOT="$SLURM_SUBMIT_DIR"
while [[ ! -d "$REPO_ROOT/.buildkite" && "$REPO_ROOT" != "/" ]]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
cd "$REPO_ROOT"
echo "Running from: $PWD"

julia --project=.buildkite experiments/calibrate_dk_sor/run_alexis_test.jl
