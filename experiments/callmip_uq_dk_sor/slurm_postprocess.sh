#!/bin/bash
#SBATCH --job-name=callmip_postprocess
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=experiments/callmip_uq_dk_sor/logs/callmip_postprocess_%j.out
#SBATCH --error=experiments/callmip_uq_dk_sor/logs/callmip_postprocess_%j.err

set -euo pipefail

# Navigate to repo root (directory containing .buildkite) from submit location.
REPO_ROOT="$SLURM_SUBMIT_DIR"
while [[ ! -d "$REPO_ROOT/.buildkite" && "$REPO_ROOT" != "/" ]]; do
      REPO_ROOT="$(dirname "$REPO_ROOT")"
done
cd "$REPO_ROOT"

module use /groups/esm/modules
module load climacommon

echo "=== Step 1: Running evaluate_calibration.jl ==="
echo "=== Step 1: Running analyze_posterior_ensemble.jl ==="
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/analyze_posterior_ensemble.jl

echo "=== Step 2: Running evaluate_calibration.jl ==="
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/evaluate_calibration.jl

echo "=== Step 3: Writing CalLMIP NetCDF files ==="
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl

echo "=== Done ==="
