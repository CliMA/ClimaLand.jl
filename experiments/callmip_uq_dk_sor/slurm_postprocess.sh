#!/bin/bash
#SBATCH --job-name=callmip_postprocess
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/callmip_postprocess_%j.out
#SBATCH --error=logs/callmip_postprocess_%j.err

set -euo pipefail

cd /resnick/home/renatob/ClimaLand.jl

module use /groups/esm/modules
module load climacommon

echo "=== Step 1: Running evaluate_calibration.jl ==="
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/evaluate_calibration.jl

echo "=== Step 2: Writing CalLMIP NetCDF files ==="
julia --project=experiments/callmip_uq_dk_sor \
      experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl

echo "=== Done ==="
