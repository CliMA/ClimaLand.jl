#!/bin/bash
#SBATCH --job-name=dk_sor_postprocess
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/postprocess_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/postprocess_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=expansion

export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor postprocessing ==="
echo "  Host: $(hostname)"
echo "  Start: $(date)"

EXP=experiments/callmip_phase1a_v2

echo "--- analyze_posterior_ensemble ---"
julia --project=$EXP $EXP/analyze_posterior_ensemble.jl || echo "WARNING: posterior ensemble analysis failed"

echo "--- evaluate_calibration ---"
julia --project=$EXP $EXP/evaluate_calibration.jl || echo "WARNING: evaluate_calibration failed"

echo "--- write_callmip_netcdf ---"
julia --project=$EXP $EXP/write_callmip_netcdf.jl

echo "--- plot_calibration_check ---"
julia --project=$EXP $EXP/plot_calibration_check.jl || echo "WARNING: calibration plots failed"

echo "--- plot_posterior_distributions ---"
julia --project=$EXP $EXP/plot_posterior_distributions.jl || echo "WARNING: posterior plots failed"

EXIT_CODE=$?
echo "Postprocess done at $(date)"
exit $EXIT_CODE
