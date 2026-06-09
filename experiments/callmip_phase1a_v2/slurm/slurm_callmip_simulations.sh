#!/bin/bash
#SBATCH --job-name=dk_sor_callmip_sims
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/callmip_sims_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/callmip_sims_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=expansion

export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor CalLMIP prior + posterior simulations ==="
echo "  Host: $(hostname)"
echo "  Start: $(date)"

julia --project=experiments/callmip_phase1a_v2 \
      experiments/callmip_phase1a_v2/run_callmip_simulations.jl

EXIT_CODE=$?
echo "Exit code: $EXIT_CODE  at $(date)"
exit $EXIT_CODE
