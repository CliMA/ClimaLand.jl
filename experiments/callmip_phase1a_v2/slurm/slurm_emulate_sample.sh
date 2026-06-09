#!/bin/bash
#SBATCH --job-name=dk_sor_emulate_sample
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/emulate_sample_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/emulate_sample_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=expansion

export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor emulate + sample ==="
echo "  Host: $(hostname)"
echo "  Start: $(date)"

julia --project=experiments/callmip_phase1a_v2 \
      experiments/callmip_phase1a_v2/emulate_sample.jl

EXIT_CODE=$?
echo "Exit code: $EXIT_CODE  at $(date)"
exit $EXIT_CODE
