#!/bin/bash
#SBATCH --job-name=dk_sor_posterior_ens
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/posterior_ens_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/posterior_ens_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=26
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --partition=expansion

export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor 50-member posterior ensemble ==="
echo "  Host: $(hostname)"
echo "  Procs: ${SLURM_NTASKS} (${SLURM_CPUS_PER_TASK} threads each)"
echo "  Start: $(date)"

# 25 workers for pmap (50 members / 2 rounds)
julia --project=experiments/callmip_phase1a_v2 \
      -p $(( SLURM_NTASKS - 1 )) \
      experiments/callmip_phase1a_v2/run_posterior_ensemble.jl

EXIT_CODE=$?
echo "Exit code: $EXIT_CODE  at $(date)"
exit $EXIT_CODE
