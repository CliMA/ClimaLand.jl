#!/bin/bash
#SBATCH --job-name=dk_sor_calibration
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/calibration_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/calibration_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --partition=expansion

# ── Environment ───────────────────────────────────────────────────────────────
export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ── Navigate to repo root ──────────────────────────────────────────────────────
cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor EKI calibration ==="
echo "  Host: $(hostname)  Nodes: $(scontrol show hostnames | wc -l)"
echo "  Julia: $(julia --version)"
echo "  Procs: ${SLURM_NTASKS} (${SLURM_CPUS_PER_TASK} threads each)"
echo "  Start: $(date)"

# 33 ensemble members: run_calibration.jl uses WorkerBackend which spawns workers.
# SLURM_NTASKS=32 leaves 1 task for the driver; each worker gets 8 threads.
julia --project=experiments/callmip_phase1a_v2 \
      -p $(( SLURM_NTASKS - 1 )) \
      experiments/callmip_phase1a_v2/run_calibration.jl

EXIT_CODE=$?
echo "Exit code: $EXIT_CODE  at $(date)"
exit $EXIT_CODE
