#!/bin/bash
#SBATCH --job-name=dk_sor_prior_check
#SBATCH --output=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/prior_check_%j.out
#SBATCH --error=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2/slurm_logs/prior_check_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=expansion

# ── Environment ───────────────────────────────────────────────────────────────
export JULIA_DEPOT_PATH=/central/scratch/renatob/julia_depot
export JULIA_PROJECT=/central/scratch/renatob/ClimaLand.jl/experiments/callmip_phase1a_v2
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ── Navigate to repo root ──────────────────────────────────────────────────────
# Never use $0 — SLURM copies the script to /var/spool before running it.
# Use SLURM_SUBMIT_DIR (always the directory where sbatch was called) or
# the fixed absolute path below.
cd /central/scratch/renatob/ClimaLand.jl

mkdir -p experiments/callmip_phase1a_v2/slurm_logs

echo "=== DK-Sor prior check ==="
echo "  Host: $(hostname)"
echo "  Julia: $(julia --version)"
echo "  Start: $(date)"

julia --project=experiments/callmip_phase1a_v2 \
      experiments/callmip_phase1a_v2/run_prior_check.jl

EXIT_CODE=$?
echo "Exit code: $EXIT_CODE  at $(date)"
exit $EXIT_CODE
