#!/bin/bash
# =============================================================================
# CalLMIP Phase 1a — Full Pipeline Submission
#
# Submits all four SLURM jobs in dependency order:
#   1. calibration   — generate_observations.jl + run_calibration.jl (EKI)
#   2. emulate       — emulate_sample.jl (CES GP + MCMC)
#   3. callmip_sims  — run_callmip_simulations.jl (prior + posterior fwd runs)
#   4. postprocess   — evaluate_calibration.jl + write_callmip_netcdf.jl
#
# Usage (from repo root):
#   bash submit_callmip_pipeline.sh
#
# Each job starts only if the previous one completed successfully (afterok).
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"   # ensure we're at the repo root

echo "=== CalLMIP Phase 1a pipeline submission ==="
echo "Repo root: $PWD"
echo ""

# ── Step 1: EKI Calibration ─────────────────────────────────────────────────
JOB1=$(sbatch --parsable \
    experiments/calibrate_dk_sor/slurm_calibration.sh)
echo "  [1/4] Calibration   submitted → job $JOB1"

# ── Step 2: CES Emulate + Sample ────────────────────────────────────────────
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh)
echo "  [2/4] Emulate+Sample submitted → job $JOB2  (after $JOB1)"

# ── Step 3: CalLMIP Forward Simulations ─────────────────────────────────────
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh)
echo "  [3/4] CalLMIP sims  submitted → job $JOB3  (after $JOB2)"

# ── Step 4: Postprocess + Write NetCDF ──────────────────────────────────────
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    experiments/callmip_uq_dk_sor/slurm_postprocess.sh)
echo "  [4/4] Postprocess   submitted → job $JOB4  (after $JOB3)"

echo ""
echo "All jobs queued. Monitor with:"
echo "  squeue -u $USER"
echo ""
echo "Expected outputs when complete:"
echo "  experiments/callmip_uq_dk_sor/callmip_output/"
echo "    ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Prior.nc"
echo "    ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Posterior.nc"