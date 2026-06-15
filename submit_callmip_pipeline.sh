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
#   bash submit_callmip_pipeline.sh --watch
#
# Each job starts only if the previous one completed successfully (afterok).
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"   # ensure we're at the repo root

WATCH=0
if [[ "${1:-}" == "--watch" ]]; then
    WATCH=1
fi

ETA_CALIBRATION="20:00:00"
ETA_EMULATE="48:00:00"
ETA_CALLMIP="24:00:00"
ETA_POSTPROCESS="03:00:00"

echo "=== CalLMIP Phase 1a pipeline submission ==="
echo "Repo root: $PWD"
echo ""

echo "Running readiness test first..."
if ! bash ./test_callmip_pipeline.sh; then
    echo ""
    echo "Pipeline readiness test failed. Not submitting."
    exit 1
fi

echo ""
echo "Planned stage ETAs:"
echo "  [1/4] Calibration   ETA ${ETA_CALIBRATION}"
echo "  [2/4] Emulate+Sample ETA ${ETA_EMULATE}"
echo "  [3/4] CalLMIP sims  ETA ${ETA_CALLMIP}"
echo "  [4/4] Postprocess   ETA ${ETA_POSTPROCESS}"
echo ""

# ── Step 1: EKI Calibration ─────────────────────────────────────────────────
JOB1=$(sbatch --parsable \
    experiments/calibrate_dk_sor/slurm_calibration.sh)
echo "  [1/4] Calibration   submitted → job $JOB1  (ETA ${ETA_CALIBRATION})"

# ── Step 2: CES Emulate + Sample ────────────────────────────────────────────
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh)
echo "  [2/4] Emulate+Sample submitted → job $JOB2  (after $JOB1, ETA ${ETA_EMULATE})"

# ── Step 3: CalLMIP Forward Simulations ─────────────────────────────────────
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh)
echo "  [3/4] CalLMIP sims  submitted → job $JOB3  (after $JOB2, ETA ${ETA_CALLMIP})"

# ── Step 4: Postprocess + Write NetCDF ──────────────────────────────────────
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    experiments/callmip_uq_dk_sor/slurm_postprocess.sh)
echo "  [4/4] Postprocess   submitted → job $JOB4  (after $JOB3, ETA ${ETA_POSTPROCESS})"

echo ""
echo "Stage job map:"
echo "  calibration    $JOB1"
echo "  emulate_sample $JOB2"
echo "  callmip_sims   $JOB3"
echo "  postprocess    $JOB4"

echo ""
echo "All jobs queued. Monitor with:"
echo "  squeue -u $USER"
echo ""
echo "Expected outputs when complete:"
echo "  experiments/callmip_uq_dk_sor/callmip_output/"
echo "    ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Prior.nc"
echo "    ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Posterior.nc"

if [[ $WATCH -eq 1 ]]; then
    echo ""
    echo "Watching jobs and reporting stage status..."
    while true; do
        TS=$(date '+%Y-%m-%d %H:%M:%S')
        S1=$(squeue -h -j "$JOB1" -o '%T' | head -1 || true)
        S2=$(squeue -h -j "$JOB2" -o '%T' | head -1 || true)
        S3=$(squeue -h -j "$JOB3" -o '%T' | head -1 || true)
        S4=$(squeue -h -j "$JOB4" -o '%T' | head -1 || true)
        echo "[$TS] calibration=${S1:-DONE/UNKNOWN} emulate=${S2:-DONE/UNKNOWN} callmip=${S3:-DONE/UNKNOWN} postprocess=${S4:-DONE/UNKNOWN}"
        if [[ -z "$S1$S2$S3$S4" ]]; then
            echo "All jobs have left the queue. Check sacct or stage logs for final states."
            break
        fi
        sleep 60
    done
fi