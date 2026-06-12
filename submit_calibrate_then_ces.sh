#!/bin/bash
# Submit DK-Sor calibration, then automatically queue the CES emulate_sample
# step to run upon successful completion.
#
# Usage (from repo root):
#   bash submit_calibrate_then_ces.sh

set -e
cd "$(dirname "$0")"

CAL_JOB=$(sbatch --parsable experiments/calibrate_dk_sor/slurm_calibration.sh)
echo "Calibration job submitted: $CAL_JOB"

CES_JOB=$(sbatch --parsable \
    --dependency=afterok:${CAL_JOB} \
    experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh)
echo "Emulate-sample job submitted: $CES_JOB (runs after job $CAL_JOB completes)"

echo ""
echo "Monitor with:"
echo "  squeue -j ${CAL_JOB},${CES_JOB}"
echo "  tail -f dk_sor_calibration_${CAL_JOB}.out"
