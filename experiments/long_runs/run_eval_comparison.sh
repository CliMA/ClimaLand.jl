#!/bin/bash
# =============================================================================
# Evaluation: Compare main_uncal vs changes_uncal (2-year global LandModel run)
#
# This script runs the snowy_land_pmodel.jl global simulation twice:
#   1. on the 'main' branch with uncalibrated parameters  → main_uncal_<device>/
#   2. on the current branch with uncalibrated parameters  → changes_uncal_<device>/
#
# After both runs complete, it generates leaderboard RMSE plots for each run
# and runs a side-by-side comparison script.
#
# Usage:
#   # Run both sequentially (takes ~2x simulation time):
#   bash experiments/long_runs/run_eval_comparison.sh
#
#   # Run only baseline (main):
#   bash experiments/long_runs/run_eval_comparison.sh baseline
#
#   # Run only changes (current branch):
#   bash experiments/long_runs/run_eval_comparison.sh changes
#
#   # Run only the comparison (after both runs are done):
#   bash experiments/long_runs/run_eval_comparison.sh compare
# =============================================================================

set -euo pipefail

CLIMALAND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$CLIMALAND_DIR"

CHANGES_BRANCH=$(git branch --show-current)
JULIA=${JULIA:-"julia --project=."}
RUN_SCRIPT="experiments/long_runs/snowy_land_pmodel.jl"
COMPARE_SCRIPT="experiments/long_runs/compare_rmse.jl"

# Detect device (will be confirmed by Julia, but use gpu as default)
DEVICE=${DEVICE:-"gpu"}

BASELINE_ROOT="main_uncal_${DEVICE}"
CHANGES_ROOT="changes_uncal_${DEVICE}"

echo "=============================================="
echo "  Energy Flux Evaluation: main vs $CHANGES_BRANCH"
echo "=============================================="
echo "  Baseline output: ${BASELINE_ROOT}/"
echo "  Changes output:  ${CHANGES_ROOT}/"
echo "  Run script:      ${RUN_SCRIPT}"
echo "=============================================="

run_baseline() {
    echo ""
    echo "[1/2] Running BASELINE (main branch, uncalibrated)"
    echo "----------------------------------------------"

    # Save any uncommitted changes
    STASH_MSG="eval_comparison_$(date +%s)"
    if ! git diff --quiet || ! git diff --cached --quiet; then
        git stash push -m "$STASH_MSG"
        STASHED=true
    else
        STASHED=false
    fi

    git checkout main

    echo "Starting baseline simulation..."
    RUN_NAME=main_uncal UNCALIBRATED="" $JULIA "$RUN_SCRIPT"

    # Return to changes branch
    git checkout "$CHANGES_BRANCH"
    if [ "$STASHED" = true ]; then
        git stash pop
    fi

    echo ""
    echo "[1/2] BASELINE complete → ${BASELINE_ROOT}/"
}

run_changes() {
    echo ""
    echo "[2/2] Running CHANGES ($CHANGES_BRANCH branch, uncalibrated)"
    echo "----------------------------------------------"

    echo "Starting changes simulation..."
    RUN_NAME=changes_uncal UNCALIBRATED="" $JULIA "$RUN_SCRIPT"

    echo ""
    echo "[2/2] CHANGES complete → ${CHANGES_ROOT}/"
}

run_compare() {
    echo ""
    echo "[Compare] Computing RMSE differences"
    echo "----------------------------------------------"

    BASELINE_DIAG="${BASELINE_ROOT}/global_diagnostics/output_active"
    CHANGES_DIAG="${CHANGES_ROOT}/global_diagnostics/output_active"

    if [ ! -d "$BASELINE_DIAG" ]; then
        echo "ERROR: Baseline diagnostics not found at $BASELINE_DIAG"
        echo "       Run: bash $0 baseline"
        exit 1
    fi
    if [ ! -d "$CHANGES_DIAG" ]; then
        echo "ERROR: Changes diagnostics not found at $CHANGES_DIAG"
        echo "       Run: bash $0 changes"
        exit 1
    fi

    $JULIA "$COMPARE_SCRIPT" "$BASELINE_DIAG" "$CHANGES_DIAG"

    echo ""
    echo "=============================================="
    echo "  Evaluation complete!"
    echo ""
    echo "  Leaderboard plots:"
    echo "    Baseline: ${BASELINE_ROOT}/ERA5_global_rmse_and_bias_graphs.png"
    echo "    Changes:  ${CHANGES_ROOT}/ERA5_global_rmse_and_bias_graphs.png"
    echo "=============================================="
}

# --- Main ---
PHASE=${1:-all}

case $PHASE in
    baseline)
        run_baseline
        ;;
    changes)
        run_changes
        ;;
    compare)
        run_compare
        ;;
    all)
        run_baseline
        run_changes
        run_compare
        ;;
    *)
        echo "Usage: $0 {baseline|changes|compare|all}"
        exit 1
        ;;
esac
