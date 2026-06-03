#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

PASS_COUNT=0
FAIL_COUNT=0

pass() {
    echo "[PASS] $1"
    PASS_COUNT=$((PASS_COUNT + 1))
}

fail() {
    echo "[FAIL] $1"
    FAIL_COUNT=$((FAIL_COUNT + 1))
}

check_file() {
    local path="$1"
    if [[ -f "$path" ]]; then
        pass "file exists: $path"
    else
        fail "missing file: $path"
    fi
}

check_bash_syntax() {
    local path="$1"
    if bash -n "$path"; then
        pass "bash syntax: $path"
    else
        fail "bash syntax error: $path"
    fi
}

echo "=== CalLMIP Pipeline Readiness Test ==="
echo "Repo root: $PWD"
echo

echo "[1/6] Checking critical files"
check_file "submit_callmip_pipeline.sh"
check_file "experiments/calibrate_dk_sor/slurm_calibration.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_postprocess.sh"
check_file "experiments/calibrate_dk_sor/generate_observations.jl"
check_file "experiments/calibrate_dk_sor/run_calibration.jl"
check_file "experiments/callmip_uq_dk_sor/emulate_sample.jl"
check_file "experiments/callmip_uq_dk_sor/run_callmip_simulations.jl"
check_file "experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl"
check_file "DK_Sor/DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc"
check_file "DK_Sor/DK-Sor_1997-2014_FLUXNET2015_Met.nc"

echo
echo "[2/6] Checking script syntax"
check_bash_syntax "submit_callmip_pipeline.sh"
check_bash_syntax "experiments/calibrate_dk_sor/slurm_calibration.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_postprocess.sh"

echo
echo "[3/6] Checking observation artifact schema compatibility"
if grep -q "y_obs =" experiments/calibrate_dk_sor/generate_observations.jl && \
   grep -q "noise_cov =" experiments/calibrate_dk_sor/generate_observations.jl && \
   grep -q "obs_dates =" experiments/calibrate_dk_sor/generate_observations.jl; then
    pass "generate_observations writes calibration+CES keys"
else
    fail "generate_observations missing calibration+CES keys"
fi

if [[ -f "experiments/calibrate_dk_sor/observations.jld2" ]]; then
    echo "[WARN] observations.jld2 exists; stage 1 will regenerate it before calibration"
else
    echo "[WARN] observations.jld2 not present yet; stage 1 will generate it"
fi

echo
echo "[4/6] Checking stage wiring expectations"
if grep -q "analyze_posterior_ensemble.jl" experiments/callmip_uq_dk_sor/slurm_postprocess.sh; then
    pass "postprocess includes posterior analysis"
else
    fail "postprocess missing analyze_posterior_ensemble stage"
fi

if grep -q "DateTime(1996, 1, 1)" experiments/callmip_uq_dk_sor/callmip_model_interface.jl && \
   grep -q "OUTPUT_START_DATE = Date(1997, 1, 1)" experiments/callmip_uq_dk_sor/callmip_model_interface.jl; then
    pass "callmip spinup/output period configured (1996 spinup, 1997-2014 output)"
else
    fail "callmip period/spinup configuration mismatch"
fi

if grep -q "const DT          = Float64(900)" experiments/callmip_uq_dk_sor/run_callmip_simulations.jl && \
   grep -q "const DT         = Float64(900)" experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl; then
    pass "forward stages DT aligned to calibration"
else
    fail "DT mismatch across stages"
fi

echo
echo "[5/6] Checking SLURM dry-run parse (if available)"
if command -v sbatch >/dev/null 2>&1; then
    DRYRUN_OK=1
    for script in \
        experiments/calibrate_dk_sor/slurm_calibration.sh \
        experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh \
        experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh \
        experiments/callmip_uq_dk_sor/slurm_postprocess.sh; do
        if sbatch --test-only "$script" >/tmp/callmip_sbatch_test.log 2>&1; then
            pass "sbatch --test-only: $script"
        else
            DRYRUN_OK=0
            fail "sbatch --test-only failed: $script"
        fi
    done
    [[ $DRYRUN_OK -eq 1 ]] || echo "See /tmp/callmip_sbatch_test.log"
else
    fail "sbatch not found in PATH"
fi

echo
echo "[6/6] Summary"
echo "Pass: $PASS_COUNT"
echo "Fail: $FAIL_COUNT"

if [[ $FAIL_COUNT -eq 0 ]]; then
    echo "READY_TO_SUBMIT=YES"
    exit 0
else
    echo "READY_TO_SUBMIT=NO"
    exit 1
fi