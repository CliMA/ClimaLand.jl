#!/usr/bin/env bash
# Sweep the column-ensemble benchmark: the old single `Column` as a baseline, then the new
# `ColumnEnsemble` at powers of two. One Julia process per configuration, run strictly one at a
# time, so compilation, GC state, ClimaUtilities' global open-file cache, and (on GPU) the device
# itself are never shared between points.
#
# Every run appends a row to
# experiments/comparison/output/benchmark_column_ensemble_<cpu|gpu>.txt -- one file per device,
# so CPU and GPU results never mix. Plot them together with plot_column_ensemble_benchmark.jl.
#
# On the cluster, submit the whole sweep as ONE slurm job. It is sequential inside, which is what
# keeps a single process on the GPU at any moment:
#
#   sbatch --job-name=colens-cpu --output=colens-cpu-%j.log --time=08:00:00 --cpus-per-task=1 \
#          --mem=32G experiments/comparison/sweep_column_ensemble.sh
#
#   sbatch --job-name=colens-gpu --output=colens-gpu-%j.log --time=08:00:00 --gpus=1 \
#          --mem=32G --export=ALL,DEVICE=gpu experiments/comparison/sweep_column_ensemble.sh
#
# Do not submit one job per configuration: several would land on the same GPU and the timings
# would measure contention. Run interactively the same way with `DEVICE=gpu ./sweep...sh` inside
# an salloc.
#
#   ./experiments/comparison/sweep_column_ensemble.sh                # old 1, then new 1..32
#   ./experiments/comparison/sweep_column_ensemble.sh 1 2 4          # only these column counts
#   NSTEPS=64 ./experiments/comparison/sweep_column_ensemble.sh 1 2  # quick smoke sweep
#   SCALE_STEPS=1 ./experiments/comparison/sweep_column_ensemble.sh  # even wall time per point
#   DEVICE=gpu MAX_COLUMNS=64 ./experiments/comparison/sweep_column_ensemble.sh
#
# The default 960 steps is five simulated days at dt = 450, chosen so the timed window is long
# enough to average over GC pauses and to cross the forcing interval many times. Cost grows with
# both the step count and N, so the default sweep is dominated by its largest N; SCALE_STEPS=1
# trades averaging window for a flat wall time per point.
#
# A failing configuration is reported and the sweep continues, so one blowup at large N does not
# discard the points already gathered.

set -uo pipefail

SCRIPT=experiments/comparison/benchmark_column_ensemble.jl

# Hard-coded, and deliberately not derived from this script's own path: `sbatch` copies the batch
# script into a spool directory and runs it from there, so `dirname $0` resolves to /var/spool
# rather than the checkout. Override for another checkout with REPO_ROOT=...
REPO_ROOT=${REPO_ROOT:-/home/kphan2/worktree/ClimaLand.jl/test-multiple-columns}
cd "$REPO_ROOT" || exit 2

DEVICE=${DEVICE:-cpu}
case "$DEVICE" in
    cpu | gpu) ;;
    *)
        echo "DEVICE must be cpu or gpu, got '$DEVICE'" >&2
        exit 2
        ;;
esac

JULIA=${JULIA:-julia}
PROJECT=${PROJECT:-.buildkite}
NSTEPS=${NSTEPS:-960}
SITE=${SITE:-US-MOz}
MAX_COLUMNS=${MAX_COLUMNS:-32}
SKIP_OLD=${SKIP_OLD:-0}
# SCALE_STEPS=1 divides the step count by N (floored at MIN_STEPS) so every configuration takes
# roughly the same wall time instead of the largest N dominating the sweep. Off by default: a
# fixed step count gives every point the same averaging window.
SCALE_STEPS=${SCALE_STEPS:-0}
MIN_STEPS=${MIN_STEPS:-64}
# Optional prefix for each run, e.g. SRUN="srun -n1 --gpus=1" to make every configuration its own
# job step. Empty by default: inside an sbatch allocation Julia already sees the allocated device.
SRUN=${SRUN:-}
MODULES=${MODULES:-climacommon}

# `module` is a shell function, so a non-login shell (which is what sbatch gives you) may not have
# it until the init script is sourced.
if [ "${SKIP_MODULE_LOAD:-0}" != "1" ]; then
    if ! type module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
        # shellcheck disable=SC1091
        . /etc/profile.d/modules.sh
    fi
    if type module >/dev/null 2>&1; then
        # shellcheck disable=SC2086
        module load $MODULES
    else
        echo "note: no module command found; using julia from PATH" >&2
    fi
fi

# ClimaComms picks the device and context from the environment, and both domain constructors
# default to ClimaComms.device(), so this is the only place the device is selected.
if [ "$DEVICE" = "gpu" ]; then
    export CLIMACOMMS_DEVICE="CUDA"
    # `--gpus=1` already narrows the allocation, but an salloc holding several devices would
    # otherwise leave the choice to CUDA.jl. Keep exactly the first entry slurm handed us so a
    # sweep can never spread across GPUs, and so the log records which one it used.
    if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
        export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES%%,*}"
    fi
else
    # `CPU` resolves to CPUMultiThreaded when Julia has more than one thread, so the thread count
    # changes what is being measured; it is recorded in the results file's `nthreads` column.
    export CLIMACOMMS_DEVICE="CPU"
fi
export CLIMACOMMS_CONTEXT="SINGLETON"

# Column counts: the positional arguments when given, else powers of two up to MAX_COLUMNS.
if [ "$#" -gt 0 ]; then
    ncols_list=("$@")
else
    ncols_list=()
    n=1
    while [ "$n" -le "$MAX_COLUMNS" ]; do
        ncols_list+=("$n")
        n=$((n * 2))
    done
fi

failed=()
succeeded=0

steps_for() {
    local ncolumns=$1
    if [ "$SCALE_STEPS" = "1" ]; then
        local scaled=$((NSTEPS / ncolumns))
        if [ "$scaled" -lt "$MIN_STEPS" ]; then scaled=$MIN_STEPS; fi
        echo "$scaled"
    else
        echo "$NSTEPS"
    fi
}

run_case() {
    local mode=$1 ncolumns=$2
    local nsteps
    nsteps=$(steps_for "$ncolumns")
    echo
    echo "=== $mode  N=$ncolumns  nsteps=$nsteps  site=$SITE  device=$DEVICE ==="
    # shellcheck disable=SC2086
    if $SRUN "$JULIA" "--project=$PROJECT" "$SCRIPT" "$mode" "$ncolumns" "$nsteps" "$SITE"; then
        succeeded=$((succeeded + 1))
    else
        echo "!!! FAILED: $mode N=$ncolumns" >&2
        failed+=("$mode N=$ncolumns")
    fi
}

started=$(date +%s)
echo "Sweep: device=$DEVICE (CLIMACOMMS_DEVICE=$CLIMACOMMS_DEVICE) site=$SITE nsteps=$NSTEPS columns=${ncols_list[*]}"
echo "repo:  $REPO_ROOT"
echo "julia: $("$JULIA" --version 2>/dev/null || echo unavailable)"
if [ "$DEVICE" = "gpu" ]; then
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset, CUDA.jl picks device 0>}"
    # nvidia-smi ignores CUDA_VISIBLE_DEVICES, so query the allocated device explicitly rather
    # than listing every GPU on the node.
    if command -v nvidia-smi >/dev/null 2>&1; then
        nvidia-smi \
            ${CUDA_VISIBLE_DEVICES:+-i "$CUDA_VISIBLE_DEVICES"} \
            --query-gpu=index,name,memory.total --format=csv,noheader
    fi
fi

# The baseline first: a single `Column` on the old 0-D CSV forcing. `old` only accepts N = 1
# because that forcing carries no spatial information.
if [ "$SKIP_OLD" = "0" ]; then
    run_case old 1
fi

for n in "${ncols_list[@]}"; do
    run_case new "$n"
done

elapsed=$(( $(date +%s) - started ))
echo
echo "Sweep finished in ${elapsed}s"
if [ "$succeeded" -gt 0 ]; then
    echo "$succeeded configuration(s) appended to experiments/comparison/output/benchmark_column_ensemble_${DEVICE}.txt"
else
    echo "No configuration succeeded; nothing was written." >&2
fi

if [ "${#failed[@]}" -gt 0 ]; then
    echo "Failed configurations:" >&2
    for f in "${failed[@]}"; do echo "  $f" >&2; done
    exit 1
fi
