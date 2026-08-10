#!/bin/bash
#SBATCH --job-name=calib_fluxnet
#SBATCH --output=experiments/integrated/generic_site/calibration_out/slurm_logs/%x_%A_%a.out
#SBATCH --error=experiments/integrated/generic_site/calibration_out/slurm_logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
# Set --array at submit time so it matches the actual site count, e.g.
#   sbatch --array=0-199%50 experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch
# `%50` caps concurrent tasks at 50; tune to your fair-share.
#
# One array task = one fluxnet site = one full TransformUnscented calibration
# (5 ensemble members × N_ITERS iterations) using ClimaCalibrate.JuliaBackend.
# Forward models run serially inside each task — no MPI, no GPU. Single-column
# simulations are too small to benefit from a GPU.
#
# Usage:
#   module load climacommon/2024_10_09          # central HPC
#   # Generate the site list once on the login node:
#   julia --project=.buildkite experiments/integrated/generic_site/list_fluxnet_sites.jl experiments/integrated/generic_site/calibration_out/site_ids.txt
#   N=$(wc -l < experiments/integrated/generic_site/calibration_out/site_ids.txt)
#   sbatch --array=0-$((N-1))%50 experiments/integrated/generic_site/calibrate_all_fluxnet.sbatch

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

SITE_LIST="experiments/integrated/generic_site/calibration_out/site_ids.txt"
mkdir -p experiments/integrated/generic_site/calibration_out/slurm_logs

if [[ ! -s "$SITE_LIST" ]]; then
    echo "Site list $SITE_LIST is missing or empty. Generate it first with:"
    echo "  julia --project=.buildkite experiments/integrated/generic_site/list_fluxnet_sites.jl > $SITE_LIST"
    exit 1
fi

# Slurm uses 1-based line numbers; arrays are 0-based, so add 1.
LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))
SITE_ID=$(sed -n "${LINE_NO}p" "$SITE_LIST")

if [[ -z "$SITE_ID" ]]; then
    echo "No site ID at line $LINE_NO of $SITE_LIST"
    exit 1
fi

echo "Array task $SLURM_ARRAY_TASK_ID -> site $SITE_ID"

# HDF5 file locking is incompatible with Lustre; matches run_calibration.sh.
export HDF5_USE_FILE_LOCKING=FALSE

# Match the AU-Cum buildkite settings; override per-submission with --export.
# Window length comes from the CALIBRATION_CONFIG file (er.jl: 365 days);
# only N_ITERS is honored as an env override in calibrate_site.jl.
export SITE_ID
export N_ITERS="${N_ITERS:-5}"
export LOCAL="false"

julia --color=yes --project=.buildkite \
    experiments/integrated/generic_site/calibrate_site.jl
