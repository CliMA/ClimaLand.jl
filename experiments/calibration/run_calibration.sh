#!/bin/bash

# See the README.md for more information on how to run a calibration.
#
# Usage: bash experiments/calibration/run_calibration.sh [OUTPUT_DIR]
#
# If OUTPUT_DIR is omitted, the default in run_calibration.jl is used.
#
# Before invoking this script, load climacommon to put Julia on PATH for the
# orchestrator process, e.g. `module load climacommon/2025_02_25` on Derecho.
# The version that forward-model PBS jobs load is pinned separately in
# run_calibration.jl via the `modules` kwarg on the backend constructor.

# Required on Derecho's Lustre scratch: HDF5 file locking is incompatible with
# the parallel filesystem and causes NC_EHDFERR (-101) when the orchestrator
# tries to read NetCDF files written by compute-job ensemble members.
export HDF5_USE_FILE_LOCKING=FALSE

OUTPUT_DIR="$1"

julia --project=.buildkite -e 'using Pkg; Pkg.update(); Pkg.instantiate()'
julia --project=.buildkite/ experiments/calibration/generate_observations.jl $OUTPUT_DIR
julia --project=.buildkite/ experiments/calibration/run_calibration.jl $OUTPUT_DIR
