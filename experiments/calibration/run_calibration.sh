#!/bin/bash

# See the README.md for more information on how to run a calibration

module load climacommon/2025_02_25

# Required on Derecho's Lustre scratch: HDF5 file locking is incompatible with
# the parallel filesystem and causes NC_EHDFERR (-101) when the orchestrator
# tries to read NetCDF files written by compute-job ensemble members.
export HDF5_USE_FILE_LOCKING=FALSE

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite/ experiments/calibration/generate_observations.jl
julia --project=.buildkite/ experiments/calibration/run_calibration.jl
