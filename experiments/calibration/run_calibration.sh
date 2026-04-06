#!/bin/bash

# See the README.md for more information on how to run a calibration

module load climacommon/2025_02_25

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=.buildkite/ experiments/calibration/generate_observations.jl
julia --project=.buildkite/ experiments/calibration/run_calibration.jl
