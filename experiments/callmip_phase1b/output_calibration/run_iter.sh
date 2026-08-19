#!/bin/bash
# Detached calibration runner. CALLMIP1B_N_ITER controls how far it goes;
# it resumes from ekp_checkpoint.jld2 automatically.
L=/home/renatob/ClimaLand.jl/experiments/callmip_phase1b/output_calibration
export CALLMIP1B_N_ITER=${CALLMIP1B_N_ITER:-1}
julia -p 23 --project=/home/renatob/ClimaLand.jl/.buildkite \
  /home/renatob/ClimaLand.jl/experiments/callmip_phase1b/run_calibration_phase1b.jl \
  > $L/calibration_iter.log 2>&1
echo "exit=$?" >> $L/calibration_iter.log
touch $L/ITER_DONE
