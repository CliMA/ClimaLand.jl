#!/bin/bash
# All calibration-site Prior deliverable runs, 20-way parallel, detached.
D=/home/renatob/ClimaLand.jl/experiments/callmip_phase1b/output_deliverables
SITES="CA-Qfo CH-Dav DE-Gri DE-Hai DE-Tha DK-Sor FI-Hyy FR-Pue IT-Lav IT-MBo IT-Noe NL-Loo RU-Fyo US-MMS US-NR1 US-SRM US-Ton US-Var US-Whs US-Wkg"
echo $SITES | tr ' ' '\n' | xargs -P 20 -I{} bash -c '
  julia --project=/home/renatob/ClimaLand.jl/.buildkite \
    /home/renatob/ClimaLand.jl/experiments/callmip_phase1b/deliverable_runs_phase1b.jl {} Prior \
    > '"$D"'/logs/{}_Prior.log 2>&1
  echo "exit=$?" >> '"$D"'/logs/{}_Prior.log'
touch $D/PRIORS_DONE
