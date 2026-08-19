#!/bin/bash
# Detached Phase 3 benchmark: GPU sweep then CPU sweep. Survives session exit.
L=/home/renatob/ClimaLand.jl/experiments/callmip_phase1b/benchmark_logs
export CLIMACOMMS_DEVICE=CUDA CUDA_VISIBLE_DEVICES=0
julia --project=/home/renatob/ClimaLand.jl/.buildkite /home/renatob/ClimaLand.jl/experiments/callmip_phase1b/benchmark_columns.jl > $L/bench_gpu.log 2>&1
echo "exit=$?" >> $L/bench_gpu.log
unset CLIMACOMMS_DEVICE CUDA_VISIBLE_DEVICES
julia --project=/home/renatob/ClimaLand.jl/.buildkite /home/renatob/ClimaLand.jl/experiments/callmip_phase1b/benchmark_columns.jl > $L/bench_cpu.log 2>&1
echo "exit=$?" >> $L/bench_cpu.log
touch $L/BENCHMARK_DONE
