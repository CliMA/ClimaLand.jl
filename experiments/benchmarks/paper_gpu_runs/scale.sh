#!/bin/bash
#PBS -A UCIT0011
#PBS -N scale
#PBS -q main
#PBS -m n
#PBS -M your_email@uni.edu
#PBS -m bae
#PBS -l walltime=00:13:00
#PBS -l select=64:ncpus=64:mem=480GB:ngpus=4

# How to use:
# To run with 1 gpu, set ngpus on line 9 to 1, and select to 1. Then uncomment line 33
# To run with 2 gpus, set ngpus on line 9 to 2, and select to 1. Then uncomment line 34
# To run with 4+ gpus, set ngpus on line 9 to 4, and select to (number of desired gpus divided by 4). Then uncomment line 35
# and set -n to number of desired gpus


# Use scratch for temporary files to avoid space limits in /tmp
export TMPDIR=${SCRATCH}/temp
mkdir -p ${TMPDIR}

export CLIMACOMMS_DEVICE=CUDA
export CLIMACOMMS_CONTEXT=MPI

# This is specific to derecho
source /glade/u/apps/derecho/23.09/spack/opt/spack/lmod/8.7.24/gcc/7.5.0/c645/lmod/lmod/init/zsh

export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
module load climacommon

# Run code

# mpiexec -n 1 -ppn 1 --cpu-bind verbose,list:0:16 set_gpu_rank julia --threads=3 --project=../../../.buildkite experiments/benchmarks/paper_gpu_runs/paper_gpu_runs.jl
# mpiexec -n 2 -ppn 2 --cpu-bind verbose,list:0:16 set_gpu_rank julia --threads=3 --project=../../../.buildkite experiments/benchmarks/paper_gpu_runs/paper_gpu_runs.jl
# mpiexec -n 256 -ppn 4 --cpu-bind verbose,list:0:16:32:48 set_gpu_rank julia --threads=3 --project=../../../.buildkite experiments/benchmarks/paper_gpu_runs/strong_scaling_gpu_runs.jl
