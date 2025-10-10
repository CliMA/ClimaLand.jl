#!/bin/bash
#PBS -A UCIT0011
#PBS -N zzz_cl_scale
#PBS -q main
#PBS -m n
#PBS -M treddy@caltech.edu
#PBS -m bae
#PBS -l walltime=00:13:00
#PBS -l select=64:ncpus=64:mem=480GB:ngpus=4

# Use scratch for temporary files to avoid space limits in /tmp
export TMPDIR=${SCRATCH}/temp
mkdir -p ${TMPDIR}

export CLIMACOMMS_DEVICE=CUDA
export CLIMACOMMS_CONTEXT=MPI

# If you are using zsh as default shell
source /glade/u/apps/derecho/23.09/spack/opt/spack/lmod/8.7.24/gcc/7.5.0/c645/lmod/lmod/init/zsh

export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
module load climacommon

# Run code

# mpiexec -n 2 -ppn 2 --cpu-bind verbose,list:0:16 set_gpu_rank julia --threads=3 --project=.buildkite experiments/benchmarks/paper_gpu_runs/paper_gpu_runs.jl
mpiexec -n 256 -ppn 4 --cpu-bind verbose,list:0:16:32:48 set_gpu_rank julia --threads=3 --project=.buildkite experiments/benchmarks/paper_gpu_runs/paper_gpu_runs.jl
# mpiexec -n 1 -ppn 1 --cpu-bind verbose,list:0:16 set_gpu_rank julia --threads=3 --project=.buildkite experiments/benchmarks/paper_gpu_runs/paper_gpu_runs.jl
