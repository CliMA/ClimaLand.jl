#!/bin/bash
#SBATCH --job-name=global_run
#SBATCH --output=global_%j.log
#SBATCH --error=global_%j.err
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --partition=beta
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --account=esm

echo "Starting job on $(hostname)"
echo "Time: $(date)"

export MODULEPATH="/groups/esm/modules:$MODULEPATH"

module purge
module load julia
module load cuda
module load climacommon/2026_02_18

export JULIA_PROJECT=/resnick/groups/esm/xwu/ClimaLand.jl
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CLIMACOMMS_DEVICE=CUDA

# Optional: Log GPU information
echo "CLIMACOMMS_DEVICE=$CLIMACOMMS_DEVICE"
nvidia-smi

julia /resnick/groups/esm/xwu/ClimaLand.jl/experiments/long_runs/snowy_land_pmodel.jl 

echo "Finished at $(date)"

