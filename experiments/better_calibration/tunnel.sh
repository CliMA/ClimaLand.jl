#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="tunnel.txt"
#SBATCH --time=24:00:00
#SBATCH --ntasks=7
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4

./code tunnel
