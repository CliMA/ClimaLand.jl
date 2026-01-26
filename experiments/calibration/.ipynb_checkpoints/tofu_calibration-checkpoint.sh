#!/bin/bash
# A script to launch the global calibration on the tofu server

echo "Starting calibration run on tofu..."
echo "Start time: $(date)"

# --- CONFIGURATION ---
PROJECT_DIR="/home/egreich/Projects/Climaexplore/ClimaLand.jl"
RESULTS_DIR="/home/egreich/Projects/Climaexplore/ClimaLand.jl/results/run_$(date +%F)_tofu_v1"
mkdir -p "$RESULTS_DIR"

echo "Results will be saved in: $RESULTS_DIR"

# --- GPU SELECTION ---
#
# 1. Manually check GPU status with the `nvidia-smi` command.
# 2. Set the variable below to the ID of a free GPU (e.g., 0 or 1).
#
export CUDA_VISIBLE_DEVICES=0

# --- CLIMA ENVIRONMENT ---
# Tell ClimaComms to use the (now visible) NVIDIA GPU for computations.
export CLIMACOMMS_DEVICE="CUDA"

# --- EXECUTION ---
cd "$PROJECT_DIR"
echo "Running on GPU: $CUDA_VISIBLE_DEVICES"

# STEP 1: GENERATE THE OBSERVATION VECTOR
# This command creates the land_observation_vector.jld2 file.
echo "Step 1: Generating observation vector..."
nice julia --project experiments/calibration/generate_observations.jl && \
\
# STEP 2: RUN THE CALIBRATION
# The '&&' ensures this step only runs if Step 1 was successful.
echo "Step 2: Starting calibration using the new observation vector..."
nice julia --project experiments/calibration/run_uspac_calibration.jl --output-dir "$RESULTS_DIR"

echo "-----------------------------------"
echo "Calibration script finished."
echo "End time: $(date)"