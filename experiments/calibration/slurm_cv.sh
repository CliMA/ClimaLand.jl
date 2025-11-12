#!/bin/bash
#SBATCH --job-name=smap_cv_processing
#SBATCH --output=smap_cv_%j.out
#SBATCH --error=smap_cv_%j.err
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

# SMAP CV Processing Script
# Processes SMAP L3_SM_P_E data to calculate coefficient of variation
# Usage: sbatch run_smap_cv.sh [START_YEAR] [STOP_YEAR]

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load modules
module load julia

# Set up environment variables
export SMAP_DATA_PATH="${SMAP_DATA_PATH:-/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E}"
export OUTPUT_DIR="/resnick/scratch/egreich/SMAP/output_cv"
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Parse command line arguments
START_YEAR=${1:-2015}
STOP_YEAR=${2:-2023}

echo "================================================================"
echo "SMAP Coefficient of Variation Processing"
echo "================================================================"
echo "Start Year:        $START_YEAR"
echo "Stop Year:         $STOP_YEAR"
echo "Data Directory:    $SMAP_DATA_PATH"
echo "Output Directory:  $OUTPUT_DIR"
echo "Julia Threads:     $JULIA_NUM_THREADS"
echo "================================================================"
echo ""

# Check if data directory exists
if [ ! -d "$SMAP_DATA_PATH" ]; then
    echo "ERROR: SMAP data directory not found: $SMAP_DATA_PATH"
    echo "Please set SMAP_DATA_PATH environment variable or create the directory"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the Julia script
echo "Starting Julia processing..."
julia --project=.buildkite \
      --threads=$JULIA_NUM_THREADS \
      experiments/calibration/process_smap_cv.jl \
      $START_YEAR \
      $STOP_YEAR \
      "$SMAP_DATA_PATH" \
      "$OUTPUT_DIR"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "================================================================"
    echo "Processing completed successfully!"
    echo "================================================================"
    echo "Output files saved to: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    ls -lh "$OUTPUT_DIR"/smap_cv_${START_YEAR}_${STOP_YEAR}.*
    echo "================================================================"
else
    echo ""
    echo "================================================================"
    echo "ERROR: Processing failed with exit code $EXIT_CODE"
    echo "================================================================"
    exit $EXIT_CODE
fi