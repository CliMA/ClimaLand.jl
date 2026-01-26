#!/bin/bash
#SBATCH --job-name=smap_obs_gen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=smap_obs_gen_%j.out
#SBATCH --error=smap_obs_gen_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=egreich@caltech.edu

# Load modules
module load julia

# Set environment variables
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Set SMAP data path
export SMAP_DATA_PATH="/resnick/scratch/egreich/SMAP/SMAP_L3_SM_P_E"

# Print job info
echo "======================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Node: ${SLURM_NODELIST}"
echo "Start Time: $(date)"
echo "======================================"
echo ""

# Print environment info
echo "Julia version:"
julia --version
echo ""
echo "SMAP data path: ${SMAP_DATA_PATH}"
echo "Julia threads: ${JULIA_NUM_THREADS}"
echo ""

# Navigate to ClimaLand directory
cd /resnick/scratch/egreich/ClimaLand.jl

# Clean up corrupted precompilation cache
echo "Cleaning Julia precompilation cache..."
rm -rf ~/.julia/compiled/v1.11/ClimaLand
rm -rf /scratch/egreich/.julia/compiled/v1.11/ClimaLand 2>/dev/null || true
echo ""

# Activate the main ClimaLand environment and add missing packages
echo "Setting up environment..."
julia --project=. -e '
using Pkg

# Add packages needed for calibration if not already present
packages = ["JLD2", "NCDatasets", "Glob", "Interpolations"]
for pkg in packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Adding package: $pkg")
        Pkg.add(pkg)
    end
end

println("Running instantiate...")
Pkg.instantiate(; verbose=true)

println("Forcing precompilation...")
Pkg.precompile()

println("Environment setup complete!")
'
echo ""

# Run the observation generation script
echo "======================================"
echo "Starting SMAP observation generation"
echo "======================================"
echo ""

julia --project=. experiments/calibration/generate_smap_observations.jl

# Check if output file was created
if [ -f "experiments/calibration/land_observation_vector_SMAP.jld2" ]; then
    echo ""
    echo "======================================"
    echo "SUCCESS: Observation vector created"
    echo "======================================"
    
    # Print file info
    ls -lh experiments/calibration/land_observation_vector_SMAP.jld2
    
    # Print summary statistics using Julia
    julia --project=. -e '
    using JLD2
    using Statistics
    
    data = JLD2.load("experiments/calibration/land_observation_vector_SMAP.jld2")
    
    println("\n========== SMAP Observation Vector Summary ==========")
    println("Number of observations: ", length(data["observation_vector"]))
    println("Number of valid indices: ", length(data["valid_indices"]))
    println("Grid elements: ", data["nelements"])
    println("Date ranges: ", data["sample_date_ranges"])
    println("Quality flag threshold: ", data["quality_flag_threshold"])
    println("Mask threshold: ", data["mask_threshold"])
    println("\nObservation statistics:")
    obs = data["observation_vector"]
    println("  Min: ", minimum(obs))
    println("  Max: ", maximum(obs))
    println("  Mean: ", mean(obs))
    println("  Median: ", median(obs))
    println("  Std: ", std(obs))
    println("\nMetadata:")
    for (i, meta) in enumerate(data["metadata"])
        println("Period ", i, ":")
        for (k, v) in meta
            println("  ", k, ": ", v)
        end
    end
    println("=====================================================")
    '
else
    echo ""
    echo "======================================"
    echo "ERROR: Observation vector NOT created"
    echo "======================================"
    exit 1
fi

echo ""
echo "======================================"
echo "End Time: $(date)"
echo "======================================"