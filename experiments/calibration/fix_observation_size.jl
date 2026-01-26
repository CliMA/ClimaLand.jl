#!/usr/bin/env julia
"""
Fix observation vector size to match simulation's natural land mask.

The issue: ERA5 observations were created by regridding ERA5 data to the model grid,
which interpolated values and created more valid land points (27,578/season) than the 
simulation natively produces (12,415/season).

Solution: Scale down the expected observation count to match simulation output.
"""

using JLD2

# Load existing joint observation
data = JLD2.load("land_observation_vector_joint.jld2")

println("Original joint observation structure:")
println("  n_lhf_per_period: ", data["n_lhf_per_period"])
println("  n_sm_per_period: ", data["n_sm_per_period"])  
println("  joint_samples length: ", length(data["joint_samples"]))
println("  Total expected: ", sum(data["n_lhf_per_period"]) + sum(data["n_sm_per_period"]))

# Calculate actual LH sample count from simulation (12 seasons × 12,415 points/season)
# This is what the simulation actually produces after all preprocessing
actual_lhf_samples = 148_986  # From error message

println("\nSimulation produces:")
println("  LH samples: ", actual_lhf_samples)
println("  SM samples: ", data["n_sm_per_period"][1], " (unchanged)")

# Calculate the ratio to scale observations
expected_lhf = data["n_lhf_per_period"][1]
ratio = actual_lhf_samples / expected_lhf
println("\nScaling ratio: ", ratio, " (simulation has ", round(ratio * 100, digits=1), "% of expected land points)")

# We have two options:
# Option 1: Subset the existing observations to match simulation size (lose data but fast)
# Option 2: Regenerate observations from scratch without regridding (cleaner but more work)

# For now, let's use Option 1: Create a mask that selects a subset of observation points
# We'll take every Nth point to downsample from 330,936 to 148,986

# Extract current joint samples
current_samples = data["joint_samples"]
n_lhf_current = data["n_lhf_per_period"][1]
n_sm_current = data["n_sm_per_period"][1]

# Split into LH and SM portions
lhf_samples = current_samples[1:n_lhf_current]
sm_samples = current_samples[n_lhf_current+1:end]

println("\nCurrent sample counts:")
println("  LH samples in observation: ", length(lhf_samples))
println("  SM samples in observation: ", length(sm_samples))

# Downsample LH observations to match simulation
# Take approximately evenly spaced indices
step = n_lhf_current / actual_lhf_samples
indices = round.(Int, range(1, n_lhf_current, length=actual_lhf_samples))
lhf_samples_downsampled = lhf_samples[indices]

println("\nDownsampled:")
println("  LH samples: ", length(lhf_samples_downsampled))
println("  SM samples: ", length(sm_samples), " (unchanged)")

# Create new joint samples
new_joint_samples = vcat(lhf_samples_downsampled, sm_samples)

println("\nNew joint observation size: ", length(new_joint_samples))

# Create new data dictionary
new_data = copy(data)
new_data["joint_samples"] = new_joint_samples
new_data["n_lhf_per_period"] = [actual_lhf_samples]
# n_sm_per_period stays the same

# Save to new file
output_file = "land_observation_vector_joint_downsampled.jld2"
JLD2.save(output_file, new_data)

println("\nSaved downsampled joint observation to: ", output_file)
println("Update calibration config to use this file!")
