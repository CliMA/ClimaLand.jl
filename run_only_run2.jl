# Load the full comparison script but skip Run 1
ENV["SKIP_RUN1"] = "true"

# Check if Run 1 is already complete
if isdir("comparison_default_height_gpu") && isdir("comparison_default_height_gpu/global_diagnostics")
    println("✓ Run 1 already complete. Skipping to Run 2...")
    
    # Include everything except the Run 1 execution
    # We'll manually run just Run 2
    include("experiments/long_runs/compare_canopy_heights.jl")
else
    println("ERROR: Run 1 not found. Please complete Run 1 first.")
    exit(1)
end
