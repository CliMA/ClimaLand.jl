"""
Extract calibrated parameter estimates from EKI output

This script loads the eki_file.jld2 from the calibration run and extracts
the final parameter estimates with their uncertainties.
"""

using JLD2
using Statistics
using Printf

# Get the ClimaLand root directory
climaland_root = dirname(dirname(@__DIR__))

# Try iteration 3 first, fall back to iteration 2 if needed
# Check both the new location and old backup location
output_dir = joinpath(climaland_root, "experiments", "calibration", "land_model")
old_output_dir = joinpath(climaland_root, "experiments", "calibration", "land_model_wrong_P50_kx_coord")

iteration = 3
eki_path = joinpath(output_dir, "iteration_$(string(iteration, pad=3))", "eki_file.jld2")

if !isfile(eki_path) && isdir(old_output_dir)
    # Try old directory
    eki_path = joinpath(old_output_dir, "iteration_$(string(iteration, pad=3))", "eki_file.jld2")
end

if !isfile(eki_path)
    iteration = 2
    eki_path = joinpath(output_dir, "iteration_$(string(iteration, pad=3))", "eki_file.jld2")
    if !isfile(eki_path) && isdir(old_output_dir)
        eki_path = joinpath(old_output_dir, "iteration_$(string(iteration, pad=3))", "eki_file.jld2")
    end
end

if !isfile(eki_path)
    error("Could not find eki_file.jld2 in iteration 002 or 003 (checked both land_model and land_model_wrong_P50_kx_coord)")
end

println("Loading EKI file from iteration $iteration:")
println("  $eki_path")

# Load the EKI file
ekp_data = JLD2.load(eki_path, "single_stored_object")

# Extract parameter ensemble (each column is a parameter, each row is an ensemble member)
# The last iteration contains the final parameter estimates
u_container = ekp_data.u[end]  # Final ensemble DataContainer
u = u_container.data  # Extract the actual array (N_params × N_ensemble)

println("\nEnsemble size: $(size(u))")
println("  Parameters: $(size(u, 1))")
println("  Ensemble members: $(size(u, 2))")

# Parameter names (from run_uspac_calibration.jl)
param_names = ["βkx_base", "βkx_coord", "βψx50_base", "βψx50_slope", "βΠR_base", "βΠR_slope"]

# Compute mean and std for each parameter across the ensemble
param_means = vec(mean(u, dims=2))
param_stds = vec(std(u, dims=2))

println("\n" * "="^70)
println("CALIBRATED PARAMETER ESTIMATES (Iteration $iteration)")
println("="^70)
println("Parameter         Mean        Std Dev     95% CI")
println("-"^70)

for i in 1:length(param_names)
    ci95 = 1.96 * param_stds[i]
    @printf("%-15s %10.4f  %10.4f  [%7.4f, %7.4f]\n", 
            param_names[i], param_means[i], param_stds[i],
            param_means[i] - ci95, param_means[i] + ci95)
end

println("="^70)

# Generate Julia code snippet for calculate_traits_from_parameters.jl
println("\n" * "="^70)
println("CODE SNIPPET FOR calculate_traits_from_parameters.jl:")
println("="^70)
println("# Load parameter estimates from calibration iteration $iteration")
println("parameters = (")
for i in 1:length(param_names)
    comment = if param_names[i] == "βkx_base"
        "# Baseline log conductivity"
    elseif param_names[i] == "βkx_coord"
        "# kx-P50 coordination (safety-efficiency trade-off)"
    elseif param_names[i] == "βψx50_base"
        "# P50 baseline"
    elseif param_names[i] == "βψx50_slope"
        "# P50 climate sensitivity"
    elseif param_names[i] == "βΠR_base"
        "# Stomatal strategy baseline (isohydric-anisohydric)"
    elseif param_names[i] == "βΠR_slope"
        "# Stomatal strategy climate sensitivity"
    else
        ""
    end
    
    # Format with proper spacing
    param_str = @sprintf("    %-15s = %7.4f", param_names[i], param_means[i])
    
    # Add comma if not last parameter
    if i < length(param_names)
        param_str *= ","
    end
    
    # Add comment
    println(param_str * "  " * comment)
end
println(")")
println("="^70)

# Save to file for easy copying
output_file = joinpath(output_dir, "calibrated_parameters_iter$(string(iteration, pad=3)).jl")
open(output_file, "w") do io
    println(io, "# Calibrated parameter estimates from iteration $iteration")
    println(io, "# Ensemble mean ± std dev")
    println(io, "")
    println(io, "parameters = (")
    for i in 1:length(param_names)
        comment = if param_names[i] == "βkx_base"
            "# Baseline log conductivity (± $(round(param_stds[i], digits=4)))"
        elseif param_names[i] == "βkx_coord"
            "# kx-P50 coordination (± $(round(param_stds[i], digits=4)))"
        elseif param_names[i] == "βψx50_base"
            "# P50 baseline (± $(round(param_stds[i], digits=4)))"
        elseif param_names[i] == "βψx50_slope"
            "# P50 climate sensitivity (± $(round(param_stds[i], digits=4)))"
        elseif param_names[i] == "βΠR_base"
            "# Stomatal strategy baseline (± $(round(param_stds[i], digits=4)))"
        elseif param_names[i] == "βΠR_slope"
            "# Stomatal strategy climate sensitivity (± $(round(param_stds[i], digits=4)))"
        else
            ""
        end
        
        param_str = @sprintf("    %-15s = %7.4f", param_names[i], param_means[i])
        
        if i < length(param_names)
            param_str *= ","
        end
        
        println(io, param_str * "  " * comment)
    end
    println(io, ")")
end

println("\n✓ Parameter estimates saved to: $output_file")
println("\nYou can now copy the 'parameters = (...)' block above")
println("into calculate_traits_from_parameters.jl (lines 39-47)")

# Also save the full ensemble for uncertainty propagation
ensemble_file = joinpath(output_dir, "parameter_ensemble_iter$(string(iteration, pad=3)).jld2")
JLD2.save(ensemble_file, 
    "parameter_ensemble", u,
    "parameter_names", param_names,
    "parameter_means", param_means,
    "parameter_stds", param_stds,
    "iteration", iteration
)

println("\n✓ Full parameter ensemble saved to: $ensemble_file")
println("  Use this for uncertainty propagation in trait calculations")
println("\nEnsemble contains:")
println("  - parameter_ensemble: $(size(u, 1)) params × $(size(u, 2)) ensemble members")
println("  - parameter_names: $(length(param_names)) names")
println("  - parameter_means: ensemble mean for each parameter")
println("  - parameter_stds: ensemble std dev for each parameter")
