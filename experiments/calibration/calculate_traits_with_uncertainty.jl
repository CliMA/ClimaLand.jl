"""Calculate ψx50, kx, and ΠR from calibrated parameter ensemble with uncertainty

This script propagates parameter uncertainty through the trait calculations
by computing traits for each ensemble member, then computing statistics.
"""

using JLD2
using Statistics
using ClimaCore
using ClimaLand
using ClimaComms

# Load Global Aridity Index data (P/ET0, dimensionless)
aridity_file = joinpath(@__DIR__, "aridity.jld2")
data = JLD2.load(aridity_file)
println("Available keys in aridity.jld2:")
println(keys(data))

# Extract CWD/aridity values
cwd_field = data["CWD_field"]
println("\nType of aridity field: $(typeof(cwd_field))")

# Get the spectral element space
field_space = ClimaCore.Fields.axes(cwd_field)

# Extract actual lat/lon coordinates from the Field
coords = ClimaCore.Fields.coordinate_field(field_space)
lats_field = coords.lat
lons_field = coords.long

# Access the underlying array directly from the Field
aridity = parent(cwd_field)  # Global Aridity Index (P/ET0)
lats = parent(lats_field)
lons = parent(lons_field)

println("\nData shape:")
println("  Aridity: $(size(aridity))")
println("  Total points: $(length(aridity))")

# Load parameter ensemble
climaland_root = dirname(dirname(@__DIR__))
output_dir = joinpath(climaland_root, "experiments", "calibration", "land_model")

# Try to find the ensemble file
global ensemble_file = nothing
for iter in [3, 2]
    candidate = joinpath(output_dir, "parameter_ensemble_iter$(string(iter, pad=3)).jld2")
    if isfile(candidate)
        global ensemble_file = candidate
        println("\nFound parameter ensemble: $ensemble_file")
        break
    end
end

if isnothing(ensemble_file)
    error("Could not find parameter_ensemble file")
end

# Load ensemble
ensemble_data = JLD2.load(ensemble_file)
param_ensemble = ensemble_data["parameter_ensemble"]  # (6 × N_ensemble)
param_names = ensemble_data["parameter_names"]
param_means = ensemble_data["parameter_means"]
param_stds = ensemble_data["parameter_stds"]

n_params, n_ensemble = size(param_ensemble)

println("\nParameter ensemble:")
println("  Parameters: $n_params")
println("  Ensemble members: $n_ensemble")
println("  Parameter names: $param_names")

# Normalize aridity
land_aridity = filter(!isnan, aridity)
println("\nAridity statistics:")
println("  Global Aridity Index (P/ET0) - min: $(minimum(land_aridity)), max: $(maximum(land_aridity)), mean: $(mean(land_aridity))")

# Normalization (matching stomatalconductance.jl)
aridity_norm = @. 1.0 - clamp(aridity / 2.0, 0.0, 1.0)

land_norm = filter(!isnan, aridity_norm)
println("  Normalized aridity - min: $(minimum(land_norm)), max: $(maximum(land_norm)), mean: $(mean(land_norm))")

# Function to compute traits for a single parameter vector
function compute_traits(params, aridity_norm)
    βkx_base, βkx_coord, βψx50_base, βψx50_slope, βΠR_base, βΠR_slope = params
    
    # ψx50 calculation (matching stomatalconductance.jl)
    ψx50_exponent = @. ifelse(isnan(aridity_norm), 0.0, 
                              clamp(βψx50_base + βψx50_slope * aridity_norm, -2.0, 4.0))
    ψx50 = @. ifelse(isnan(aridity_norm), NaN, -exp(ψx50_exponent))
    
    # kx calculation with CORRECTED sign (-)
    kx_exponent = @. βkx_base - βkx_coord * log(clamp(-ψx50, 0.1, 100.0))
    kx = @. ifelse(isnan(aridity_norm), NaN,
                   exp(clamp(kx_exponent, -20.0, 10.0)))
    
    # ΠR calculation
    ΠR_logit = @. βΠR_base + βΠR_slope * aridity_norm
    ΠR = @. ifelse(isnan(aridity_norm), NaN,
                   1 / (1 + exp(-clamp(ΠR_logit, -20.0, 20.0))))
    
    return ψx50, kx, ΠR
end

println("\n" * "="^70)
println("Computing traits for each ensemble member...")
println("="^70)

# Preallocate arrays for ensemble results
n_points = length(aridity)
ψx50_ensemble = Array{Float64}(undef, n_points, n_ensemble)
kx_ensemble = Array{Float64}(undef, n_points, n_ensemble)
ΠR_ensemble = Array{Float64}(undef, n_points, n_ensemble)

# Compute traits for each ensemble member
for i in 1:n_ensemble
    params = param_ensemble[:, i]
    ψx50, kx, ΠR = compute_traits(params, aridity_norm)
    ψx50_ensemble[:, i] = vec(ψx50)
    kx_ensemble[:, i] = vec(kx)
    ΠR_ensemble[:, i] = vec(ΠR)
    
    if i % 5 == 0 || i == n_ensemble
        print("\r  Processed $i / $n_ensemble ensemble members")
    end
end
println()

# Compute statistics (mean and std) for each trait at each point
println("\nComputing ensemble statistics...")

ψx50_mean = vec(mean(ψx50_ensemble, dims=2))
ψx50_std = vec(std(ψx50_ensemble, dims=2))

kx_mean = vec(mean(kx_ensemble, dims=2))
kx_std = vec(std(kx_ensemble, dims=2))

ΠR_mean = vec(mean(ΠR_ensemble, dims=2))
ΠR_std = vec(std(ΠR_ensemble, dims=2))

# Filter to land points for statistics
land_mask = .!isnan.(ψx50_mean)

println("\nTrait statistics (land only, ensemble mean ± std):")
println("\nψx50 (MPa):")
println("  Range: [$(minimum(ψx50_mean[land_mask])), $(maximum(ψx50_mean[land_mask]))]")
println("  Mean uncertainty (std): $(mean(ψx50_std[land_mask]))")

println("\nkx:")
println("  Range: [$(minimum(kx_mean[land_mask])), $(maximum(kx_mean[land_mask]))]")
println("  Mean uncertainty (std): $(mean(kx_std[land_mask]))")
println("  Mean relative uncertainty: $(mean(kx_std[land_mask] ./ kx_mean[land_mask]))")

println("\nΠR:")
println("  Range: [$(minimum(ΠR_mean[land_mask])), $(maximum(ΠR_mean[land_mask]))]")
println("  Mean uncertainty (std): $(mean(ΠR_std[land_mask]))")

# Flatten coordinates
lats_flat = vec(lats)
lons_flat = vec(lons)
aridity_flat = vec(aridity)

# Save results with uncertainties
output_file = joinpath(@__DIR__, "calculated_traits_with_uncertainty.jld2")
JLD2.save(output_file,
    # Means
    "ψx50_mean", ψx50_mean,
    "kx_mean", kx_mean,
    "ΠR_mean", ΠR_mean,
    # Standard deviations
    "ψx50_std", ψx50_std,
    "kx_std", kx_std,
    "ΠR_std", ΠR_std,
    # Coordinates
    "lats", lats_flat,
    "lons", lons_flat,
    "aridity", aridity_flat,
    # Metadata
    "parameter_ensemble", param_ensemble,
    "parameter_names", param_names,
    "parameter_means", param_means,
    "parameter_stds", param_stds
)

println("\n✓ Results saved to: $output_file")
println("✓ Includes mean and std for: ψx50, kx, ΠR")
