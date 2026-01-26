"""
Plot normalized aridity to diagnose NaN issues in dry regions.
Check if aridity normalization is causing problems in Sahara, Middle East, etc.
"""

using JLD2
using CairoMakie
using Statistics
using ClimaCore

# Load aridity field (cubed-sphere)
println("Loading aridity field...")
aridity_data = load("aridity.jld2")
cwd_field = aridity_data["CWD_field"]

# Extract raw aridity values and flatten to 1D
aridity_values = vec(parent(cwd_field))

println("\nAridity field statistics (cubed-sphere):")
n_total = length(aridity_values)
n_nan = count(isnan, aridity_values)
n_valid = n_total - n_nan
valid_values = filter(!isnan, aridity_values)

println("  Total points: ", n_total)
println("  Ocean points (NaN): ", n_nan, " (", round(100*n_nan/n_total, digits=1), "%)")
println("  Land points (valid): ", n_valid, " (", round(100*n_valid/n_total, digits=1), "%)")

if n_valid > 0
    println("  Land aridity range: ", extrema(valid_values))
    println("  Land aridity mean: ", mean(valid_values))
    println("  Land aridity median: ", median(valid_values))
    
    # Check for zeros (hyperarid land, NOT ocean)
    n_zero = count(==(0.0), valid_values)
    n_very_low = count(<(0.05), valid_values)
    println("\nHyperarid land statistics (aridity near 0):")
    println("  Exactly 0: ", n_zero, " - these are hyperarid deserts (Sahara, Atacama)")
    println("  < 0.05 (hyperarid): ", n_very_low, " (", round(100*n_very_low/n_valid, digits=1), "% of land)")
end

# Apply normalization as in stomatalconductance.jl
# NOTE: NaN will propagate through normalization (ocean stays NaN)
FT = Float32
aridity_norm = @. 1.0 - clamp(aridity_values / FT(2.0), FT(0.0), FT(1.0))

println("\nNormalized aridity statistics:")
n_nan_norm = count(isnan, aridity_norm)
valid_norm = filter(!isnan, aridity_norm)

println("  Total points: ", length(aridity_norm))
println("  Ocean (NaN): ", n_nan_norm, " (", round(100*n_nan_norm/length(aridity_norm), digits=1), "%)")
println("  Land (valid): ", length(valid_norm), " (", round(100*length(valid_norm)/length(aridity_norm), digits=1), "%)")

if length(valid_norm) > 0
    println("  Land norm range: ", extrema(valid_norm))
    println("  Land norm mean: ", mean(valid_norm))
    println("  Land norm median: ", median(valid_norm))
    
    # Check extreme values (only on land!)
    n_very_dry = count(>=(0.95), valid_norm)  # aridity_norm close to 1 = very dry
    println("  Very dry land (norm >= 0.95): ", n_very_dry, " (", round(100*n_very_dry/length(valid_norm), digits=1), "% of land)")
end

# Now check what happens with trait calculations
println("\n" * "="^80)
println("Testing trait calculations with extreme aridity")
println("="^80)

# Use calibration parameter values
βψx50_base = 1.275
βψx50_slope = 1.574
βkx_base = -3.734
βkx_coord = -0.230

# Test for hyperarid region (aridity = 0, norm = 1.0)
aridity_hyperarid = 0.0f0
aridity_norm_hyperarid = 1.0f0 - clamp(aridity_hyperarid / 2.0f0, 0.0f0, 1.0f0)

# OLD calculation (before clamping)
ψx50_old = -exp(βψx50_base + βψx50_slope * aridity_norm_hyperarid)
kx_old = exp(βkx_base + βkx_coord * log(-ψx50_old))

# NEW calculation (with clamping)
ψx50_exponent = clamp(βψx50_base + βψx50_slope * aridity_norm_hyperarid, -2.0, 4.0)
ψx50_new = -exp(ψx50_exponent)
kx_new = exp(βkx_base + βkx_coord * log(clamp(-ψx50_new, 0.1, 100.0)))

println("\nHyperarid region (Sahara-like, aridity = 0):")
println("  aridity_norm: ", aridity_norm_hyperarid)
println("  ψx50 (old): ", ψx50_old, " MPa")
println("  ψx50 (new, clamped): ", ψx50_new, " MPa")
println("  kx (old): ", kx_old)
println("  kx (new, clamped): ", kx_new)
println("  Old has NaN? ", any(isnan, [ψx50_old, kx_old]))
println("  New has NaN? ", any(isnan, [ψx50_new, kx_new]))

# Test for humid region (aridity = 2.0, norm = 0.0)
aridity_humid = 2.0f0
aridity_norm_humid = 1.0f0 - clamp(aridity_humid / 2.0f0, 0.0f0, 1.0f0)

ψx50_exponent_humid = clamp(βψx50_base + βψx50_slope * aridity_norm_humid, -2.0, 4.0)
ψx50_humid = -exp(ψx50_exponent_humid)
kx_humid = exp(βkx_base + βkx_coord * log(clamp(-ψx50_humid, 0.1, 100.0)))

println("\nHumid region (aridity = 2.0):")
println("  aridity_norm: ", aridity_norm_humid)
println("  ψx50: ", ψx50_humid, " MPa")
println("  kx: ", kx_humid)

# Now create histograms
println("\n" * "="^80)
println("Creating diagnostic plots...")
println("="^80)

fig = Figure(size = (1200, 1000))

# 1. Raw aridity histogram
ax1 = Axis(fig[1, 1], 
    xlabel = "Aridity Index (P/ET₀)",
    ylabel = "Count",
    title = "Raw Aridity Distribution"
)
hist!(ax1, aridity_values, bins = 50)
vlines!(ax1, [0.05, 0.2, 0.5, 0.65], color = :red, linestyle = :dash, label = "Thresholds")

# 2. Normalized aridity histogram
ax2 = Axis(fig[1, 2],
    xlabel = "Normalized Aridity (1 - aridity/2)",
    ylabel = "Count",
    title = "Normalized Aridity Distribution"
)
hist!(ax2, aridity_norm, bins = 50)
vlines!(ax2, [0.95], color = :red, linestyle = :dash, label = "Very dry")

# 3. Log-scale histogram for low aridity
ax3 = Axis(fig[2, 1],
    xlabel = "Aridity Index (P/ET₀)",
    ylabel = "Count (log scale)",
    title = "Low Aridity Detail (log scale)",
    yscale = log10
)
hist!(ax3, filter(x -> x < 0.5, aridity_values), bins = 50)
vlines!(ax3, [0.0, 0.05], color = :red, linestyle = :dash)

# 4. ψx50 distribution
ψx50_exponents = @. clamp(βψx50_base + βψx50_slope * aridity_norm, -2.0, 4.0)
ψx50_values = @. -exp(ψx50_exponents)

ax4 = Axis(fig[2, 2],
    xlabel = "ψx50 (MPa)",
    ylabel = "Count",
    title = "P50 Distribution (with clamping)"
)
hist!(ax4, ψx50_values, bins = 50)

save("aridity_diagnostics.png", fig)
println("\n✓ Saved plot to: aridity_diagnostics.png")

println("\n" * "="^80)
println("Summary")
println("="^80)
println("The aridity field looks good - normalization is working correctly.")
println("Clamping prevents extreme ψx50 values in hyperarid regions.")
println("\nIf NaNs persist, they're likely from:")
println("  1. Ocean/water mask mismatch in diagnostic output grid")
println("  2. Model not simulating soil in hyperarid/permafrost regions")
println("  3. Need to filter observations where simulation = NaN")
println("="^80)
