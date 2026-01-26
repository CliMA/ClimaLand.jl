"""Calculate ψx50, kx, and ΠR from calibrated parameter estimates

This script uses the parameter estimates from iteration_004/member_001/parameters.toml
and the Global Aridity Index from aridity.jld2 to compute the functional traits
according to the equations in the uSPAC conductance model.
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
# CWD_field is a ClimaCore Field
cwd_field = data["CWD_field"]
println("\nType of aridity field: $(typeof(cwd_field))")

# Get the spectral element space
field_space = ClimaCore.Fields.axes(cwd_field)

# Extract actual lat/lon coordinates from the Field
coords = ClimaCore.Fields.coordinate_field(field_space)
lats_field = coords.lat
lons_field = coords.long

# Access the underlying array directly from the Field
# The field structure has: values -> array (the actual data)
aridity = parent(cwd_field)  # Global Aridity Index (P/ET0)
lats = parent(lats_field)
lons = parent(lons_field)

# Load parameter estimates
parameters = (
    βkx_base = -6.5,
    βΠR_base = -4.4,      # Adjusted for percentile-based normalization
    βψx50_slope = 1.2,    # Reduced slope for realistic P50 range (-1 to -10 MPa)
    βkx_coord = -2.5,     # Increased coordination for more kx variation (2-3 orders of magnitude)
    βψx50_base = 0.2,     # Baseline for P50 equation
    βΠR_slope = 6.0       # Slope for isohydric-anisohydric transition
)

println("\nCalibrated parameters:")
for (k, v) in pairs(parameters)
    println("  $k = $v")
end

println("\nData shape:")
println("  Aridity: $(size(aridity))")
println("  Lats: $(size(lats)), range: $(extrema(lats))")
println("  Lons: $(size(lons)), range: $(extrema(lons))")
println("  Total points: $(length(aridity))")

println("\nAridity statistics:")
land_aridity = filter(!isnan, aridity)
println("  Global Aridity Index (P/ET0) - min: $(minimum(land_aridity)), max: $(maximum(land_aridity)), mean: $(mean(land_aridity))")

# Normalization: For Global Aridity Index, typical values range 0-6.5
# Ocean points are NaN and will be filtered out from statistics
# We normalize to get aridity_norm in a reasonable range for the trait equations

# Normalize aridity for trait calculations
# Raw aridity: high = wet (more P/ET0), low = dry (less P/ET0)
# We INVERT so that: high norm = dry, low norm = wet
# Use 95th percentile as reference to avoid skewing by extreme humid values
using Statistics
ref_aridity = quantile(land_aridity, 0.95)  # 95th percentile = "wet" reference
aridity_norm = @. (ref_aridity - aridity) / ref_aridity
# Result: wet (high aridity) → low norm (~0), dry (low aridity) → high norm (~1)
# NaN (ocean) stays NaN and propagates through calculations

land_norm = filter(!isnan, aridity_norm)
#println("  Reference aridity (95th percentile): $(ref_aridity)")
println("  Inverted and normalized - min: $(minimum(land_norm)), max: $(maximum(land_norm)), mean: $(mean(land_norm))")
println("  Median: $(median(land_norm)), 25%: $(quantile(land_norm, 0.25)), 75%: $(quantile(land_norm, 0.75))")

# --- Calculate functional traits from calibrated parameters ---

# ψx50: Negative water potential using exponential
# Physical meaning: Dry plants have more cavitation-resistant xylem (more negative P50)
# WET (low aridity_norm) ←───── aridity_norm ─────→ DRY (high aridity_norm)
# With POSITIVE slope: high aridity_norm → large positive exponent → very negative P50
# ψx50:  ╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲  (more negative in dry)
#     -1 MPa (wet) → -10^36 MPa (dry)
ψx50 = @. -exp(parameters.βψx50_base + parameters.βψx50_slope * aridity_norm)

land_psi50 = filter(!isnan, ψx50)
println("\nψx50 (MPa) statistics (land only):")
println("  min: $(minimum(land_psi50)), max: $(maximum(land_psi50)), mean: $(mean(land_psi50)), median: $(median(land_psi50))")

# kx: Hydraulic conductance with safety-efficiency coordination
# Physical meaning: Encodes the universal trade-off in the safety-efficiency spectrum
# Literature shows kx and P50 are coordinated (Liu et al. 2019)
# High kx → vulnerable to cavitation → requires less negative P50 (wet climates)
# Low kx → safer from cavitation → can have very negative P50 (dry climates)
# kx:    ╱╱╱╱╱╱╱╱╱╱╱╱╱╱  (increases with aridity through P50 coordination)
#     low → high
kx = @. exp(parameters.βkx_base + parameters.βkx_coord * log(-ψx50))

land_kx = filter(!isnan, kx)
println("\nkx (hydraulic conductance) statistics (land only):")
println("  min: $(minimum(land_kx)), max: $(maximum(land_kx)), mean: $(mean(land_kx)), median: $(median(land_kx))")

# ΠR: Isohydric-anisohydric regulation strategy
# Physical meaning: Dry plants are more anisohydric 
# (tolerate lower leaf water potential, less stomatal regulation)
# WET (low aridity_norm) ←───── aridity_norm ─────→ DRY (high aridity_norm)
# With POSITIVE slope: high aridity_norm → high ΠR (anisohydric)
# ΠR:        ╭──────  (sigmoid, increases with aridity_norm)
#     ───╯
#     0 (isohydric, wet) → 1 (anisohydric, dry)
ΠR = @. 1 / (1 + exp(-(parameters.βΠR_base + parameters.βΠR_slope * aridity_norm)))

land_pr = filter(!isnan, ΠR)
println("\nΠR (regulation strategy) statistics (land only):")
println("  min: $(minimum(land_pr)), max: $(maximum(land_pr)), mean: $(mean(land_pr)), median: $(median(land_pr))")

# --- Create proper land/sea mask using topography data ---
# Don't use CWD threshold - use actual topographic land/sea mask
println("\n--- Creating land/sea mask from topography ---")

# Get the spectral element grid configuration
# The CWD field was created with this grid: nelements = (101, 15)
# This is from generate_cwd_field.jl line 11
ne = (101, 15)

println("  Grid: ne=$ne ($(ne[1]) × $(ne[2]) elements)")

# Create domain for mask calculation (CPU only, don't waste GPU)
domain = ClimaLand.Domains.SphericalShell(;
    radius = 0.1,
    depth = 0.1,
    nelements = ne,
    context = ClimaComms.SingletonCommsContext{ClimaComms.CPUSingleThreaded}(
        ClimaComms.CPUSingleThreaded(),
    ),
)

# Get proper land/sea mask based on topography (threshold = 0.99 means >99% land)
landsea_mask_field = ClimaLand.Domains.landsea_mask(domain; threshold = 0.99)
land_mask = vec(parent(landsea_mask_field)) .== 1.0

n_land = count(land_mask)
n_ocean = count(.!land_mask)

println("  Land points (topography >99% land): $n_land ($(round(100*n_land/length(land_mask), digits=1))%)")
println("  Ocean points (topography >1% water): $n_ocean ($(round(100*n_ocean/length(land_mask), digits=1))%)")

# Compare with CWD-based mask
cwd_based_mask = vec(aridity) .> 0.1
n_cwd_land = count(cwd_based_mask)
println("\n  For comparison:")
println("  CWD>0.1 would give: $n_cwd_land land points ($(round(100*n_cwd_land/length(cwd_based_mask), digits=1))%)")
println("  Difference: $(abs(n_land - n_cwd_land)) points")

# Set ocean points to NaN
ocean_mask = .!land_mask
ψx50[ocean_mask] .= NaN
kx[ocean_mask] .= NaN
ΠR[ocean_mask] .= NaN

println("\nAfter masking - Land statistics only:")
println("  ψx50 range: $(extrema(ψx50[land_mask]))")
println("  kx range: $(extrema(kx[land_mask]))")
println("  ΠR range: $(extrema(ΠR[land_mask]))")

# Flatten all arrays to 1D vectors for easier analysis
aridity_flat = vec(aridity)
lats_flat = vec(lats)
lons_flat = vec(lons)
ψx50_flat = vec(ψx50)
kx_flat = vec(kx)
ΠR_flat = vec(ΠR)

println("\nFlattened to $(length(aridity_flat)) points ($n_land land, $n_ocean ocean)")

# Save results with coordinates
output_file = joinpath(@__DIR__, "calculated_traits.jld2")
JLD2.save(output_file, 
    "ψx50", ψx50_flat,
    "kx", kx_flat,
    "ΠR", ΠR_flat,
    "aridity", aridity_flat,
    "lats", lats_flat,
    "lons", lons_flat,
    "parameters", parameters
)

println("\n✓ Results saved to: $output_file")
println("✓ Done! Use Python or other tools for visualization.")
