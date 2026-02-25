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

# Load parameter estimates from calibration iteration_002
parameters = (
    βkx_base        =  1.916,   # Baseline log conductivity
    βkx_coord       =  0.745,   # kx-P50 coordination (safety-efficiency trade-off)
    βψx50_base      =  1.371,   # P50 baseline
    βψx50_slope     = -1.765,   # P50 climate sensitivity
    βΠR_base        =  0.334,   # Stomatal strategy baseline (isohydric-anisohydric)
    βΠR_slope       = -0.892    # Stomatal strategy climate sensitivity
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

# Normalize aridity for trait calculations
# For Global Aridity Index (P/ET0):
# 0 = hyper-arid, 0.03-0.2 = arid, 0.2-0.5 = semi-arid, 0.5-0.65 = dry sub-humid
# 0.65-1.0 = humid, >1.0 = very humid, >2.0 = extremely humid
# Normalization: scale by reference value (2.0 = very humid), clamp to [0,1]
# This gives: low norm = dry, high norm = wet
# NOTE: Ocean points are already NaN in the aridity field
aridity_norm = @. clamp(aridity / 2.0, 0.0, 1.0)
# Result: 
#   aridity=0 (hyper-arid) → norm=0.0 (very dry)
#   aridity=2.0 (very humid) → norm=1.0 (very wet)
#   NaN (ocean) stays NaN and propagates through calculations

land_norm = filter(!isnan, aridity_norm)
println("  Normalized aridity - min: $(minimum(land_norm)), max: $(maximum(land_norm)), mean: $(mean(land_norm))")
println("  Median: $(median(land_norm)), 25%: $(quantile(land_norm, 0.25)), 75%: $(quantile(land_norm, 0.75))")

# --- Calculate functional traits from calibrated parameters ---
# This matches the computation in src/standalone/Vegetation/stomatalconductance.jl (lines 393-415)

#     DRY ←─────────── aridity_norm ─────────→ WET
# 0.0                                       1.0
# (aridity_norm: 0 = dry, 1 = wet)

# ψx50:  ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱  (exponential, always negative, INCREASES with aridity_norm)
#     -10 MPa (dry) → -0.5 MPa (wet)

# kx:    ╱╱╱╱╱╱╱╱╱╱╱╱╱╱  (increases with aridity_norm via coordination with ψx50)
#     low → high

# ΠR:    ──────╮
#              ╰───  (sigmoid S-curve)
#     1 (anisohydric, dry) → 0 (isohydric, wet)

# ψx50: Negative water potential using exponential
# Higher aridity_norm (wet) → less negative P50 (less resistant)
# Lower aridity_norm (dry) → more negative P50 (cavitation-resistant)
# Clamp the exponent to prevent extreme values in hyperarid regions
# Use ifelse to skip ocean points (aridity_norm = NaN)
ψx50_exponent = @. ifelse(isnan(aridity_norm), 0.0, 
                          clamp(parameters.βψx50_base + parameters.βψx50_slope * aridity_norm, -2.0, 4.0))
ψx50 = @. ifelse(isnan(aridity_norm), NaN, -exp(ψx50_exponent))

land_psi50 = filter(!isnan, ψx50)
println("\nψx50 (MPa) statistics (land only):")
println("  min: $(minimum(land_psi50)), max: $(maximum(land_psi50)), mean: $(mean(land_psi50)), median: $(median(land_psi50))")

# kx: Hydraulic conductance with safety-efficiency coordination
# Encode the universal trade-off: safety-efficiency spectrum
# literature shows kx and P50 are coordinated
# High kx → vulnerable to cavitation → requires less negative P50 (wet climates)
# Low kx → safer from cavitation → can have very negative P50 (dry climates)
# kx's direct climate response is weak once you account for P50 (Liu et al. 2019)
# 
# COORDINATION EQUATION: kx = exp(βkx_base - βkx_coord * log(-ψx50))
# Note the NEGATIVE sign: as ψx50 becomes more negative (larger -ψx50),
# log(-ψx50) increases, so kx DECREASES (with positive βkx_coord)
# This creates: less negative P50 → high kx (wet), more negative P50 → low kx (dry)
# Clamp log argument to prevent issues with extreme ψx50 values
# Also clamp the full exponent to prevent overflow/underflow in exp()
kx_exponent = @. parameters.βkx_base - parameters.βkx_coord * log(clamp(-ψx50, 0.1, 100.0))
kx = @. ifelse(isnan(aridity_norm), NaN,
               exp(clamp(kx_exponent, -20.0, 10.0)))

land_kx = filter(!isnan, kx)
println("\nkx (hydraulic conductance) statistics (land only):")
println("  min: $(minimum(land_kx)), max: $(maximum(land_kx)), mean: $(mean(land_kx)), median: $(median(land_kx))")

# ΠR: Logistic sigmoid (0 to 1)
# Higher aridity_norm (wet) → lower ΠR (isohydric strategy)
# Lower aridity_norm (dry) → higher ΠR (anisohydric strategy)
# Clamp the sigmoid argument to prevent overflow in exp()
ΠR_logit = @. parameters.βΠR_base + parameters.βΠR_slope * aridity_norm
ΠR = @. ifelse(isnan(aridity_norm), NaN,
               1 / (1 + exp(-clamp(ΠR_logit, -20.0, 20.0))))

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
