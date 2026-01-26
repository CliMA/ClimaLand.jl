"""
Plot normalized aridity on a world map to visualize spatial distribution.
"""

using JLD2
using CairoMakie
using GeoMakie
using Statistics
using ClimaCore
using ClimaCore.Fields

# Load aridity field (cubed-sphere)
println("Loading aridity field...")
aridity_data = load("aridity.jld2")
cwd_field = aridity_data["CWD_field"]

# The field includes space information - extract it properly
space = axes(cwd_field)

println("Extracting coordinates...")
# Get coordinate field from the space
coord_field = Fields.coordinate_field(space)
lons_raw = Array(parent(coord_field.long))
lats_raw = Array(parent(coord_field.lat))

# Flatten
lons_flat = vec(lons_raw)
lats_flat = vec(lats_raw)
aridity_flat = vec(parent(cwd_field))

println("  Points: ", length(aridity_flat))
println("  Lon range: ", extrema(lons_flat))
println("  Lat range: ", extrema(lats_flat))

# Apply normalization
FT = Float32
aridity_norm = @. 1.0 - clamp(aridity_flat / FT(2.0), FT(0.0), FT(1.0))

println("\nCreating map plots...")

# Create figure with multiple panels
fig = Figure(size = (1400, 1000))

# 1. Raw aridity map
ax1 = GeoAxis(fig[1, 1]; 
    dest = "+proj=robin",
    title = "Raw Aridity Index (P/ET₀)"
)

sc1 = scatter!(ax1, lons_flat, lats_flat; 
    color = aridity_flat,
    markersize = 3,
    colormap = :viridis,
    colorrange = (0, 2.0)
)
Colorbar(fig[1, 2], sc1, label = "Aridity")

# 2. Normalized aridity map
ax2 = GeoAxis(fig[2, 1];
    dest = "+proj=robin", 
    title = "Normalized Aridity (1 - aridity/2)"
)

sc2 = scatter!(ax2, lons_flat, lats_flat;
    color = aridity_norm,
    markersize = 3,
    colormap = :RdYlBu,  # Red = dry, Blue = wet
    colorrange = (0, 1)
)
Colorbar(fig[2, 2], sc2, label = "Normalized (1=dry, 0=wet)")

# 3. Hyperarid regions (aridity = 0)
ax3 = GeoAxis(fig[3, 1];
    dest = "+proj=robin",
    title = "Hyperarid Regions (aridity = 0)"
)

hyperarid_mask = aridity_flat .== 0.0
scatter!(ax3, lons_flat[hyperarid_mask], lats_flat[hyperarid_mask];
    color = :red,
    markersize = 3
)
scatter!(ax3, lons_flat[.!hyperarid_mask], lats_flat[.!hyperarid_mask];
    color = :lightgray,
    markersize = 1,
    alpha = 0.3
)

# 4. Very dry regions (norm >= 0.95)
ax4 = GeoAxis(fig[3, 2];
    dest = "+proj=robin",
    title = "Very Dry Regions (norm ≥ 0.95)"
)

very_dry_mask = aridity_norm .>= 0.95
scatter!(ax4, lons_flat[very_dry_mask], lats_flat[very_dry_mask];
    color = :orange,
    markersize = 3
)
scatter!(ax4, lons_flat[.!very_dry_mask], lats_flat[.!very_dry_mask];
    color = :lightgray,
    markersize = 1,
    alpha = 0.3
)

save("aridity_maps.png", fig)
println("\n✓ Saved map to: aridity_maps.png")

# Print statistics by region
println("\n" * "="^80)
println("Regional statistics")
println("="^80)

function region_stats(lons, lats, values, name)
    println("\n$name:")
    println("  Points: ", length(values))
    println("  Mean aridity_norm: ", round(mean(values), digits=3))
    println("  Median: ", round(median(values), digits=3))
    println("  Range: ", round.(extrema(values), digits=3))
end

# Define regions
africa = (lats_flat .>= -35) .& (lats_flat .<= 35) .& (lons_flat .>= -20) .& (lons_flat .<= 55)
middle_east = (lats_flat .>= 15) .& (lats_flat .<= 45) .& (lons_flat .>= 30) .& (lons_flat .<= 70)
arctic = lats_flat .>= 60
tropics = (lats_flat .>= -23.5) .& (lats_flat .<= 23.5)

region_stats(lons_flat[africa], lats_flat[africa], aridity_norm[africa], "Africa")
region_stats(lons_flat[middle_east], lats_flat[middle_east], aridity_norm[middle_east], "Middle East")
region_stats(lons_flat[arctic], lats_flat[arctic], aridity_norm[arctic], "Arctic (>60°N)")
region_stats(lons_flat[tropics], lats_flat[tropics], aridity_norm[tropics], "Tropics (±23.5°)")

println("\n" * "="^80)
println("The map should show:")
println("  - Red/orange in Sahara, Arabian Desert, Central Asia (hyperarid)")
println("  - Blue in Amazon, Congo, Southeast Asia (humid)")
println("  - Yellow/green in temperate regions")
println("  - Check if the distribution matches expected climate zones")
println("="^80)
