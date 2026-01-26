#!/usr/bin/env python3
"""
Visualize hydraulic trait maps (ψx50, kx, ΠR) on the model grid
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Load all three trait datasets
print("Loading trait data...")
psi50_data = np.genfromtxt('traits_psi50_data.csv', delimiter=',', skip_header=1)
kx_data = np.genfromtxt('traits_kx_data.csv', delimiter=',', skip_header=1)
pr_data = np.genfromtxt('traits_pr_data.csv', delimiter=',', skip_header=1)

lons = psi50_data[:, 0]
lats = psi50_data[:, 1]
psi50 = psi50_data[:, 2]
kx = kx_data[:, 2]
pr = pr_data[:, 2]

print(f"Total land points: {len(psi50)}")

# Create figure with 3 rows (one per trait)
fig = plt.figure(figsize=(20, 15))

# --- ψx50 (P50) ---
ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.Robinson())
ax1.coastlines(linewidth=0.5)
ax1.set_global()

# Use log scale for ψx50 to handle extreme values
# Filter out extreme values for visualization
psi50_clipped = np.clip(psi50, -10, -0.1)  # MPa range for visualization

scatter1 = ax1.scatter(lons, lats, c=psi50_clipped, s=2, 
                       cmap='RdYlBu', vmin=-5, vmax=-0.5,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax1.set_title('ψx50 (P50) - Cavitation Resistance [MPa]', fontsize=14, fontweight='bold')
cbar1 = plt.colorbar(scatter1, ax=ax1, orientation='horizontal', pad=0.05, shrink=0.6)
cbar1.set_label('P50 [MPa] (more negative = more resistant)', fontsize=10)

# ψx50 - Americas zoom
ax2 = fig.add_subplot(3, 2, 2, projection=ccrs.PlateCarree())
ax2.coastlines(linewidth=0.8)
ax2.set_extent([-120, -40, -35, 15], crs=ccrs.PlateCarree())

scatter2 = ax2.scatter(lons, lats, c=psi50_clipped, s=8, 
                       cmap='RdYlBu', vmin=-5, vmax=-0.5,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax2.set_title('ψx50 - Americas (Zoomed)', fontsize=13, fontweight='bold')
ax2.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
cbar2 = plt.colorbar(scatter2, ax=ax2, orientation='horizontal', pad=0.08, shrink=0.8)
cbar2.set_label('P50 [MPa]', fontsize=10)

# --- kx (Hydraulic Conductance) ---
ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.Robinson())
ax3.coastlines(linewidth=0.5)
ax3.set_global()

# Use log scale for better visualization across 4 orders of magnitude
kx_log = np.log10(kx + 1e-10)
kx_95 = np.percentile(kx, 95)
kx_log_95 = np.log10(kx_95)

scatter3 = ax3.scatter(lons, lats, c=kx_log, s=2, 
                       cmap='viridis', vmin=-5, vmax=kx_log_95,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax3.set_title('kx - Hydraulic Conductance [log₁₀]', fontsize=14, fontweight='bold')
cbar3 = plt.colorbar(scatter3, ax=ax3, orientation='horizontal', pad=0.05, shrink=0.6)
cbar3.set_label('log₁₀(kx) (higher = more conductive, saturated at 95th percentile)', fontsize=10)

# kx - Americas zoom
ax4 = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree())
ax4.coastlines(linewidth=0.8)
ax4.set_extent([-120, -40, -35, 15], crs=ccrs.PlateCarree())

scatter4 = ax4.scatter(lons, lats, c=kx_log, s=8, 
                       cmap='viridis', vmin=-5, vmax=kx_log_95,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax4.set_title('kx - Americas (Zoomed)', fontsize=13, fontweight='bold')
ax4.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
cbar4 = plt.colorbar(scatter4, ax=ax4, orientation='horizontal', pad=0.08, shrink=0.8)
cbar4.set_label('log₁₀(kx)', fontsize=10)

# --- ΠR (Regulation Strategy) ---
ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.Robinson())
ax5.coastlines(linewidth=0.5)
ax5.set_global()

scatter5 = ax5.scatter(lons, lats, c=pr, s=2, 
                       cmap='RdYlGn_r', vmin=0, vmax=1,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax5.set_title('ΠR - Stomatal Regulation Strategy', fontsize=14, fontweight='bold')
cbar5 = plt.colorbar(scatter5, ax=ax5, orientation='horizontal', pad=0.05, shrink=0.6)
cbar5.set_label('ΠR (0=isohydric, 1=anisohydric)', fontsize=10)

# ΠR - Americas zoom
ax6 = fig.add_subplot(3, 2, 6, projection=ccrs.PlateCarree())
ax6.coastlines(linewidth=0.8)
ax6.set_extent([-120, -40, -35, 15], crs=ccrs.PlateCarree())

scatter6 = ax6.scatter(lons, lats, c=pr, s=8, 
                       cmap='RdYlGn_r', vmin=0, vmax=1,
                       transform=ccrs.PlateCarree(), rasterized=True)
ax6.set_title('ΠR - Americas (Zoomed)', fontsize=13, fontweight='bold')
ax6.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
cbar6 = plt.colorbar(scatter6, ax=ax6, orientation='horizontal', pad=0.08, shrink=0.8)
cbar6.set_label('ΠR', fontsize=10)

plt.tight_layout()
plt.savefig('hydraulic_traits_global.png', dpi=150, bbox_inches='tight')
print(f"\n✓ Saved: hydraulic_traits_global.png")

print("\n" + "="*60)
print("INTERPRETATION:")
print("="*60)
print("ψx50 (P50): Water potential at 50% loss of conductance")
print("  - More negative = more drought-resistant xylem")
print("  - Wet climates: ~-1 MPa, Dry climates: -5 to -10 MPa")
print("\nkx: Hydraulic conductance")
print("  - Higher = more efficient water transport (vulnerable)")
print("  - Lower = safer but less efficient (drought-adapted)")
print("\nΠR: Stomatal regulation strategy")
print("  - 0 (isohydric): Tight stomatal control, maintains leaf water potential")
print("  - 1 (anisohydric): Allows leaf water potential to drop, less regulation")
print("="*60)
