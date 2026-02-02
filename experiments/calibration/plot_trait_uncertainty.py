"""Create global maps showing trait uncertainty from parameter ensemble

This script creates global maps showing the uncertainty in each trait
that comes from parameter uncertainty in the calibration ensemble.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')

print("Loading trait data with uncertainty...")
df = pd.read_csv('traits_with_uncertainty.csv')

# Filter to land points only
df = df[~np.isnan(df['psi50_mean'])].copy()
n_land = len(df)
print(f"Total land points: {n_land}")

print(f"\nTrait statistics:")
print(f"ψx50 mean: {df['psi50_mean'].min():.2f} to {df['psi50_mean'].max():.2f} MPa")
print(f"  Mean uncertainty: {df['psi50_std'].mean():.2f} MPa")
print(f"kx mean: {df['kx_mean'].min():.2f} to {df['kx_mean'].max():.2f}")
print(f"  Mean uncertainty: {df['kx_std'].mean():.2f}")
print(f"ΠR mean: {df['pr_mean'].min():.2f} to {df['pr_mean'].max():.2f}")
print(f"  Mean uncertainty: {df['pr_std'].mean():.2f}")

def create_single_trait_plot(df, trait_col, trait_name, units, cmap='viridis',
                             output_file=None, vmin=None, vmax=None):
    """Create a simple global map of a single trait"""
    
    fig = plt.figure(figsize=(14, 7))
    ax = plt.subplot(111, projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    
    # Plot trait
    if vmin is not None and vmax is not None:
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=df[trait_col].quantile(0.02), 
                        vmax=df[trait_col].quantile(0.98))
    
    sc = ax.scatter(df['lon'], df['lat'], c=df[trait_col],
                   s=0.5, alpha=0.6, cmap=cmap, norm=norm,
                   transform=ccrs.PlateCarree())
    
    ax.set_title(f'{trait_name}', fontsize=16, weight='bold')
    
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, shrink=0.7)
    cbar.set_label(f'{trait_name} {units}', fontsize=12)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    plt.close()

# Generate individual trait maps showing UNCERTAINTY (std dev from ensemble)
print("\nGenerating uncertainty maps...")

# P50 uncertainty
create_single_trait_plot(
    df, 'psi50_std',
    trait_name='P50 (ψ50) Uncertainty',
    units='(MPa)',
    cmap='Reds',
    output_file='trait_uncertainty_psi50.png'
)

# kx uncertainty
create_single_trait_plot(
    df, 'kx_std',
    trait_name='Hydraulic Conductivity (kx) Uncertainty',
    units='',
    cmap='Reds',
    output_file='trait_uncertainty_kx.png'
)

# ΠR uncertainty
create_single_trait_plot(
    df, 'pr_std',
    trait_name='Stomatal Regulation (ΠR) Uncertainty',
    units='',
    cmap='Reds',
    output_file='trait_uncertainty_pr.png'
)

print("\n✓ All uncertainty maps complete!")
print("\nGenerated files:")
print("  - trait_uncertainty_psi50.png")
print("  - trait_uncertainty_kx.png")
print("  - trait_uncertainty_pr.png")
print("\nThese maps show the UNCERTAINTY (std dev) in each trait")
print("due to parameter uncertainty from the 13-member ensemble.")
