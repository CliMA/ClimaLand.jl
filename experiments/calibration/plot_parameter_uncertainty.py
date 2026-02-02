"""Visualize trait distributions with parameter uncertainty

This script loads the calculated traits with uncertainty from the Julia script
and creates maps showing both the mean trait values and their uncertainties.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
print(f"ψx50: {df['psi50_mean'].min():.2f} to {df['psi50_mean'].max():.2f} MPa")
print(f"  Mean uncertainty: {df['psi50_std'].mean():.2f} MPa")
print(f"kx: {df['kx_mean'].min():.2f} to {df['kx_mean'].max():.2f}")
print(f"  Mean uncertainty: {df['kx_std'].mean():.2f}")
print(f"ΠR: {df['pr_mean'].min():.2f} to {df['pr_mean'].max():.2f}")
print(f"  Mean uncertainty: {df['pr_std'].mean():.2f}")

def create_trait_map_with_uncertainty(df, trait_mean_col, trait_std_col, 
                                       trait_name, units, cmap_mean='RdYlBu_r',
                                       output_file=None, vmin=None, vmax=None):
    """Create a 2-panel map: mean and uncertainty"""
    
    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace=0.15)
    
    # Panel 1: Mean trait value
    ax1 = plt.subplot(gs[0], projection=ccrs.Robinson())
    ax1.set_global()
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    
    # Plot mean
    if vmin is not None and vmax is not None:
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=df[trait_mean_col].quantile(0.02), 
                        vmax=df[trait_mean_col].quantile(0.98))
    
    sc1 = ax1.scatter(df['lon'], df['lat'], c=df[trait_mean_col],
                     s=0.5, alpha=0.6, cmap=cmap_mean, norm=norm,
                     transform=ccrs.PlateCarree())
    
    ax1.set_title(f'{trait_name} (Ensemble Mean)', fontsize=14, weight='bold')
    
    # Colorbar for mean
    cbar1 = plt.colorbar(sc1, ax=ax1, orientation='horizontal', pad=0.05, shrink=0.8)
    cbar1.set_label(f'{trait_name} {units}', fontsize=12)
    
    # Panel 2: Uncertainty (std)
    ax2 = plt.subplot(gs[1], projection=ccrs.Robinson())
    ax2.set_global()
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    
    # Plot uncertainty
    norm_std = Normalize(vmin=0, vmax=df[trait_std_col].quantile(0.98))
    
    sc2 = ax2.scatter(df['lon'], df['lat'], c=df[trait_std_col],
                     s=0.5, alpha=0.6, cmap='Reds', norm=norm_std,
                     transform=ccrs.PlateCarree())
    
    ax2.set_title(f'{trait_name} Uncertainty (Std Dev)', fontsize=14, weight='bold')
    
    # Colorbar for uncertainty
    cbar2 = plt.colorbar(sc2, ax=ax2, orientation='horizontal', pad=0.05, shrink=0.8)
    cbar2.set_label(f'Uncertainty {units}', fontsize=12)
    
    plt.suptitle(f'{trait_name} from Parameter Ensemble (13 members)', 
                fontsize=16, weight='bold', y=0.98)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    plt.close()

def create_coefficient_of_variation_map(df, trait_mean_col, trait_std_col,
                                         trait_name, output_file=None):
    """Create map of coefficient of variation (CV = std/mean)"""
    
    # Calculate CV (relative uncertainty)
    cv = np.abs(df[trait_std_col] / df[trait_mean_col])
    
    fig = plt.figure(figsize=(14, 7))
    ax = plt.subplot(111, projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    
    # Plot CV
    norm = Normalize(vmin=0, vmax=np.percentile(cv, 95))
    
    sc = ax.scatter(df['lon'], df['lat'], c=cv,
                   s=0.5, alpha=0.6, cmap='YlOrRd', norm=norm,
                   transform=ccrs.PlateCarree())
    
    ax.set_title(f'{trait_name} Coefficient of Variation (σ/μ)', 
                fontsize=16, weight='bold')
    
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, shrink=0.7)
    cbar.set_label('Relative Uncertainty (CV)', fontsize=12)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    plt.close()

# Generate visualizations
print("\nGenerating visualizations...")

# P50 (ψx50)
create_trait_map_with_uncertainty(
    df, 'psi50_mean', 'psi50_std',
    trait_name='P50 (ψ50)',
    units='(MPa)',
    cmap_mean='RdYlBu',  # Blue = less negative (wet), Red = more negative (dry)
    output_file='trait_uncertainty_psi50_map.png'
)

create_coefficient_of_variation_map(
    df, 'psi50_mean', 'psi50_std',
    trait_name='P50 (ψ50)',
    output_file='trait_uncertainty_psi50_cv.png'
)

# kx
create_trait_map_with_uncertainty(
    df, 'kx_mean', 'kx_std',
    trait_name='Hydraulic Conductivity (kx)',
    units='',
    cmap_mean='viridis',
    output_file='trait_uncertainty_kx_map.png'
)

create_coefficient_of_variation_map(
    df, 'kx_mean', 'kx_std',
    trait_name='Hydraulic Conductivity (kx)',
    output_file='trait_uncertainty_kx_cv.png'
)

# ΠR
create_trait_map_with_uncertainty(
    df, 'pr_mean', 'pr_std',
    trait_name='Stomatal Regulation (ΠR)',
    units='',
    cmap_mean='RdYlGn_r',  # Red = high regulation, Green = low regulation
    output_file='trait_uncertainty_pr_map.png',
    vmin=0, vmax=1
)

create_coefficient_of_variation_map(
    df, 'pr_mean', 'pr_std',
    trait_name='Stomatal Regulation (ΠR)',
    output_file='trait_uncertainty_pr_cv.png'
)

print("\n✓ All visualizations complete!")
print("\nGenerated files:")
print("  - trait_uncertainty_psi50_map.png (mean + uncertainty)")
print("  - trait_uncertainty_psi50_cv.png (coefficient of variation)")
print("  - trait_uncertainty_kx_map.png (mean + uncertainty)")
print("  - trait_uncertainty_kx_cv.png (coefficient of variation)")
print("  - trait_uncertainty_pr_map.png (mean + uncertainty)")
print("  - trait_uncertainty_pr_cv.png (coefficient of variation)")
