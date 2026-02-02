"""Create distribution plots showing trait uncertainty across latitude bands

This script creates violin/KDE plots showing how trait uncertainty varies
with latitude, using the parameter ensemble to show uncertainty.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("Loading trait data with uncertainty...")
df = pd.read_csv('traits_with_uncertainty.csv')

# Filter to land points only
df = df[~np.isnan(df['psi50_mean'])].copy()
n_land = len(df)
print(f"Total land points: {n_land}")

# Create latitude bins
lat_bins = np.linspace(-90, 90, 19)  # ~10 degree bins
df['lat_bin'] = pd.cut(df['lat'], lat_bins, labels=lat_bins[:-1] + 5)

# Filter to bins with sufficient data
min_points = 50
grouped = df.groupby('lat_bin').size()
valid_bins = grouped[grouped >= min_points].index
print(f"Number of latitude bins with >{min_points} points: {len(valid_bins)}")

# Filter to valid bins
df_filtered = df[df['lat_bin'].isin(valid_bins)].copy()

def create_uncertainty_violins(df, lat_bins, trait_mean_col, trait_std_col,
                               trait_name, units, output_file=None,
                               ylim=None):
    """Create violin plots showing mean ± std for each latitude bin"""
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Prepare data
    lat_centers = sorted([float(x) for x in lat_bins])
    positions = np.arange(len(lat_centers))
    
    # Panel 1: Distribution of MEAN values at each latitude
    data_means = []
    for lat in lat_centers:
        subset = df[df['lat_bin'] == lat]
        data_means.append(subset[trait_mean_col].values)
    
    parts1 = ax1.violinplot(data_means, positions=positions, widths=0.7,
                            showmeans=True, showmedians=False)
    
    for pc in parts1['bodies']:
        pc.set_facecolor('#1f77b4')
        pc.set_alpha(0.7)
    
    ax1.set_xticks(positions)
    ax1.set_xticklabels([f'{int(lat)}°' for lat in lat_centers], rotation=45)
    ax1.set_xlabel('Latitude', fontsize=12)
    ax1.set_ylabel(f'{trait_name} (Mean) {units}', fontsize=12)
    ax1.set_title(f'{trait_name} Distribution by Latitude (Spatial Variability)', 
                 fontsize=14, weight='bold')
    ax1.grid(axis='y', alpha=0.3)
    if ylim:
        ax1.set_ylim(ylim)
    
    # Panel 2: Mean UNCERTAINTY at each latitude
    data_stds = []
    for lat in lat_centers:
        subset = df[df['lat_bin'] == lat]
        data_stds.append(subset[trait_std_col].values)
    
    parts2 = ax2.violinplot(data_stds, positions=positions, widths=0.7,
                            showmeans=True, showmedians=False)
    
    for pc in parts2['bodies']:
        pc.set_facecolor('#d62728')
        pc.set_alpha(0.7)
    
    ax2.set_xticks(positions)
    ax2.set_xticklabels([f'{int(lat)}°' for lat in lat_centers], rotation=45)
    ax2.set_xlabel('Latitude', fontsize=12)
    ax2.set_ylabel(f'{trait_name} Uncertainty (Std) {units}', fontsize=12)
    ax2.set_title(f'{trait_name} Parameter Uncertainty by Latitude', 
                 fontsize=14, weight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    plt.suptitle(f'{trait_name}: Spatial Variability vs Parameter Uncertainty',
                fontsize=16, weight='bold', y=0.995)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    plt.close()

def create_combined_distribution(df, lat_bins, trait_mean_col, trait_std_col,
                                 trait_name, units, output_file=None):
    """Create a single plot showing trait value with uncertainty bands by latitude"""
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Prepare data
    lat_centers = sorted([float(x) for x in lat_bins])
    positions = np.arange(len(lat_centers))
    
    means_by_lat = []
    lower_by_lat = []
    upper_by_lat = []
    
    for lat in lat_centers:
        subset = df[df['lat_bin'] == lat]
        # Mean of means
        mean_val = subset[trait_mean_col].mean()
        # Mean of uncertainties
        mean_std = subset[trait_std_col].mean()
        
        means_by_lat.append(mean_val)
        lower_by_lat.append(mean_val - mean_std)
        upper_by_lat.append(mean_val + mean_std)
    
    # Plot mean with uncertainty band
    ax.plot(positions, means_by_lat, 'o-', linewidth=2, markersize=8, 
           color='#1f77b4', label='Ensemble Mean')
    ax.fill_between(positions, lower_by_lat, upper_by_lat, 
                    alpha=0.3, color='#1f77b4', 
                    label='Mean ± Std (Parameter Uncertainty)')
    
    # Also show spatial variability (std of means at each latitude)
    spatial_std_by_lat = []
    for lat in lat_centers:
        subset = df[df['lat_bin'] == lat]
        spatial_std_by_lat.append(subset[trait_mean_col].std())
    
    spatial_lower = np.array(means_by_lat) - np.array(spatial_std_by_lat)
    spatial_upper = np.array(means_by_lat) + np.array(spatial_std_by_lat)
    
    ax.fill_between(positions, spatial_lower, spatial_upper,
                    alpha=0.2, color='green',
                    label='Spatial Variability (std within lat bin)')
    
    ax.set_xticks(positions)
    ax.set_xticklabels([f'{int(lat)}°' for lat in lat_centers], rotation=45)
    ax.set_xlabel('Latitude', fontsize=12)
    ax.set_ylabel(f'{trait_name} {units}', fontsize=12)
    ax.set_title(f'{trait_name} Latitudinal Pattern with Uncertainty', 
                fontsize=14, weight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
    
    plt.close()

# Generate visualizations
print("\nGenerating distribution visualizations...")

# P50 (ψx50)
create_uncertainty_violins(
    df_filtered, valid_bins, 'psi50_mean', 'psi50_std',
    trait_name='P50 (ψ50)',
    units='(MPa)',
    output_file='trait_dist_psi50_uncertainty.png'
)

create_combined_distribution(
    df_filtered, valid_bins, 'psi50_mean', 'psi50_std',
    trait_name='P50 (ψ50)',
    units='(MPa)',
    output_file='trait_dist_psi50_latitudinal.png'
)

# kx
create_uncertainty_violins(
    df_filtered, valid_bins, 'kx_mean', 'kx_std',
    trait_name='Hydraulic Conductivity (kx)',
    units='',
    output_file='trait_dist_kx_uncertainty.png'
)

create_combined_distribution(
    df_filtered, valid_bins, 'kx_mean', 'kx_std',
    trait_name='Hydraulic Conductivity (kx)',
    units='',
    output_file='trait_dist_kx_latitudinal.png'
)

# ΠR
create_uncertainty_violins(
    df_filtered, valid_bins, 'pr_mean', 'pr_std',
    trait_name='Stomatal Regulation (ΠR)',
    units='',
    output_file='trait_dist_pr_uncertainty.png',
    ylim=(0, 1)
)

create_combined_distribution(
    df_filtered, valid_bins, 'pr_mean', 'pr_std',
    trait_name='Stomatal Regulation (ΠR)',
    units='',
    output_file='trait_dist_pr_latitudinal.png'
)

print("\n✓ All distribution visualizations complete!")
print("\nGenerated files:")
print("  For each trait (P50, kx, ΠR):")
print("    - trait_dist_*_uncertainty.png (violin plots: spatial variability vs parameter uncertainty)")
print("    - trait_dist_*_latitudinal.png (latitudinal pattern with uncertainty bands)")
