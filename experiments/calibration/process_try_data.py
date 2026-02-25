#!/usr/bin/env python3
"""
Process TRY database trait data and compare with model estimates

TRY traits to download:
- P50 (Ψ50): Water potential at 50% loss of conductivity (MPa)
  TRY trait IDs: 163, 3113, 3114
- Hydraulic conductance: Stem or leaf hydraulic conductance
  TRY trait IDs: 233, 234, 235, 236
- Turgor loss point (πtlp): Related to stomatal regulation
  TRY trait IDs: 162, 3115

Instructions:
1. Go to https://www.try-db.org/
2. Request data for the following trait IDs:
   - 163: Leaf water potential at 50% loss of conductivity
   - 3113: Stem P50
   - 3114: Root P50
   - 233: Leaf hydraulic conductance
   - 234: Stem hydraulic conductance
   - 162: Turgor loss point
3. Download the data as CSV
4. Place the file in this directory as 'try_data.txt' or 'try_data.csv'
5. Run this script to process and compare

TRY vs Model P50 Comparison

TRY Database (3,904 observations):
Mean: -3.08 MPa
Median: -2.49 MPa
Range: -22.00 to 0.00 MPa
Std Dev: 2.31 MPa


uSPAC Model (17,912 grid points):
Mean: -2.47 MPa
Median: -2.42 MPa
Range: -3.94 to -0.67 MPa
Std Dev: 1.08 MPa
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

def load_try_data(filepath):
    """Load and clean TRY database export"""
    print(f"Loading TRY data from {filepath}...")
    
    # TRY data is typically tab-delimited
    try:
        df = pd.read_csv(filepath, sep='\t', encoding='ISO-8859-1', low_memory=False)
    except:
        df = pd.read_csv(filepath, encoding='ISO-8859-1', low_memory=False)
    
    print(f"Total records: {len(df)}")
    print(f"Columns: {df.columns.tolist()}")
    
    return df

def extract_p50_data(df):
    """Extract P50 measurements (water potential at 50% loss of conductivity)"""
    print("\n--- Extracting P50 data ---")
    
    # Find trait name column
    trait_col = 'TraitName' if 'TraitName' in df.columns else None
    
    if trait_col is None:
        print("Could not find TraitName column")
        return pd.DataFrame()
    
    # Filter for P50 measurements using actual trait names from TRY
    # These are the specific trait names in the data
    p50_trait_names = [
        'Xylem hydraulic vulnerability curve (P20, P50, P80)',
        'Xylem hydraulic vulnerability, xylem cavitation vulnerability, embolism vulnerability, (P20, P50, P80)'
    ]
    
    p50_df = df[df[trait_col].isin(p50_trait_names)].copy()
    
    print(f"Found {len(p50_df)} P50-related measurements")
    
    if len(p50_df) == 0:
        print("No P50 data found")
        print("Available trait names:")
        print(df[trait_col].value_counts().head(20))
        return pd.DataFrame()
    
    # TRY data structure: DataName contains P20, P50, or P80
    # We only want P50 values
    if 'DataName' in p50_df.columns:
        # Filter for P50 specifically
        p50_df = p50_df[p50_df['DataName'].str.contains('P50', case=False, na=False)].copy()
        print(f"After filtering for P50 specifically: {len(p50_df)} measurements")
    
    # Extract relevant columns
    required_cols = ['StdValue']  # Standard value column
    optional_cols = ['Latitude', 'Longitude', 'AccSpeciesName', 'UnitName']
    
    # Check which columns exist
    available = [col for col in required_cols + optional_cols if col in p50_df.columns]
    
    if 'StdValue' not in available:
        print("ERROR: StdValue column not found")
        return pd.DataFrame()
    
    p50_clean = p50_df[available].copy()
    
    # Rename columns to standard names
    rename_map = {
        'StdValue': 'value',
        'Latitude': 'lat',
        'Longitude': 'lon',
        'AccSpeciesName': 'species',
        'UnitName': 'unit'
    }
    
    p50_clean = p50_clean.rename(columns={k: v for k, v in rename_map.items() if k in p50_clean.columns})
    
    # Convert to numeric
    p50_clean['value'] = pd.to_numeric(p50_clean['value'], errors='coerce')
    
    # Filter: require value
    p50_clean = p50_clean.dropna(subset=['value'])
    
    # If we have lat/lon, convert them too
    if 'lat' in p50_clean.columns and 'lon' in p50_clean.columns:
        p50_clean['lat'] = pd.to_numeric(p50_clean['lat'], errors='coerce')
        p50_clean['lon'] = pd.to_numeric(p50_clean['lon'], errors='coerce')
        p50_clean = p50_clean.dropna(subset=['lat', 'lon'])
        
        print(f"After filtering (with coordinates): {len(p50_clean)} P50 measurements")
        print(f"  Value range: {p50_clean['value'].min():.2f} to {p50_clean['value'].max():.2f}")
        
        if 'unit' in p50_clean.columns:
            print(f"  Units: {p50_clean['unit'].unique()}")
    else:
        print(f"After filtering: {len(p50_clean)} P50 measurements")
        print(f"  WARNING: No lat/lon coordinates found")
        print(f"  Value range: {p50_clean['value'].min():.2f} to {p50_clean['value'].max():.2f}")
    
    return p50_clean

def load_model_traits():
    """Load model-estimated traits"""
    print("\nLoading model trait estimates...")
    
    psi50_df = pd.read_csv('traits_psi50_data.csv')
    kx_df = pd.read_csv('traits_kx_data.csv')
    
    print(f"Model grid points: {len(psi50_df)}")
    print(f"Model P50 range: {psi50_df['psi50'].min():.2f} to {psi50_df['psi50'].max():.2f} MPa")
    
    return psi50_df, kx_df

def bin_observations_to_grid(obs_df, grid_resolution=5):
    """Bin point observations to a regular lat-lon grid"""
    
    # Create grid bins
    lat_bins = np.arange(-90, 90 + grid_resolution, grid_resolution)
    lon_bins = np.arange(-180, 180 + grid_resolution, grid_resolution)
    
    # Assign observations to grid cells
    obs_df['lat_bin'] = pd.cut(obs_df['lat'], bins=lat_bins, labels=lat_bins[:-1] + grid_resolution/2)
    obs_df['lon_bin'] = pd.cut(obs_df['lon'], bins=lon_bins, labels=lon_bins[:-1] + grid_resolution/2)
    
    # Average within grid cells
    gridded = obs_df.groupby(['lat_bin', 'lon_bin']).agg({
        'p50': ['mean', 'std', 'count']
    }).reset_index()
    
    gridded.columns = ['lat', 'lon', 'p50_mean', 'p50_std', 'n_obs']
    gridded['lat'] = pd.to_numeric(gridded['lat'])
    gridded['lon'] = pd.to_numeric(gridded['lon'])
    
    return gridded

def compare_traits(try_data, model_data, trait_name='P50'):
    """Create comparison visualization"""
    
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.2)
    
    # Panel 1: TRY observations
    ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.Robinson())
    ax1.coastlines(linewidth=0.5)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    ax1.set_global()
    
    sc1 = ax1.scatter(try_data['lon'], try_data['lat'], 
                      c=try_data['p50_mean'], s=100, 
                      cmap='RdYlBu', vmin=-5, vmax=-0.5,
                      transform=ccrs.PlateCarree(), 
                      edgecolors='black', linewidth=0.5, alpha=0.8)
    ax1.set_title(f'TRY Database {trait_name} (n={len(try_data)})', fontsize=14, weight='bold')
    plt.colorbar(sc1, ax=ax1, orientation='horizontal', pad=0.05, shrink=0.6, label=f'{trait_name} [MPa]')
    
    # Panel 2: Model estimates
    ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.Robinson())
    ax2.coastlines(linewidth=0.5)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.3)
    ax2.set_global()
    
    sc2 = ax2.scatter(model_data['lon'], model_data['lat'], 
                      c=model_data['psi50'], s=1, 
                      cmap='RdYlBu', vmin=-5, vmax=-0.5,
                      transform=ccrs.PlateCarree(), alpha=0.6)
    ax2.set_title(f'Model Estimated {trait_name}', fontsize=14, weight='bold')
    plt.colorbar(sc2, ax=ax2, orientation='horizontal', pad=0.05, shrink=0.6, label=f'{trait_name} [MPa]')
    
    # Panel 3: Histogram comparison
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.hist(try_data['p50_mean'], bins=30, alpha=0.6, label='TRY observations', 
             color='blue', edgecolor='black')
    ax3.hist(model_data['psi50'], bins=30, alpha=0.6, label='Model estimates', 
             color='red', edgecolor='black')
    ax3.axvline(try_data['p50_mean'].median(), color='blue', linestyle='--', 
                linewidth=2, label=f'TRY median: {try_data["p50_mean"].median():.2f} MPa')
    ax3.axvline(model_data['psi50'].median(), color='red', linestyle='--', 
                linewidth=2, label=f'Model median: {model_data["psi50"].median():.2f} MPa')
    ax3.set_xlabel(f'{trait_name} [MPa]', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Distribution Comparison', fontsize=14, weight='bold')
    ax3.legend()
    ax3.grid(alpha=0.3)
    
    # Panel 4: Statistics table
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    stats_text = f"""
    COMPARISON STATISTICS
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    
    TRY Database:
      N observations: {len(try_data)}
      Mean: {try_data['p50_mean'].mean():.2f} MPa
      Median: {try_data['p50_mean'].median():.2f} MPa
      Std dev: {try_data['p50_mean'].std():.2f} MPa
      Range: [{try_data['p50_mean'].min():.2f}, {try_data['p50_mean'].max():.2f}] MPa
    
    Model Estimates:
      N grid points: {len(model_data)}
      Mean: {model_data['psi50'].mean():.2f} MPa
      Median: {model_data['psi50'].median():.2f} MPa
      Std dev: {model_data['psi50'].std():.2f} MPa
      Range: [{model_data['psi50'].min():.2f}, {model_data['psi50'].max():.2f}] MPa
    
    Difference (Model - TRY):
      Mean difference: {model_data['psi50'].mean() - try_data['p50_mean'].mean():.2f} MPa
      Median difference: {model_data['psi50'].median() - try_data['p50_mean'].median():.2f} MPa
    """
    
    ax4.text(0.1, 0.5, stats_text, fontsize=11, family='monospace',
             verticalalignment='center', bbox=dict(boxstyle='round', 
             facecolor='wheat', alpha=0.3))
    
    plt.suptitle(f'{trait_name} Comparison: TRY Database vs. Model Estimates', 
                 fontsize=16, weight='bold', y=0.98)
    
    plt.savefig(f'try_vs_model_{trait_name.lower()}.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: try_vs_model_{trait_name.lower()}.png")
    
    plt.close()

def main():
    # Check if TRY data exists
    import os
    try_files = ['try_data.txt', 'try_data.csv', 'TRY_data.txt', 'TRY_data.csv']
    try_file = None
    
    for f in try_files:
        if os.path.exists(f):
            try_file = f
            break
    
    if try_file is None:
        print("="*70)
        print("TRY DATABASE DATA NOT FOUND")
        print("="*70)
        print("\nTo compare with TRY data:")
        print("1. Go to https://www.try-db.org/TryWeb/Home.php")
        print("2. Request a data extract with the following traits:")
        print("   - Trait 163: Leaf water potential at 50% loss of conductivity")
        print("   - Trait 3113: Stem P50")
        print("   - Trait 3114: Root P50")
        print("   - Trait 233: Leaf hydraulic conductance")
        print("   - Trait 234: Stem hydraulic conductance")
        print("   - Trait 162: Turgor loss point")
        print("\n3. Download as TXT or CSV")
        print("4. Save as 'try_data.txt' or 'try_data.csv' in:")
        print(f"   {os.getcwd()}/")
        print("\n5. Re-run this script")
        print("="*70)
        return
    
    # Process TRY data
    try_df = load_try_data(try_file)
    p50_obs = extract_p50_data(try_df)
    
    if p50_obs is not None:
        # Grid the observations
        p50_gridded = bin_observations_to_grid(p50_obs, grid_resolution=5)
        
        # Load model estimates
        model_psi50, model_kx = load_model_traits()
        
        # Create comparison
        compare_traits(p50_gridded, model_psi50, trait_name='P50')
        
        print("\n✓ Comparison complete!")
    else:
        print("\nCould not extract P50 data from TRY file.")
        print("Check the file format and column names.")

if __name__ == '__main__':
    main()
