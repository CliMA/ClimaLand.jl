#!/usr/bin/env python3
"""
Compare TRY P50 distribution with model estimates (no spatial matching)

Since the TRY data lacks coordinates, we'll compare:
1. Overall distributions
2. Statistical moments
3. Range comparison

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
from matplotlib.gridspec import GridSpec

print("="*70)
print("TRY P50 Distribution vs Model Estimates")
print("="*70)

# Load TRY data
print("\nLoading TRY data...")
try_df =pd.read_csv('TRY_data.txt', sep='\t', low_memory=False, encoding='ISO-8859-1')

# Extract P50
p50_trait_names = [
    'Xylem hydraulic vulnerability curve (P20, P50, P80)',
    'Xylem hydraulic vulnerability, xylem cavitation vulnerability, embolism vulnerability, (P20, P50, P80)'
]

p50_df = try_df[try_df['TraitName'].isin(p50_trait_names)]
p50_df = p50_df[p50_df['DataName'].str.contains('P50', case=False, na=False)]
p50_df['value'] = pd.to_numeric(p50_df['StdValue'], errors='coerce')
p50_df = p50_df.dropna(subset=['value'])

# Filter outliers (some unrealistic values in TRY)
p50_df = p50_df[(p50_df['value'] >= -25) & (p50_df['value'] <= 0)]

print(f"  TRY P50 observations: {len(p50_df)}")
print(f"  Range: {p50_df['value'].min():.2f} to {p50_df['value'].max():.2f} MPa")
print(f"  Mean: {p50_df['value'].mean():.2f} MPa")
print(f"  Median: {p50_df['value'].median():.2f} MPa")
print(f"  Std: {p50_df['value'].std():.2f} MPa")

# Load model data
print("\nLoading model estimates...")
model_df = pd.read_csv('traits_psi50_data.csv')

print(f"  Model grid points: {len(model_df)}")
print(f"  Range: {model_df['psi50'].min():.2f} to {model_df['psi50'].max():.2f} MPa")
print(f"  Mean: {model_df['psi50'].mean():.2f} MPa") 
print(f"  Median: {model_df['psi50'].median():.2f} MPa")
print(f"  Std: {model_df['psi50'].std():.2f} MPa")

# Create comparison figure
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3, bottom=0.35, top=0.92)

# Panel 1: Overlapping histograms
ax1 = fig.add_subplot(gs[0, :])
bins = np.linspace(-10, 0, 40)
ax1.hist(p50_df['value'], bins=bins, alpha=0.6, label=f'TRY observations (n={len(p50_df)})', 
         color='steelblue', edgecolor='navy', linewidth=0.5, density=True)
ax1.hist(model_df['psi50'], bins=bins, alpha=0.6, label=f'Model estimates (n={len(model_df)})', 
         color='coral', edgecolor='darkred', linewidth=0.5, density=True)

# Add median lines
ax1.axvline(p50_df['value'].median(), color='navy', linestyle='--', linewidth=2,
            label=f'TRY median: {p50_df["value"].median():.2f} MPa')
ax1.axvline(model_df['psi50'].median(), color='darkred', linestyle='--', linewidth=2,
            label=f'Model median: {model_df["psi50"].median():.2f} MPa')

ax1.set_xlabel('P50 (MPa)', fontsize=13)
ax1.set_ylabel('Probability Density', fontsize=13)
ax1.set_title('P50 Distribution Comparison', fontsize=15, fontweight='bold')
ax1.legend(fontsize=11, framealpha=0.9)
ax1.grid(alpha=0.3, axis='y')

# Panel 2: Cumulative distributions
ax2 = fig.add_subplot(gs[1, 0])

try_sorted = np.sort(p50_df['value'].values)
try_cumsum = np.arange(1, len(try_sorted)+1) / len(try_sorted)

model_sorted = np.sort(model_df['psi50'].values)
model_cumsum = np.arange(1, len(model_sorted)+1) / len(model_sorted)

ax2.plot(try_sorted, try_cumsum, linewidth=2.5, color='steelblue', label='TRY observations')
ax2.plot(model_sorted, model_cumsum, linewidth=2.5, color='coral', label='Model estimates')

ax2.set_xlabel('P50 (MPa)', fontsize=12)
ax2.set_ylabel('Cumulative Probability', fontsize=12)
ax2.set_title('Cumulative Distribution Function', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

# Panel 3: Box plots
ax3 = fig.add_subplot(gs[1, 1])

box_data = [p50_df['value'].values, model_df['psi50'].values]
bp = ax3.boxplot(box_data, labels=['TRY\nObservations', 'Model\nEstimates'],
                 patch_artist=True, widths=0.6,
                 medianprops=dict(color='black', linewidth=2),
                 boxprops=dict(facecolor='lightblue', edgecolor='navy', linewidth=1.5),
                 whiskerprops=dict(color='navy', linewidth=1.5),
                 capprops=dict(color='navy', linewidth=1.5))

# Color the boxes differently
bp['boxes'][0].set_facecolor('steelblue')
bp['boxes'][1].set_facecolor('coral')

ax3.set_ylabel('P50 (MPa)', fontsize=12)
ax3.set_title('Distribution Statistics', fontsize=14, fontweight='bold')
ax3.grid(alpha=0.3, axis='y')

# Add text with statistics
stats_text = f"""
COMPARISON SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TRY Database (n={len(p50_df)}):
  Mean:     {p50_df['value'].mean():6.2f} MPa
  Median:   {p50_df['value'].median():6.2f} MPa
  Std Dev:  {p50_df['value'].std():6.2f} MPa
  Range:    [{p50_df['value'].min():.2f}, {p50_df['value'].max():.2f}] MPa
  
  Percentiles:
    10th:   {p50_df['value'].quantile(0.10):6.2f} MPa
    25th:   {p50_df['value'].quantile(0.25):6.2f} MPa
    75th:   {p50_df['value'].quantile(0.75):6.2f} MPa
    90th:   {p50_df['value'].quantile(0.90):6.2f} MPa

Model Estimates (n={len(model_df)}):
  Mean:     {model_df['psi50'].mean():6.2f} MPa
  Median:   {model_df['psi50'].median():6.2f} MPa
  Std Dev:  {model_df['psi50'].std():6.2f} MPa
  Range:    [{model_df['psi50'].min():.2f}, {model_df['psi50'].max():.2f}] MPa
  
  Percentiles:
    10th:   {model_df['psi50'].quantile(0.10):6.2f} MPa
    25th:   {model_df['psi50'].quantile(0.25):6.2f} MPa
    75th:   {model_df['psi50'].quantile(0.75):6.2f} MPa
    90th:   {model_df['psi50'].quantile(0.90):6.2f} MPa

Difference (Model - TRY):
  Mean difference:    {model_df['psi50'].mean() - p50_df['value'].mean():6.2f} MPa
  Median difference:  {model_df['psi50'].median() - p50_df['value'].median():6.2f} MPa
"""

fig.text(0.5, 0.17, stats_text, fontsize=9, family='monospace',
         verticalalignment='center', horizontalalignment='center',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, pad=1.0))

plt.suptitle('TRY P50 Database vs. Model Estimates', 
             fontsize=16, fontweight='bold', y=0.975)

# Add subtitle with counts
fig.text(0.5, 0.945, f'Model: {len(model_df):,} grid points  |  TRY: {len(p50_df):,} observations',
         ha='center', fontsize=12, style='italic', color='#333333')

plt.savefig('try_vs_model_p50_distribution.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: try_vs_model_p50_distribution.png")

# Save data summary
summary = pd.DataFrame({
    'Statistic': ['Count', 'Mean (MPa)', 'Median (MPa)', 'Std Dev (MPa)', 
                  'Min (MPa)', 'Max (MPa)', 'Q25 (MPa)', 'Q75 (MPa)'],
    'TRY_Observed': [
        len(p50_df),
        p50_df['value'].mean(),
        p50_df['value'].median(),
        p50_df['value'].std(),
        p50_df['value'].min(),
        p50_df['value'].max(),
        p50_df['value'].quantile(0.25),
        p50_df['value'].quantile(0.75)
    ],
    'Model_Estimated': [
        len(model_df),
        model_df['psi50'].mean(),
        model_df['psi50'].median(),
        model_df['psi50'].std(),
        model_df['psi50'].min(),
        model_df['psi50'].max(),
        model_df['psi50'].quantile(0.25),
        model_df['psi50'].quantile(0.75)
    ]
})

summary['Difference'] = summary['Model_Estimated'] - summary['TRY_Observed']
summary.to_csv('try_vs_model_p50_summary.csv', index=False)
print("✓ Saved: try_vs_model_p50_summary.csv")

print("\n" + "="*70)
print("Comparison complete!")
print("="*70)
print("\nNOTE: TRY data lacks coordinates, so no spatial comparison possible.")
print("To enable spatial comparison, request TRY data with Latitude/Longitude fields.")
print("="*70)
