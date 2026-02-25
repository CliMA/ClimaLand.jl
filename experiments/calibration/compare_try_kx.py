#!/usr/bin/env python3
"""
Compare TRY hydraulic conductivity distribution with model estimates

Since the TRY data lacks coordinates, we'll compare:
1. Overall distributions
2. Statistical moments
3. Range comparison

NOTE: TRY and model may use different units/definitions of hydraulic conductivity
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

print("="*70)
print("TRY Hydraulic Conductivity vs Model Estimates")
print("="*70)

# Load TRY data
print("\nLoading TRY data...")
try_df = pd.read_csv('TRY_data.txt', sep='\t', low_memory=False, encoding='ISO-8859-1')

# Extract hydraulic conductivity/conductance
kx_trait_names = [
    'Branch hydraulic conductance',
    'Plant hydraulic conductance',
    'Leaf hydraulic conductance',
    'Shoot: plant above ground hydraulic conductance'
]

kx_df = try_df[try_df['TraitName'].isin(kx_trait_names)]

# Use OrigValueStr since StdValue is not populated for these traits
kx_df['value'] = pd.to_numeric(kx_df['OrigValueStr'], errors='coerce')
kx_df = kx_df.dropna(subset=['value'])

# Filter positive values only (conductivity should be positive)
kx_df = kx_df[kx_df['value'] > 0]

# Log transform for better visualization (conductivity spans many orders of magnitude)
kx_df['log_value'] = np.log10(kx_df['value'])

print(f"  TRY hydraulic conductivity observations: {len(kx_df)}")
print(f"  Range: {kx_df['value'].min():.2e} to {kx_df['value'].max():.2e}")
print(f"  Log10 range: {kx_df['log_value'].min():.2f} to {kx_df['log_value'].max():.2f}")
print(f"  Median: {kx_df['value'].median():.2e}")

if 'OrigUnitStr' in kx_df.columns:
    print(f"  Units in TRY data:")
    print(f"    {kx_df['OrigUnitStr'].value_counts().head()}")

# Load model data
print("\nLoading model estimates...")
model_df = pd.read_csv('traits_kx_data.csv')

# Model kx is dimensionless - take log for comparison
model_df['log_kx'] = np.log10(model_df['kx'])

print(f"  Model grid points: {len(model_df)}")
print(f"  Range: {model_df['kx'].min():.2f} to {model_df['kx'].max():.2f}")
print(f"  Log10 range: {model_df['log_kx'].min():.2f} to {model_df['log_kx'].max():.2f}")
print(f"  Median: {model_df['kx'].median():.2f}")

# Create comparison figure
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3, bottom=0.35, top=0.92)

# Panel 1: Overlapping histograms (log scale)
ax1 = fig.add_subplot(gs[0, :])
bins = np.linspace(-6, 3, 40)

ax1.hist(kx_df['log_value'], bins=bins, alpha=0.6, 
         label=f'TRY observations (n={len(kx_df)})', 
         color='steelblue', edgecolor='navy', linewidth=0.5, density=True)
ax1.hist(model_df['log_kx'], bins=bins, alpha=0.6, 
         label=f'Model estimates (n={len(model_df)})', 
         color='coral', edgecolor='darkred', linewidth=0.5, density=True)

# Add median lines
ax1.axvline(kx_df['log_value'].median(), color='navy', linestyle='--', linewidth=2,
            label=f'TRY median: {kx_df["log_value"].median():.2f}')
ax1.axvline(model_df['log_kx'].median(), color='darkred', linestyle='--', linewidth=2,
            label=f'Model median: {model_df["log_kx"].median():.2f}')

ax1.set_xlabel('log₁₀(Hydraulic Conductivity)', fontsize=13)
ax1.set_ylabel('Probability Density', fontsize=13)
ax1.set_title('Hydraulic Conductivity Distribution Comparison (Log Scale)', fontsize=15, fontweight='bold')
ax1.legend(fontsize=11, framealpha=0.9)
ax1.grid(alpha=0.3, axis='y')

# Panel 2: Cumulative distributions
ax2 = fig.add_subplot(gs[1, 0])

try_sorted = np.sort(kx_df['log_value'].values)
try_cumsum = np.arange(1, len(try_sorted)+1) / len(try_sorted)

model_sorted = np.sort(model_df['log_kx'].values)
model_cumsum = np.arange(1, len(model_sorted)+1) / len(model_sorted)

ax2.plot(try_sorted, try_cumsum, linewidth=2.5, color='steelblue', label='TRY observations')
ax2.plot(model_sorted, model_cumsum, linewidth=2.5, color='coral', label='Model estimates')

ax2.set_xlabel('log₁₀(Hydraulic Conductivity)', fontsize=12)
ax2.set_ylabel('Cumulative Probability', fontsize=12)
ax2.set_title('Cumulative Distribution Function', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

# Panel 3: Box plots
ax3 = fig.add_subplot(gs[1, 1])

box_data = [kx_df['log_value'].values, model_df['log_kx'].values]
bp = ax3.boxplot(box_data, tick_labels=['TRY\nObservations', 'Model\nEstimates'],
                 patch_artist=True, widths=0.6,
                 medianprops=dict(color='black', linewidth=2),
                 boxprops=dict(facecolor='lightblue', edgecolor='navy', linewidth=1.5),
                 whiskerprops=dict(color='navy', linewidth=1.5),
                 capprops=dict(color='navy', linewidth=1.5))

# Color the boxes differently
bp['boxes'][0].set_facecolor('steelblue')
bp['boxes'][1].set_facecolor('coral')

ax3.set_ylabel('log₁₀(Hydraulic Conductivity)', fontsize=12)
ax3.set_title('Distribution Statistics', fontsize=14, fontweight='bold')
ax3.grid(alpha=0.3, axis='y')

# Add text with statistics
stats_text = f"""
COMPARISON SUMMARY (Log₁₀ Scale)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TRY Database (n={len(kx_df)}):
  Mean:     {kx_df['log_value'].mean():6.2f}
  Median:   {kx_df['log_value'].median():6.2f}
  Std Dev:  {kx_df['log_value'].std():6.2f}
  Range:    [{kx_df['log_value'].min():.2f}, {kx_df['log_value'].max():.2f}]
  
  Percentiles:
    10th:   {kx_df['log_value'].quantile(0.10):6.2f}
    25th:   {kx_df['log_value'].quantile(0.25):6.2f}
    75th:   {kx_df['log_value'].quantile(0.75):6.2f}
    90th:   {kx_df['log_value'].quantile(0.90):6.2f}

Model Estimates (n={len(model_df)}):
  Mean:     {model_df['log_kx'].mean():6.2f}
  Median:   {model_df['log_kx'].median():6.2f}
  Std Dev:  {model_df['log_kx'].std():6.2f}
  Range:    [{model_df['log_kx'].min():.2f}, {model_df['log_kx'].max():.2f}]
  
  Percentiles:
    10th:   {model_df['log_kx'].quantile(0.10):6.2f}
    25th:   {model_df['log_kx'].quantile(0.25):6.2f}
    75th:   {model_df['log_kx'].quantile(0.75):6.2f}
    90th:   {model_df['log_kx'].quantile(0.90):6.2f}

Difference (Model - TRY):
  Mean difference:    {model_df['log_kx'].mean() - kx_df['log_value'].mean():6.2f}
  Median difference:  {model_df['log_kx'].median() - kx_df['log_value'].median():6.2f}

NOTE: TRY and model may use different units/definitions.
      Direct quantitative comparison should be interpreted cautiously.
      Comparison focuses on distributional patterns.
"""

fig.text(0.5, 0.17, stats_text, fontsize=9, family='monospace',
         verticalalignment='center', horizontalalignment='center',
         bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8, pad=1.0))

plt.suptitle('TRY Hydraulic Conductivity vs. Model Estimates', 
             fontsize=16, fontweight='bold', y=0.975)

# Add subtitle with counts
fig.text(0.5, 0.945, f'Model: {len(model_df):,} grid points  |  TRY: {len(kx_df):,} observations',
         ha='center', fontsize=12, style='italic', color='#333333')

plt.savefig('try_vs_model_kx_distribution.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: try_vs_model_kx_distribution.png")

# Save data summary
summary = pd.DataFrame({
    'Statistic': ['Count', 'Mean (log10)', 'Median (log10)', 'Std Dev (log10)', 
                  'Min (log10)', 'Max (log10)', 'Q25 (log10)', 'Q75 (log10)'],
    'TRY_Observed': [
        len(kx_df),
        kx_df['log_value'].mean(),
        kx_df['log_value'].median(),
        kx_df['log_value'].std(),
        kx_df['log_value'].min(),
        kx_df['log_value'].max(),
        kx_df['log_value'].quantile(0.25),
        kx_df['log_value'].quantile(0.75)
    ],
    'Model_Estimated': [
        len(model_df),
        model_df['log_kx'].mean(),
        model_df['log_kx'].median(),
        model_df['log_kx'].std(),
        model_df['log_kx'].min(),
        model_df['log_kx'].max(),
        model_df['log_kx'].quantile(0.25),
        model_df['log_kx'].quantile(0.75)
    ]
})

summary['Difference'] = summary['Model_Estimated'] - summary['TRY_Observed']
summary.to_csv('try_vs_model_kx_summary.csv', index=False)
print("✓ Saved: try_vs_model_kx_summary.csv")

print("\n" + "="*70)
print("Comparison complete!")
print("="*70)
print("\nNOTE: TRY data lacks coordinates, so no spatial comparison possible.")
print("TRY and model hydraulic conductivity may use different units/definitions.")
print("Comparison focuses on distributional patterns rather than absolute values.")
print("="*70)
