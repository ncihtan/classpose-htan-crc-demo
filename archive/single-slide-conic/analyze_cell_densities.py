#!/usr/bin/env python3
"""
Analyze cell density CSV files and create stacked bar plots showing
cell type fractions and densities per sample/tissue.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Find all cell density CSV files
csv_pattern = "./*/*.csv"
csv_files = sorted(glob.glob(csv_pattern))

print(f"Found {len(csv_files)} CSV files")

# Load all data
all_data = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    # Extract sample ID from directory name
    sample_id = os.path.basename(os.path.dirname(csv_file))
    # Shorten the sample ID for display (take last 8 chars of UUID)
    short_id = sample_id.split('_')[-1][:8]
    df['sample_id'] = short_id
    df['full_sample_id'] = sample_id
    all_data.append(df)

# Combine all data
combined_df = pd.concat(all_data, ignore_index=True)
print(f"Total rows: {len(combined_df)}")
print(f"Unique samples: {combined_df['sample_id'].nunique()}")
print(f"Cell types: {combined_df['cell_class'].unique()}")

# Pivot data for plotting - counts
counts_pivot = combined_df.pivot_table(
    index='sample_id',
    columns='cell_class',
    values='count',
    aggfunc='sum'
).fillna(0)

# Pivot data for plotting - densities
density_pivot = combined_df.pivot_table(
    index='sample_id',
    columns='cell_class',
    values='density',
    aggfunc='sum'
).fillna(0)

# Calculate fractions from counts
fractions_pivot = counts_pivot.div(counts_pivot.sum(axis=1), axis=0) * 100

# Calculate density fractions (percentage of total density)
density_fractions_pivot = density_pivot.div(density_pivot.sum(axis=1), axis=0) * 100

# Sort all pivots by lymphocyte fraction (immune hot to cold)
sort_order = fractions_pivot['Lymphocyte'].sort_values(ascending=False).index
counts_pivot = counts_pivot.loc[sort_order]
density_pivot = density_pivot.loc[sort_order]
fractions_pivot = fractions_pivot.loc[sort_order]
density_fractions_pivot = density_fractions_pivot.loc[sort_order]

# Define a color palette for cell types
cell_types = counts_pivot.columns.tolist()
colors = plt.cm.Set2(np.linspace(0, 1, len(cell_types)))
color_dict = dict(zip(cell_types, colors))

# Create figure with 4 subplots
fig, axes = plt.subplots(2, 2, figsize=(20, 14))

# Plot 1: Stacked bar - Cell Counts (absolute)
ax1 = axes[0, 0]
bottom = np.zeros(len(counts_pivot))
for cell_type in cell_types:
    ax1.bar(counts_pivot.index, counts_pivot[cell_type], bottom=bottom,
            label=cell_type, color=color_dict[cell_type])
    bottom += counts_pivot[cell_type].values
ax1.set_xlabel('Sample ID')
ax1.set_ylabel('Cell Count')
ax1.set_title('Cell Counts by Type (Absolute)')
ax1.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax1.tick_params(axis='x', rotation=90)

# Plot 2: Stacked bar - Cell Fractions (100% scaled)
ax2 = axes[0, 1]
bottom = np.zeros(len(fractions_pivot))
for cell_type in cell_types:
    ax2.bar(fractions_pivot.index, fractions_pivot[cell_type], bottom=bottom,
            label=cell_type, color=color_dict[cell_type])
    bottom += fractions_pivot[cell_type].values
ax2.set_xlabel('Sample ID')
ax2.set_ylabel('Percentage (%)')
ax2.set_title('Cell Fractions by Type (100% Scaled)')
ax2.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax2.tick_params(axis='x', rotation=90)
ax2.set_ylim(0, 100)

# Plot 3: Stacked bar - Cell Densities (absolute)
ax3 = axes[1, 0]
bottom = np.zeros(len(density_pivot))
for cell_type in cell_types:
    ax3.bar(density_pivot.index, density_pivot[cell_type], bottom=bottom,
            label=cell_type, color=color_dict[cell_type])
    bottom += density_pivot[cell_type].values
ax3.set_xlabel('Sample ID')
ax3.set_ylabel('Cell Density')
ax3.set_title('Cell Densities by Type (Absolute)')
ax3.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax3.tick_params(axis='x', rotation=90)

# Plot 4: Stacked bar - Cell Density Fractions (100% scaled)
ax4 = axes[1, 1]
bottom = np.zeros(len(density_fractions_pivot))
for cell_type in cell_types:
    ax4.bar(density_fractions_pivot.index, density_fractions_pivot[cell_type], bottom=bottom,
            label=cell_type, color=color_dict[cell_type])
    bottom += density_fractions_pivot[cell_type].values
ax4.set_xlabel('Sample ID')
ax4.set_ylabel('Percentage (%)')
ax4.set_title('Cell Density Fractions by Type (100% Scaled)')
ax4.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
ax4.tick_params(axis='x', rotation=90)
ax4.set_ylim(0, 100)

plt.tight_layout()
plt.savefig('./cell_type_stacked_bars.png',
            dpi=150, bbox_inches='tight')
plt.savefig('./cell_type_stacked_bars.pdf',
            bbox_inches='tight')
print("Saved: cell_type_stacked_bars.png and cell_type_stacked_bars.pdf")

# Also create summary statistics
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)

print("\nMean cell counts per sample:")
print(counts_pivot.mean().sort_values(ascending=False))

print("\nMean cell densities per sample:")
print(density_pivot.mean().sort_values(ascending=False))

print("\nMean cell fractions (%):")
print(fractions_pivot.mean().sort_values(ascending=False))

plt.show()
