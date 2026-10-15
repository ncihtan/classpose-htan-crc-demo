#!/usr/bin/env python3
"""
Enhanced visualizations for CRC cell density data.
Includes clustered heatmap, violin plots, and PCA.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Find all cell density CSV files
csv_pattern = "./*/*.csv"
csv_files = sorted(glob.glob(csv_pattern))

print(f"Found {len(csv_files)} CSV files")

# Load all data
all_data = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    sample_id = os.path.basename(os.path.dirname(csv_file))
    short_id = sample_id.split('_')[-1][:8]
    df['sample_id'] = short_id
    all_data.append(df)

combined_df = pd.concat(all_data, ignore_index=True)

# Create pivot tables
fractions_pivot = combined_df.pivot_table(
    index='sample_id', columns='cell_class', values='count', aggfunc='sum'
).fillna(0)
fractions_pivot = fractions_pivot.div(fractions_pivot.sum(axis=1), axis=0) * 100

density_pivot = combined_df.pivot_table(
    index='sample_id', columns='cell_class', values='density', aggfunc='sum'
).fillna(0)

# Reorder columns by mean abundance
col_order = fractions_pivot.mean().sort_values(ascending=False).index
fractions_pivot = fractions_pivot[col_order]
density_pivot = density_pivot[col_order]

# Define colors for cell types
cell_colors = {
    'Epithelial': '#E41A1C',
    'Connective': '#377EB8',
    'Lymphocyte': '#4DAF4A',
    'Plasma cell': '#984EA3',
    'Neutrophil': '#FF7F00',
    'Eosinophil': '#FFFF33'
}

# =============================================================================
# FIGURE 1: Clustered Heatmap
# =============================================================================
fig1, ax1 = plt.subplots(figsize=(12, 10))

# Z-score normalize for clustering
fractions_zscore = (fractions_pivot - fractions_pivot.mean()) / fractions_pivot.std()

# Create clustermap
g = sns.clustermap(
    fractions_zscore,
    cmap='RdBu_r',
    center=0,
    figsize=(10, 12),
    dendrogram_ratio=(0.15, 0.15),
    cbar_pos=(0.02, 0.8, 0.03, 0.15),
    yticklabels=True,
    xticklabels=True,
    row_cluster=True,
    col_cluster=True,
    linewidths=0.5,
    method='ward'
)
g.ax_heatmap.set_xlabel('Cell Type', fontsize=12)
g.ax_heatmap.set_ylabel('Sample ID', fontsize=12)
g.fig.suptitle('Clustered Heatmap: Cell Type Composition (Z-score)', y=1.02, fontsize=14)
plt.savefig('./cell_type_clustermap.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved: cell_type_clustermap.png")

# =============================================================================
# FIGURE 2: Violin/Box plots
# =============================================================================
fig2, axes = plt.subplots(1, 2, figsize=(14, 6))

# Melt data for seaborn
fractions_melted = fractions_pivot.reset_index().melt(
    id_vars='sample_id', var_name='Cell Type', value_name='Fraction (%)'
)
density_melted = density_pivot.reset_index().melt(
    id_vars='sample_id', var_name='Cell Type', value_name='Density (cells/mm²)'
)

# Violin plot - fractions
ax1 = axes[0]
sns.violinplot(data=fractions_melted, x='Cell Type', y='Fraction (%)',
               palette=cell_colors, ax=ax1, inner='box', cut=0)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
ax1.set_title('Cell Type Fractions Distribution', fontsize=12)
ax1.set_xlabel('')

# Violin plot - densities (log scale)
ax2 = axes[1]
density_melted['Log Density'] = np.log10(density_melted['Density (cells/mm²)'] + 1)
sns.violinplot(data=density_melted, x='Cell Type', y='Log Density',
               palette=cell_colors, ax=ax2, inner='box', cut=0)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
ax2.set_title('Cell Type Densities Distribution (log10)', fontsize=12)
ax2.set_ylabel('log₁₀(Density + 1)')
ax2.set_xlabel('')

plt.tight_layout()
plt.savefig('./cell_type_violins.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved: cell_type_violins.png")

# =============================================================================
# FIGURE 3: PCA with immune hot/cold annotation
# =============================================================================
fig3, axes = plt.subplots(1, 2, figsize=(14, 6))

# PCA on fractions
scaler = StandardScaler()
fractions_scaled = scaler.fit_transform(fractions_pivot)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(fractions_scaled)

# Color by lymphocyte fraction (immune status)
lymph_frac = fractions_pivot['Lymphocyte'].values

ax1 = axes[0]
scatter = ax1.scatter(pca_result[:, 0], pca_result[:, 1],
                      c=lymph_frac, cmap='RdYlGn', s=80, edgecolor='black', linewidth=0.5)
cbar = plt.colorbar(scatter, ax=ax1)
cbar.set_label('Lymphocyte %', fontsize=10)
ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=11)
ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=11)
ax1.set_title('PCA: Samples by Cell Composition\n(color = Lymphocyte fraction)', fontsize=12)
ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.3)

# Add loading vectors
ax2 = axes[1]
loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
for i, cell_type in enumerate(fractions_pivot.columns):
    ax2.arrow(0, 0, loadings[i, 0]*3, loadings[i, 1]*3,
              head_width=0.1, head_length=0.05, fc=cell_colors.get(cell_type, 'gray'),
              ec=cell_colors.get(cell_type, 'gray'), linewidth=2)
    ax2.text(loadings[i, 0]*3.3, loadings[i, 1]*3.3, cell_type,
             fontsize=10, ha='center', va='center', fontweight='bold')

ax2.set_xlim(-2, 2)
ax2.set_ylim(-2, 2)
ax2.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=11)
ax2.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=11)
ax2.set_title('PCA Loadings: Cell Type Contributions', fontsize=12)
ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
ax2.set_aspect('equal')

plt.tight_layout()
plt.savefig('./cell_type_pca.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved: cell_type_pca.png")

# =============================================================================
# FIGURE 4: Immune Score Summary - Dot plot style
# =============================================================================
fig4, ax = plt.subplots(figsize=(10, 8))

# Sort samples by lymphocyte fraction
sort_idx = fractions_pivot['Lymphocyte'].sort_values(ascending=True).index

# Create dot plot (size = fraction, color = log density)
for i, sample in enumerate(sort_idx):
    for j, cell_type in enumerate(col_order):
        frac = fractions_pivot.loc[sample, cell_type]
        dens = np.log10(density_pivot.loc[sample, cell_type] + 1)
        ax.scatter(j, i, s=frac*10, c=[dens], cmap='viridis',
                   vmin=0, vmax=4, edgecolor='black', linewidth=0.3)

ax.set_xticks(range(len(col_order)))
ax.set_xticklabels(col_order, rotation=45, ha='right', fontsize=10)
ax.set_yticks(range(len(sort_idx)))
ax.set_yticklabels(sort_idx, fontsize=7)
ax.set_xlabel('Cell Type', fontsize=12)
ax.set_ylabel('Sample (sorted by Lymphocyte %: cold → hot)', fontsize=12)
ax.set_title('Cell Composition Dot Plot\n(size = fraction %, color = log₁₀ density)', fontsize=13)

# Add colorbar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, 4))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.5)
cbar.set_label('log₁₀(density + 1)', fontsize=10)

# Add size legend
for s in [5, 20, 50]:
    ax.scatter([], [], s=s*10, c='gray', edgecolor='black', linewidth=0.3, label=f'{s}%')
ax.legend(title='Fraction', loc='upper left', bbox_to_anchor=(1.15, 0.5), frameon=True)

plt.tight_layout()
plt.savefig('./cell_type_dotplot.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved: cell_type_dotplot.png")

# =============================================================================
# FIGURE 5: Immune infiltration summary bar
# =============================================================================
fig5, ax = plt.subplots(figsize=(12, 4))

# Calculate total immune fraction (Lymphocyte + Plasma cell + Neutrophil + Eosinophil)
immune_cells = ['Lymphocyte', 'Plasma cell', 'Neutrophil', 'Eosinophil']
available_immune = [c for c in immune_cells if c in fractions_pivot.columns]
immune_fraction = fractions_pivot[available_immune].sum(axis=1).sort_values(ascending=False)

colors_bar = ['#2ca02c' if v > 15 else '#ff7f0e' if v > 8 else '#d62728' for v in immune_fraction]

bars = ax.bar(range(len(immune_fraction)), immune_fraction, color=colors_bar, edgecolor='black', linewidth=0.5)
ax.axhline(y=15, color='green', linestyle='--', alpha=0.7, label='Immune Hot (>15%)')
ax.axhline(y=8, color='orange', linestyle='--', alpha=0.7, label='Intermediate (>8%)')
ax.set_xticks(range(len(immune_fraction)))
ax.set_xticklabels(immune_fraction.index, rotation=90, fontsize=7)
ax.set_xlabel('Sample ID', fontsize=11)
ax.set_ylabel('Total Immune Fraction (%)', fontsize=11)
ax.set_title('Tumor Immune Infiltration Score\n(Lymphocytes + Plasma cells + Granulocytes)', fontsize=13)
ax.legend(loc='upper right')

# Add immune status annotation
n_hot = (immune_fraction > 15).sum()
n_inter = ((immune_fraction > 8) & (immune_fraction <= 15)).sum()
n_cold = (immune_fraction <= 8).sum()
ax.text(0.02, 0.95, f'Hot: {n_hot}  |  Intermediate: {n_inter}  |  Cold: {n_cold}',
        transform=ax.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('./immune_infiltration_score.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved: immune_infiltration_score.png")

print("\nAll visualizations complete!")
