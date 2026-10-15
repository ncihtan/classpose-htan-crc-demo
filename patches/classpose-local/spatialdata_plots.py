#!/usr/bin/env python
"""Generate useful plots from ClassPose SpatialData output."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
from pathlib import Path
from scipy import ndimage
from matplotlib.patches import Patch
import openslide

# Paths
sdata_path = Path("output/CMU-1-Small-Region_spatialdata.zarr")
slide_path = Path("test_data/CMU-1-Small-Region.svs")
output_dir = Path("output")

# Load spatialdata
print(f"Loading SpatialData from {sdata_path}")
sdata = sd.read_zarr(sdata_path)
print(f"SpatialData contents: {sdata}")

# Get metadata
metadata = sdata.attrs.get("metadata", {})
mpp_x = metadata.get("mpp_x", 0.499)
mpp_y = metadata.get("mpp_y", 0.499)
print(f"MPP: ({mpp_x}, {mpp_y})")

# Extract cell centroids as pandas DataFrame
cell_centroids = sdata["cell_centroids"].compute()
print(f"Cell centroids shape: {cell_centroids.shape}")
print(f"Columns: {cell_centroids.columns.tolist()}")

# Get cell contours for plotting
cell_contours = sdata["cell_contours"]
print(f"Cell contours: {len(cell_contours)} polygons")

# Load slide for background
slide = openslide.OpenSlide(str(slide_path))
dims = slide.dimensions
thumbnail = slide.get_thumbnail((800, int(800 * dims[1] / dims[0])))
thumbnail_arr = np.array(thumbnail)

# Define colors for cell types
cell_type_colors = {
    'Neutrophil': '#8dd3c7',
    'Epithelial': '#ffffb3',
    'Lymphocyte': '#bebada',
    'Plasma cell': '#fb8072',
    'Eosinophil': '#80b1d3',
    'Connective': '#fdb462',
}

# Create figure with multiple subplots
fig = plt.figure(figsize=(16, 12))

# ============================================================================
# Plot 1: Cell type distribution (pie chart)
# ============================================================================
ax1 = fig.add_subplot(2, 3, 1)

cell_counts = cell_centroids['classification_name'].value_counts()
colors = [cell_type_colors.get(ct, '#999999') for ct in cell_counts.index]

wedges, texts, autotexts = ax1.pie(
    cell_counts.values,
    labels=cell_counts.index,
    colors=colors,
    autopct='%1.1f%%',
    pctdistance=0.75
)
ax1.set_title('Cell Type Distribution', fontsize=12, fontweight='bold')

# ============================================================================
# Plot 2: Cell counts bar chart
# ============================================================================
ax2 = fig.add_subplot(2, 3, 2)

bars = ax2.bar(cell_counts.index, cell_counts.values,
               color=[cell_type_colors.get(ct, '#999999') for ct in cell_counts.index])
ax2.set_ylabel('Count')
ax2.set_title('Cell Counts by Type', fontsize=12, fontweight='bold')
ax2.tick_params(axis='x', rotation=45)

for bar, count in zip(bars, cell_counts.values):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
             str(count), ha='center', va='bottom', fontsize=9)

# ============================================================================
# Plot 3: Spatial scatter plot of all cells
# ============================================================================
ax3 = fig.add_subplot(2, 3, 3)

ax3.imshow(thumbnail_arr, extent=[0, dims[0], dims[1], 0], alpha=0.5)

for cell_type in cell_centroids['classification_name'].unique():
    mask = cell_centroids['classification_name'] == cell_type
    ax3.scatter(
        cell_centroids.loc[mask, 'x'],
        cell_centroids.loc[mask, 'y'],
        c=cell_type_colors.get(cell_type, '#999999'),
        label=cell_type,
        s=5,
        alpha=0.7
    )

ax3.set_xlim(0, dims[0])
ax3.set_ylim(dims[1], 0)
ax3.set_aspect('equal')
ax3.set_title('Spatial Distribution of Cells', fontsize=12, fontweight='bold')
ax3.legend(loc='upper right', fontsize=8, markerscale=2)
ax3.axis('off')

# ============================================================================
# Plot 4: Cell density heatmap (all cells)
# ============================================================================
ax4 = fig.add_subplot(2, 3, 4)

# Create 2D histogram for density
bin_size = 100  # pixels per bin
x_bins = np.arange(0, dims[0] + bin_size, bin_size)
y_bins = np.arange(0, dims[1] + bin_size, bin_size)

density, x_edges, y_edges = np.histogram2d(
    cell_centroids['x'],
    cell_centroids['y'],
    bins=[x_bins, y_bins]
)

# Smooth the density map
density_smooth = ndimage.gaussian_filter(density.T, sigma=1.5)

ax4.imshow(thumbnail_arr, extent=[0, dims[0], dims[1], 0], alpha=0.3)
im = ax4.imshow(
    density_smooth,
    extent=[0, dims[0], dims[1], 0],
    cmap='hot',
    alpha=0.7,
    interpolation='bilinear'
)
plt.colorbar(im, ax=ax4, label='Cell density', shrink=0.8)
ax4.set_title('Cell Density Heatmap', fontsize=12, fontweight='bold')
ax4.set_xlim(0, dims[0])
ax4.set_ylim(dims[1], 0)
ax4.axis('off')

# ============================================================================
# Plot 5: Lymphocyte vs Epithelial spatial comparison
# ============================================================================
ax5 = fig.add_subplot(2, 3, 5)

ax5.imshow(thumbnail_arr, extent=[0, dims[0], dims[1], 0], alpha=0.4)

# Plot epithelial cells
epi_mask = cell_centroids['classification_name'] == 'Epithelial'
ax5.scatter(
    cell_centroids.loc[epi_mask, 'x'],
    cell_centroids.loc[epi_mask, 'y'],
    c=cell_type_colors['Epithelial'],
    label='Epithelial',
    s=8,
    alpha=0.6,
    edgecolors='orange',
    linewidths=0.3
)

# Plot lymphocytes
lymph_mask = cell_centroids['classification_name'] == 'Lymphocyte'
ax5.scatter(
    cell_centroids.loc[lymph_mask, 'x'],
    cell_centroids.loc[lymph_mask, 'y'],
    c=cell_type_colors['Lymphocyte'],
    label='Lymphocyte',
    s=12,
    alpha=0.8,
    edgecolors='purple',
    linewidths=0.5
)

ax5.set_xlim(0, dims[0])
ax5.set_ylim(dims[1], 0)
ax5.set_aspect('equal')
ax5.set_title('Epithelial vs Lymphocyte Distribution', fontsize=12, fontweight='bold')
ax5.legend(loc='upper right', fontsize=9, markerscale=1.5)
ax5.axis('off')

# ============================================================================
# Plot 6: Cell area distribution by type
# ============================================================================
ax6 = fig.add_subplot(2, 3, 6)

# Get area data for each cell type
area_data = []
labels = []
colors_box = []

for cell_type in cell_counts.index:
    mask = cell_centroids['classification_name'] == cell_type
    if 'area' in cell_centroids.columns:
        areas = cell_centroids.loc[mask, 'area'].values
        if len(areas) > 0:
            area_data.append(areas)
            labels.append(cell_type)
            colors_box.append(cell_type_colors.get(cell_type, '#999999'))

if area_data:
    bp = ax6.boxplot(area_data, labels=labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax6.set_ylabel('Area (pixels²)')
    ax6.set_title('Cell Area Distribution by Type', fontsize=12, fontweight='bold')
    ax6.tick_params(axis='x', rotation=45)
else:
    ax6.text(0.5, 0.5, 'Area data not available', ha='center', va='center', transform=ax6.transAxes)
    ax6.set_title('Cell Area Distribution by Type', fontsize=12, fontweight='bold')

plt.suptitle('ClassPose SpatialData Analysis\nCMU-1-Small-Region.svs',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = output_dir / "spatialdata_analysis.png"
plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\nSaved analysis plots to: {output_path}")

# ============================================================================
# Additional plot: Density comparison by cell type
# ============================================================================
fig2, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

cell_types = ['Epithelial', 'Connective', 'Lymphocyte', 'Plasma cell', 'Neutrophil', 'Eosinophil']

for idx, cell_type in enumerate(cell_types):
    ax = axes[idx]

    mask = cell_centroids['classification_name'] == cell_type
    type_cells = cell_centroids[mask]

    ax.imshow(thumbnail_arr, extent=[0, dims[0], dims[1], 0], alpha=0.4)

    if len(type_cells) > 0:
        # Create density map for this cell type
        density, _, _ = np.histogram2d(
            type_cells['x'],
            type_cells['y'],
            bins=[x_bins, y_bins]
        )
        density_smooth = ndimage.gaussian_filter(density.T, sigma=1.5)

        im = ax.imshow(
            density_smooth,
            extent=[0, dims[0], dims[1], 0],
            cmap='YlOrRd',
            alpha=0.7,
            interpolation='bilinear'
        )
        plt.colorbar(im, ax=ax, shrink=0.6)

    ax.set_title(f'{cell_type}\n(n={len(type_cells)})', fontsize=11, fontweight='bold',
                color=cell_type_colors.get(cell_type, 'black'))
    ax.set_xlim(0, dims[0])
    ax.set_ylim(dims[1], 0)
    ax.axis('off')

plt.suptitle('Cell Density by Type', fontsize=14, fontweight='bold')
plt.tight_layout()

output_path2 = output_dir / "spatialdata_density_by_type.png"
plt.savefig(output_path2, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved density by type plots to: {output_path2}")

plt.show()
