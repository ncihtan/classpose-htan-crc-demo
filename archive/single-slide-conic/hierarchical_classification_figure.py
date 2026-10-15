"""
Publication-quality figure for hierarchical cell classification
Combines CoNIC (immune-focused) and ConSep (epithelial-focused) models
"""

import spatialdata as sd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
import openslide

# Configuration
BASE_DIR = "."
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
UM_PER_PIXEL = 0.325
FOV_UM = 300  # Smaller FOV for detail
FOV_PX = int(FOV_UM / UM_PER_PIXEL)

# Paths
conic_zarr_path = f"{BASE_DIR}/{SAMPLE_ID}/{SAMPLE_ID}_spatialdata.zarr"
consep_zarr_path = f"{BASE_DIR}/dual_model_output/consep/{SAMPLE_ID}_consep_spatialdata.zarr"
he_slide_path = f"{BASE_DIR}/{SAMPLE_ID}.tif"

# Color palettes - carefully chosen for publication
conic_colors = {
    'Epithelial': '#4363d8',     # blue
    'Connective': '#f58231',     # orange
    'Lymphocyte': '#3cb44b',     # green
    'Plasma cell': '#e6194B',    # red
    'Eosinophil': '#911eb4',     # purple
    'Neutrophil': '#9A6324',     # brown
}

consep_colors = {
    'Healthy epithelial': '#42d4f4',   # cyan
    'Malignant epithelial': '#e6194B', # red
    'Inflammatory': '#3cb44b',          # green
    'Stroma': '#f58231',                # orange
    'Muscle': '#911eb4',                # purple
    'Other': '#a9a9a9',                 # gray
}

hierarchical_colors = {
    'Malignant Epithelial': '#e6194B',    # red
    'Healthy Epithelial': '#42d4f4',       # cyan
    'Epithelial (other)': '#aaffc3',       # mint
    'Lymphocyte': '#3cb44b',               # green
    'Plasma cell': '#f032e6',              # magenta
    'Eosinophil': '#911eb4',               # purple
    'Neutrophil': '#9A6324',               # brown
    'Stroma': '#f58231',                   # orange
    'Stroma (other)': '#ffd8b1',           # apricot
    'Muscle': '#800000',                   # maroon
    'Inflammatory (untyped)': '#bfef45',   # lime
    'Other': '#a9a9a9',                    # gray
}


def add_scale_bar(ax, length_um, pixel_scale, position='lower right', color='black', fontsize=8):
    """Add a scale bar to the axis."""
    length_px = length_um / pixel_scale

    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Calculate position (account for inverted y-axis)
    if 'lower' in position:
        y_pos = max(ylim) - (max(ylim) - min(ylim)) * 0.08
    else:
        y_pos = min(ylim) + (max(ylim) - min(ylim)) * 0.08

    if 'right' in position:
        x_pos = max(xlim) - (max(xlim) - min(xlim)) * 0.05 - length_px
    else:
        x_pos = min(xlim) + (max(xlim) - min(xlim)) * 0.05

    # Draw scale bar
    ax.plot([x_pos, x_pos + length_px], [y_pos, y_pos], color=color, linewidth=2, solid_capstyle='butt')

    # Add label
    ax.text(x_pos + length_px/2, y_pos - (max(ylim) - min(ylim)) * 0.03,
            f'{length_um} µm', ha='center', va='top', fontsize=fontsize, color=color, fontweight='bold')


def plot_contours(ax, contours, cell_type_colors, x_min, x_max, y_min, y_max):
    """Plot cell contours for a region."""
    centroids = contours.geometry.centroid
    mask = (
        (centroids.x >= x_min) & (centroids.x <= x_max) &
        (centroids.y >= y_min) & (centroids.y <= y_max)
    )
    contours_in_region = contours[mask]

    cell_counts = {}
    for cell_type, color in cell_type_colors.items():
        ct_contours = contours_in_region[contours_in_region['classification_name'] == cell_type]
        n_cells = len(ct_contours)
        cell_counts[cell_type] = n_cells

        if n_cells > 0:
            patches = []
            for geom in ct_contours.geometry:
                if geom.geom_type == 'Polygon':
                    coords = np.array(geom.exterior.coords)
                    patches.append(Polygon(coords, closed=True))
                elif geom.geom_type == 'MultiPolygon':
                    for poly in geom.geoms:
                        coords = np.array(poly.exterior.coords)
                        patches.append(Polygon(coords, closed=True))

            if patches:
                collection = PatchCollection(
                    patches, facecolor='none', edgecolor=color,
                    alpha=0.9, linewidth=0.6
                )
                ax.add_collection(collection)

    return cell_counts, len(contours_in_region)


def plot_hierarchical_points(ax, matched_df, x_min, x_max, y_min, y_max):
    """Plot hierarchical classification points for a region."""
    mask = (
        (matched_df['x'] >= x_min) & (matched_df['x'] <= x_max) &
        (matched_df['y'] >= y_min) & (matched_df['y'] <= y_max)
    )
    region_df = matched_df[mask]

    cell_counts = {}
    for cell_class in hierarchical_colors.keys():
        subset = region_df[region_df['simplified_class'] == cell_class]
        n_cells = len(subset)
        cell_counts[cell_class] = n_cells

        if n_cells > 0:
            ax.scatter(
                subset['x'], subset['y'],
                c=hierarchical_colors[cell_class],
                s=8, alpha=0.8, edgecolors='none'
            )

    return cell_counts, len(region_df)


def create_legend(ax, colors_dict, counts_dict, title, ncol=1, fontsize=7, loc='upper right'):
    """Create a clean legend with counts."""
    handles = []
    labels = []

    # Sort by count descending
    sorted_items = sorted(counts_dict.items(), key=lambda x: x[1], reverse=True)

    for cell_type, count in sorted_items:
        if count > 0 and cell_type in colors_dict:
            handles.append(plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=colors_dict[cell_type], markersize=6))
            pct = 100 * count / sum(counts_dict.values()) if sum(counts_dict.values()) > 0 else 0
            labels.append(f'{cell_type} ({pct:.0f}%)')

    legend = ax.legend(handles, labels, loc=loc, fontsize=fontsize,
                       framealpha=0.9, title=title, title_fontsize=fontsize+1,
                       ncol=ncol, handletextpad=0.3, columnspacing=0.5)
    legend.get_frame().set_linewidth(0.5)
    return legend


# Load data
print("Loading data...")
he_slide = openslide.OpenSlide(he_slide_path)
slide_dims = he_slide.dimensions

sdata_conic = sd.read_zarr(conic_zarr_path)
contours_conic = sdata_conic.shapes['cell_contours']

sdata_consep = sd.read_zarr(consep_zarr_path)
contours_consep = sdata_consep.shapes['cell_contours']

# Load or create hierarchical classification
try:
    matched = pd.read_parquet(f"{BASE_DIR}/hierarchical_cell_classification.parquet")
    print(f"Loaded hierarchical classification: {len(matched):,} cells")
except FileNotFoundError:
    print("Creating hierarchical classification...")
    # [Include matching code here if needed]
    raise FileNotFoundError("Run hierarchical classification first")

# Define regions of interest
regions = [
    (58000, 18000, "Immune-rich"),
    (15000, 38000, "Glandular"),
]

print("Creating figure...")

# Create figure with GridSpec for precise layout
fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(3, 4, height_ratios=[1.2, 1, 1], width_ratios=[1, 1, 1, 0.5],
                      hspace=0.08, wspace=0.08)

# Point size for overview panels (very small for 1M+ cells)
OVERVIEW_POINT_SIZE = 0.0001  # 10x smaller than before

# Row 1: Overview panels (CoNIC, ConSep, Hierarchical)
# Get H&E thumbnail
thumbnail_size = (1500, int(1500 * slide_dims[1] / slide_dims[0]))
he_thumbnail = he_slide.get_thumbnail(thumbnail_size)
he_thumbnail_arr = np.array(he_thumbnail.convert('RGB'))

# Panel A: CoNIC overview
ax_conic_ov = fig.add_subplot(gs[0, 0])
ax_conic_ov.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=1.0)

points_conic = sdata_conic.points["cell_centroids"].compute()
conic_counts_ov = {}
for cell_type, color in conic_colors.items():
    mask = points_conic["classification_name"] == cell_type
    subset = points_conic[mask]
    conic_counts_ov[cell_type] = len(subset)
    if len(subset) > 0:
        ax_conic_ov.scatter(subset["x"], subset["y"], c=color, s=OVERVIEW_POINT_SIZE, alpha=0.7, rasterized=True)

# Draw region boxes (yellow for visibility)
for idx, (cx, cy, _) in enumerate(regions):
    half = FOV_PX // 2
    rect = Rectangle((cx - half, cy - half), FOV_PX, FOV_PX,
                     fill=False, edgecolor='#ffff00', linewidth=2.5)
    ax_conic_ov.add_patch(rect)
    ax_conic_ov.text(cx, cy - half - 600, f'{idx+1}', color='yellow', fontsize=9,
                    ha='center', va='bottom', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='black', alpha=0.6, edgecolor='none'))

ax_conic_ov.set_xlim(0, slide_dims[0])
ax_conic_ov.set_ylim(slide_dims[1], 0)
ax_conic_ov.axis('off')
ax_conic_ov.set_title('CoNIC\n(Immune-focused)', fontsize=11, fontweight='bold', pad=5)
ax_conic_ov.text(-0.05, 1.05, 'A', transform=ax_conic_ov.transAxes, fontsize=14, fontweight='bold', va='top')
add_scale_bar(ax_conic_ov, 5000, UM_PER_PIXEL, position='lower right', color='white', fontsize=8)

# Panel B: ConSep overview
ax_consep_ov = fig.add_subplot(gs[0, 1])
ax_consep_ov.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=1.0)

points_consep = sdata_consep.points["cell_centroids"].compute()
consep_counts_ov = {}
for cell_type, color in consep_colors.items():
    mask = points_consep["classification_name"] == cell_type
    subset = points_consep[mask]
    consep_counts_ov[cell_type] = len(subset)
    if len(subset) > 0:
        ax_consep_ov.scatter(subset["x"], subset["y"], c=color, s=OVERVIEW_POINT_SIZE, alpha=0.7, rasterized=True)

for idx, (cx, cy, _) in enumerate(regions):
    half = FOV_PX // 2
    rect = Rectangle((cx - half, cy - half), FOV_PX, FOV_PX,
                     fill=False, edgecolor='#ffff00', linewidth=2.5)
    ax_consep_ov.add_patch(rect)

ax_consep_ov.set_xlim(0, slide_dims[0])
ax_consep_ov.set_ylim(slide_dims[1], 0)
ax_consep_ov.axis('off')
ax_consep_ov.set_title('ConSep\n(Epithelial-focused)', fontsize=11, fontweight='bold', pad=5)
ax_consep_ov.text(-0.05, 1.05, 'B', transform=ax_consep_ov.transAxes, fontsize=14, fontweight='bold', va='top')
add_scale_bar(ax_consep_ov, 5000, UM_PER_PIXEL, position='lower right', color='white', fontsize=8)

# Panel C: Hierarchical overview
ax_hier_ov = fig.add_subplot(gs[0, 2])
ax_hier_ov.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=1.0)

hier_counts_ov = {}
for cell_class in hierarchical_colors.keys():
    subset = matched[matched['simplified_class'] == cell_class]
    hier_counts_ov[cell_class] = len(subset)
    if len(subset) > 0:
        ax_hier_ov.scatter(subset['x'], subset['y'], c=hierarchical_colors[cell_class],
                          s=OVERVIEW_POINT_SIZE, alpha=0.7, rasterized=True)

for idx, (cx, cy, _) in enumerate(regions):
    half = FOV_PX // 2
    rect = Rectangle((cx - half, cy - half), FOV_PX, FOV_PX,
                     fill=False, edgecolor='#ffff00', linewidth=2.5)
    ax_hier_ov.add_patch(rect)

ax_hier_ov.set_xlim(0, slide_dims[0])
ax_hier_ov.set_ylim(slide_dims[1], 0)
ax_hier_ov.axis('off')
ax_hier_ov.set_title('Hierarchical\n(Combined)', fontsize=11, fontweight='bold', pad=5)
ax_hier_ov.text(-0.05, 1.05, 'C', transform=ax_hier_ov.transAxes, fontsize=14, fontweight='bold', va='top')
add_scale_bar(ax_hier_ov, 5000, UM_PER_PIXEL, position='lower right', color='white', fontsize=8)

# Panel D: Legend panel - All three classification schemes
ax_legend = fig.add_subplot(gs[0, 3])
ax_legend.axis('off')
ax_legend.set_xlim(0, 1)
ax_legend.set_ylim(0, 1)

# Draw color swatches - very compact to fit all three
y_start = 0.98
y_step = 0.028

# CoNIC legend
ax_legend.text(0.02, y_start, 'CoNIC', fontsize=7, fontweight='bold', transform=ax_legend.transAxes)
y_pos = y_start - 0.018
for ct, color in conic_colors.items():
    ax_legend.add_patch(FancyBboxPatch((0.02, y_pos - 0.008), 0.05, 0.015,
                        boxstyle="round,pad=0.003", facecolor=color, edgecolor='none',
                        transform=ax_legend.transAxes))
    ax_legend.text(0.08, y_pos, ct, fontsize=5, va='center', transform=ax_legend.transAxes)
    y_pos -= y_step

# ConSep legend
y_pos -= 0.008
ax_legend.text(0.02, y_pos, 'ConSep', fontsize=7, fontweight='bold', transform=ax_legend.transAxes)
y_pos -= 0.018
for ct, color in consep_colors.items():
    ax_legend.add_patch(FancyBboxPatch((0.02, y_pos - 0.008), 0.05, 0.015,
                        boxstyle="round,pad=0.003", facecolor=color, edgecolor='none',
                        transform=ax_legend.transAxes))
    display = ct[:14] if len(ct) > 14 else ct
    ax_legend.text(0.08, y_pos, display, fontsize=5, va='center', transform=ax_legend.transAxes)
    y_pos -= y_step

# Hierarchical legend (most important classes only)
y_pos -= 0.008
ax_legend.text(0.02, y_pos, 'Hierarchical', fontsize=7, fontweight='bold', transform=ax_legend.transAxes)
y_pos -= 0.018

# Show main hierarchical classes (sorted by typical importance)
hier_legend_order = [
    'Malignant Epithelial', 'Healthy Epithelial', 'Lymphocyte', 'Plasma cell',
    'Stroma', 'Muscle', 'Eosinophil', 'Neutrophil'
]
for ct in hier_legend_order:
    if ct in hierarchical_colors:
        color = hierarchical_colors[ct]
        ax_legend.add_patch(FancyBboxPatch((0.02, y_pos - 0.008), 0.05, 0.015,
                            boxstyle="round,pad=0.003", facecolor=color, edgecolor='none',
                            transform=ax_legend.transAxes))
        display = ct[:14] if len(ct) > 14 else ct
        ax_legend.text(0.08, y_pos, display, fontsize=5, va='center', transform=ax_legend.transAxes)
        y_pos -= y_step

# Rows 2-3: Regional crops
panel_labels = ['D', 'E', 'F', 'G', 'H', 'I']
label_idx = 0

for row_idx, (cx, cy, region_name) in enumerate(regions):
    half = FOV_PX // 2
    x_min, x_max = cx - half, cx + half
    y_min, y_max = cy - half, cy + half

    # Get H&E crop
    he_region = he_slide.read_region((x_min, y_min), 0, (FOV_PX, FOV_PX))
    he_crop = np.array(he_region.convert('RGB'))

    # CoNIC crop
    ax_conic = fig.add_subplot(gs[row_idx + 1, 0])
    ax_conic.imshow(he_crop, extent=[x_min, x_max, y_max, y_min])
    conic_counts, conic_total = plot_contours(ax_conic, contours_conic, conic_colors, x_min, x_max, y_min, y_max)
    ax_conic.set_xlim(x_min, x_max)
    ax_conic.set_ylim(y_max, y_min)
    ax_conic.axis('off')
    if row_idx == 0:
        ax_conic.set_title('CoNIC', fontsize=9, fontweight='bold', pad=2)
    ax_conic.text(-0.02, 1.02, panel_labels[label_idx], transform=ax_conic.transAxes,
                  fontsize=11, fontweight='bold', va='top')
    ax_conic.text(0.02, 0.98, f'{row_idx+1}. {region_name}', transform=ax_conic.transAxes,
                  fontsize=6, fontweight='bold', va='top', color='white',
                  bbox=dict(boxstyle='round,pad=0.15', facecolor='black', alpha=0.7, edgecolor='none'))
    add_scale_bar(ax_conic, 50, UM_PER_PIXEL, position='lower right', color='white', fontsize=6)
    label_idx += 1

    # ConSep crop
    ax_consep = fig.add_subplot(gs[row_idx + 1, 1])
    ax_consep.imshow(he_crop, extent=[x_min, x_max, y_max, y_min])
    consep_counts, consep_total = plot_contours(ax_consep, contours_consep, consep_colors, x_min, x_max, y_min, y_max)
    ax_consep.set_xlim(x_min, x_max)
    ax_consep.set_ylim(y_max, y_min)
    ax_consep.axis('off')
    if row_idx == 0:
        ax_consep.set_title('ConSep', fontsize=9, fontweight='bold', pad=2)
    ax_consep.text(-0.02, 1.02, panel_labels[label_idx], transform=ax_consep.transAxes,
                   fontsize=11, fontweight='bold', va='top')
    add_scale_bar(ax_consep, 50, UM_PER_PIXEL, position='lower right', color='white', fontsize=6)
    label_idx += 1

    # Hierarchical crop
    ax_hier = fig.add_subplot(gs[row_idx + 1, 2])
    ax_hier.imshow(he_crop, extent=[x_min, x_max, y_max, y_min])
    hier_counts, hier_total = plot_hierarchical_points(ax_hier, matched, x_min, x_max, y_min, y_max)
    ax_hier.set_xlim(x_min, x_max)
    ax_hier.set_ylim(y_max, y_min)
    ax_hier.axis('off')
    if row_idx == 0:
        ax_hier.set_title('Hierarchical', fontsize=9, fontweight='bold', pad=2)
    ax_hier.text(-0.02, 1.02, panel_labels[label_idx], transform=ax_hier.transAxes,
                 fontsize=11, fontweight='bold', va='top')
    add_scale_bar(ax_hier, 50, UM_PER_PIXEL, position='lower right', color='white', fontsize=6)
    label_idx += 1

    # Regional legend in the last column
    ax_reg_legend = fig.add_subplot(gs[row_idx + 1, 3])
    ax_reg_legend.axis('off')
    ax_reg_legend.set_xlim(0, 1)
    ax_reg_legend.set_ylim(0, 1)

    # Show top classes for this region - more compact
    y_pos = 0.92
    ax_reg_legend.text(0.02, y_pos, f'Region {row_idx+1}', fontsize=7, fontweight='bold',
                       transform=ax_reg_legend.transAxes)
    y_pos -= 0.06
    ax_reg_legend.text(0.02, y_pos, f'n={hier_total:,}', fontsize=6,
                       transform=ax_reg_legend.transAxes, style='italic')
    y_pos -= 0.09

    # Show hierarchical class breakdown - compact
    sorted_hier = sorted([(k, v) for k, v in hier_counts.items() if v > 0],
                         key=lambda x: x[1], reverse=True)[:6]
    for ct, count in sorted_hier:
        pct = 100 * count / hier_total if hier_total > 0 else 0
        color = hierarchical_colors.get(ct, '#888')
        ax_reg_legend.add_patch(FancyBboxPatch((0.02, y_pos - 0.015), 0.08, 0.03,
                                boxstyle="round,pad=0.005", facecolor=color, edgecolor='none',
                                transform=ax_reg_legend.transAxes))
        # Truncate long names
        display_name = ct[:10] if len(ct) > 10 else ct
        ax_reg_legend.text(0.12, y_pos, f'{display_name} ({pct:.0f}%)', fontsize=5,
                           va='center', transform=ax_reg_legend.transAxes)
        y_pos -= 0.09

# Main title
fig.suptitle('Hierarchical Cell Classification: Combining CoNIC and ConSep Models',
             fontsize=14, fontweight='bold', y=0.98)

# Save main figure
output_path = f"{BASE_DIR}/hierarchical_classification_figure.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nFigure saved to {output_path}")

pdf_path = f"{BASE_DIR}/hierarchical_classification_figure.pdf"
plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"PDF saved to {pdf_path}")

plt.close()

# ============================================================
# SUPPLEMENTARY FIGURE: Stacked bar comparison
# ============================================================
print("\nCreating stacked bar comparison figure...")

fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5))

# Prepare data for stacked bars
conic_data = pd.Series(conic_counts_ov)
consep_data = pd.Series(consep_counts_ov)
hier_data = pd.Series(hier_counts_ov)

# Panel 1: CoNIC stacked bar
ax_bar1 = axes2[0]
conic_sorted = conic_data.sort_values(ascending=True)
colors1 = [conic_colors[c] for c in conic_sorted.index]
y_pos = 0
bars1 = []
for ct, count in conic_sorted.items():
    bar = ax_bar1.barh(0, count, left=y_pos, color=conic_colors[ct], edgecolor='white', linewidth=0.5)
    bars1.append((ct, count, y_pos))
    y_pos += count

# Add labels inside bars
total_conic = conic_sorted.sum()
for ct, count, left in bars1:
    pct = 100 * count / total_conic
    if pct > 5:  # Only label if > 5%
        ax_bar1.text(left + count/2, 0, f'{ct}\n{pct:.1f}%', ha='center', va='center',
                    fontsize=7, fontweight='bold', color='white')
    elif pct > 2:
        ax_bar1.text(left + count/2, 0, f'{pct:.0f}%', ha='center', va='center',
                    fontsize=6, color='white')

ax_bar1.set_xlim(0, total_conic)
ax_bar1.set_ylim(-0.5, 0.5)
ax_bar1.set_yticks([])
ax_bar1.set_xlabel('Cell Count', fontsize=10)
ax_bar1.set_title(f'CoNIC Classification\n({total_conic:,} cells)', fontsize=11, fontweight='bold')
ax_bar1.spines['top'].set_visible(False)
ax_bar1.spines['right'].set_visible(False)
ax_bar1.spines['left'].set_visible(False)
ax_bar1.text(-0.05, 1.1, 'A', transform=ax_bar1.transAxes, fontsize=14, fontweight='bold', va='top')

# Panel 2: ConSep stacked bar
ax_bar2 = axes2[1]
consep_sorted = consep_data.sort_values(ascending=True)
y_pos = 0
bars2 = []
for ct, count in consep_sorted.items():
    ax_bar2.barh(0, count, left=y_pos, color=consep_colors[ct], edgecolor='white', linewidth=0.5)
    bars2.append((ct, count, y_pos))
    y_pos += count

total_consep = consep_sorted.sum()
for ct, count, left in bars2:
    pct = 100 * count / total_consep
    if pct > 5:
        ax_bar2.text(left + count/2, 0, f'{ct}\n{pct:.1f}%', ha='center', va='center',
                    fontsize=7, fontweight='bold', color='white')
    elif pct > 2:
        ax_bar2.text(left + count/2, 0, f'{pct:.0f}%', ha='center', va='center',
                    fontsize=6, color='white')

ax_bar2.set_xlim(0, total_consep)
ax_bar2.set_ylim(-0.5, 0.5)
ax_bar2.set_yticks([])
ax_bar2.set_xlabel('Cell Count', fontsize=10)
ax_bar2.set_title(f'ConSep Classification\n({total_consep:,} cells)', fontsize=11, fontweight='bold')
ax_bar2.spines['top'].set_visible(False)
ax_bar2.spines['right'].set_visible(False)
ax_bar2.spines['left'].set_visible(False)
ax_bar2.text(-0.05, 1.1, 'B', transform=ax_bar2.transAxes, fontsize=14, fontweight='bold', va='top')

# Panel 3: Hierarchical stacked bar
ax_bar3 = axes2[2]
hier_sorted = hier_data[hier_data > 0].sort_values(ascending=True)
y_pos = 0
bars3 = []
for ct, count in hier_sorted.items():
    ax_bar3.barh(0, count, left=y_pos, color=hierarchical_colors.get(ct, '#888'),
                 edgecolor='white', linewidth=0.5)
    bars3.append((ct, count, y_pos))
    y_pos += count

total_hier = hier_sorted.sum()
for ct, count, left in bars3:
    pct = 100 * count / total_hier
    if pct > 5:
        # Truncate name
        display_name = ct[:10] + '..' if len(ct) > 12 else ct
        ax_bar3.text(left + count/2, 0, f'{display_name}\n{pct:.1f}%', ha='center', va='center',
                    fontsize=6, fontweight='bold', color='white')
    elif pct > 2:
        ax_bar3.text(left + count/2, 0, f'{pct:.0f}%', ha='center', va='center',
                    fontsize=5, color='white')

ax_bar3.set_xlim(0, total_hier)
ax_bar3.set_ylim(-0.5, 0.5)
ax_bar3.set_yticks([])
ax_bar3.set_xlabel('Cell Count', fontsize=10)
ax_bar3.set_title(f'Hierarchical Classification\n({total_hier:,} matched cells)', fontsize=11, fontweight='bold')
ax_bar3.spines['top'].set_visible(False)
ax_bar3.spines['right'].set_visible(False)
ax_bar3.spines['left'].set_visible(False)
ax_bar3.text(-0.05, 1.1, 'C', transform=ax_bar3.transAxes, fontsize=14, fontweight='bold', va='top')

fig2.suptitle('Cell Type Distribution Comparison', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

bar_output = f"{BASE_DIR}/classification_stacked_bars.png"
plt.savefig(bar_output, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Stacked bar figure saved to {bar_output}")

bar_pdf = f"{BASE_DIR}/classification_stacked_bars.pdf"
plt.savefig(bar_pdf, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Stacked bar PDF saved to {bar_pdf}")

plt.show()
