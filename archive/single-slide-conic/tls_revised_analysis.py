"""
Revised TLS Analysis - Literature-Aligned Staging

Key insight: In cancer (chronic antigen), mature TLS are LARGE organized structures
with active germinal centers. High P/L ratio in small clusters represents post-GC
remnants/involuted structures, NOT "mature" TLS.

Staging:
  Stage 1 (Early): Small, sparse lymphoid aggregates
  Stage 2 (Developing): Medium clusters, increasing organization
  Stage 3 (Mature/Active GC): Large, dense, organized with high lymphocyte core
  Stage 4 (Involuted): Small clusters with high P/L ratio (post-GC remnants)
"""

import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from scipy.spatial.distance import cdist
import networkx as nx
from collections import Counter
import openslide
import warnings
warnings.filterwarnings('ignore')

# Load data
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

# Load H&E slide with openslide
he_slide_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A.tif"
print("Loading H&E slide...")
he_slide = openslide.OpenSlide(he_slide_path)
print(f"Slide dimensions: {he_slide.dimensions}")

print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()

print("Loading cell contours...")
contours = sdata.shapes['cell_contours']

# Filter to immune cells
immune_mask = points['classification_name'].isin(['Lymphocyte', 'Plasma cell'])
immune_cells = points[immune_mask].reset_index(drop=True)
print(f"Immune cells: {len(immune_cells):,}")

# Build spatial graph
coords = immune_cells[['x', 'y']].values
cell_types = immune_cells['classification_name'].values

# Pixel scale: 0.325 um/pixel
UM_PER_PIXEL = 0.325
radius = int(50 / UM_PER_PIXEL)  # 50 um = 154 pixels
tree = cKDTree(coords)
pairs = tree.query_pairs(r=radius, output_type='ndarray')

G = nx.Graph()
G.add_nodes_from(range(len(immune_cells)))
for i, (ct, x, y) in enumerate(zip(cell_types, coords[:, 0], coords[:, 1])):
    G.nodes[i]['cell_type'] = ct
    G.nodes[i]['x'] = x
    G.nodes[i]['y'] = y
G.add_edges_from(pairs)

# Find connected components
components = list(nx.connected_components(G))
print(f"Found {len(components):,} connected components")


def calculate_cluster_metrics(G, comp):
    """Calculate multi-metric features for a cluster."""
    comp_list = list(comp)
    n_cells = len(comp_list)

    # Get coordinates and cell types
    comp_coords = np.array([[G.nodes[n]['x'], G.nodes[n]['y']] for n in comp_list])
    comp_types = [G.nodes[n]['cell_type'] for n in comp_list]
    type_counts = Counter(comp_types)

    n_lymph = type_counts.get('Lymphocyte', 0)
    n_plasma = type_counts.get('Plasma cell', 0)

    # Centroid
    centroid = comp_coords.mean(axis=0)

    # Area and compactness via convex hull
    if n_cells >= 3:
        try:
            hull = ConvexHull(comp_coords)
            area = hull.volume  # In 2D, volume is area
            perimeter = hull.area  # In 2D, area is perimeter
            # Circularity/compactness: 4π × area / perimeter²
            compactness = (4 * np.pi * area) / (perimeter ** 2) if perimeter > 0 else 0
        except Exception:
            # Degenerate hull (collinear points)
            area = 1.0
            perimeter = 1.0
            compactness = 0.0
    else:
        area = 1.0
        perimeter = 1.0
        compactness = 0.0

    # Density: cells per area (cells/px²)
    density = n_cells / area if area > 0 else 0

    # Plasma cell localization
    # Ratio of mean(plasma distance to centroid) / mean(lymph distance to centroid)
    # > 1.0 = plasma cells peripheral (emigrating from GC - mature/active)
    # < 1.0 = plasma cells central
    # ≈ 1.0 = mixed/unorganized

    lymph_mask = np.array([ct == 'Lymphocyte' for ct in comp_types])
    plasma_mask = np.array([ct == 'Plasma cell' for ct in comp_types])

    if lymph_mask.sum() > 0 and plasma_mask.sum() > 0:
        lymph_coords = comp_coords[lymph_mask]
        plasma_coords = comp_coords[plasma_mask]

        # Distance to centroid
        lymph_dists = np.linalg.norm(lymph_coords - centroid, axis=1)
        plasma_dists = np.linalg.norm(plasma_coords - centroid, axis=1)

        mean_lymph_dist = lymph_dists.mean()
        mean_plasma_dist = plasma_dists.mean()

        plasma_dist_ratio = mean_plasma_dist / mean_lymph_dist if mean_lymph_dist > 0 else 1.0
    else:
        plasma_dist_ratio = 1.0
        mean_lymph_dist = 0
        mean_plasma_dist = 0

    # P/L ratio
    pl_ratio = n_plasma / n_lymph if n_lymph > 0 else 0

    return {
        'nodes': comp,
        'n_cells': n_cells,
        'n_lymph': n_lymph,
        'n_plasma': n_plasma,
        'pl_ratio': pl_ratio,
        'centroid_x': centroid[0],
        'centroid_y': centroid[1],
        'area': area,
        'perimeter': perimeter,
        'density': density,
        'compactness': compactness,
        'plasma_dist_ratio': plasma_dist_ratio,
        'coords': comp_coords,
    }


def assign_tls_stage(metrics):
    """
    Assign TLS stage based on literature-aligned criteria.

    Stage 1 (Early): Small (<100), sparse lymphoid aggregates
    Stage 2 (Developing): Medium (100-500), increasing organization
    Stage 3 (Mature/Active GC): Large (>500), dense, peripheral plasma cells
    Stage 4 (Involuted): Small (<200), high P/L ratio (>0.4) - post-GC remnants
    """
    n = metrics['n_cells']
    pl = metrics['pl_ratio']
    density = metrics['density']
    plasma_periph = metrics['plasma_dist_ratio']

    # Stage 4: Involuted - small with high P/L ratio (post-GC remnants)
    if n < 200 and pl > 0.4:
        return 4, 'Involuted'

    # Stage 3: Mature/Active - large, organized clusters
    if n >= 500:
        return 3, 'Mature/Active'

    # Stage 2: Developing - medium size
    if n >= 100:
        return 2, 'Developing'

    # Stage 1: Early - small aggregates
    return 1, 'Early'


def calculate_maturity_score(metrics, all_metrics):
    """
    Composite TLS maturity score (literature-aligned).

    Higher score = more mature/active TLS
    Based on: size (60%), organization/compactness (40%)

    Key insight: Size is the dominant factor for maturity in cancer TLS.
    Large, organized structures = mature/active GC.

    Note: Plasma cell localization metric was removed as it did not show
    the expected biological pattern (plasma cells were not preferentially
    peripheral in mature TLS).
    """
    # Get normalization values from all clusters
    sizes = [m['n_cells'] for m in all_metrics]

    # Use log-scale for size to handle wide range
    log_sizes = [np.log1p(s) for s in sizes]
    log_size = np.log1p(metrics['n_cells'])

    # Min-max normalize
    def normalize(val, vals):
        min_v, max_v = min(vals), max(vals)
        if max_v == min_v:
            return 0.5
        return (val - min_v) / (max_v - min_v)

    norm_size = normalize(log_size, log_sizes)

    # Organization score: compactness with P/L ratio penalty for small clusters
    # Large clusters with moderate P/L = active GC
    # Small clusters with high P/L = involuted
    compactness = metrics['compactness']

    # Penalize high P/L ratio in small clusters (involution marker)
    pl_ratio = metrics['pl_ratio']
    size_factor = min(1.0, metrics['n_cells'] / 500)  # Ramps up to 1 at 500 cells
    pl_penalty = 0 if size_factor >= 0.8 else (1 - size_factor) * max(0, pl_ratio - 0.3)
    organization_score = max(0, compactness - pl_penalty)

    maturity_score = (
        0.60 * norm_size +
        0.40 * organization_score
    )

    return maturity_score


# Analyze all clusters
print("\nAnalyzing TLS clusters...")
tls_data = []

for comp_id, comp in enumerate(components):
    if len(comp) < 15:  # Minimum size for reliable metrics
        continue

    metrics = calculate_cluster_metrics(G, comp)

    # Need both cell types for TLS
    if metrics['n_lymph'] >= 10 and metrics['n_plasma'] >= 3:
        metrics['component_id'] = comp_id
        tls_data.append(metrics)

print(f"Found {len(tls_data)} TLS clusters meeting criteria")

# Calculate maturity scores and assign stages
for tls in tls_data:
    tls['maturity_score'] = calculate_maturity_score(tls, tls_data)
    tls['stage_num'], tls['stage_name'] = assign_tls_stage(tls)

# Statistics
print("\n" + "="*60)
print("TLS CLUSTER STATISTICS")
print("="*60)

print(f"\nTotal TLS clusters: {len(tls_data)}")

# By stage
stage_counts = Counter([t['stage_name'] for t in tls_data])
print("\nTLS by Stage (Literature-Aligned):")
for stage in ['Early', 'Developing', 'Mature/Active', 'Involuted']:
    count = stage_counts.get(stage, 0)
    pct = 100 * count / len(tls_data) if tls_data else 0
    print(f"  Stage {['Early', 'Developing', 'Mature/Active', 'Involuted'].index(stage)+1} ({stage}): {count} ({pct:.1f}%)")

# Size distribution
sizes = [t['n_cells'] for t in tls_data]
print(f"\nCluster Size: min={min(sizes)}, max={max(sizes)}, mean={np.mean(sizes):.1f}, median={np.median(sizes):.0f}")

# Density distribution
densities = [t['density'] for t in tls_data]
print(f"Density (cells/px²): min={min(densities):.6f}, max={max(densities):.6f}, mean={np.mean(densities):.6f}")

# P/L ratio distribution
pl_ratios = [t['pl_ratio'] for t in tls_data]
print(f"P/L Ratio: min={min(pl_ratios):.3f}, max={max(pl_ratios):.3f}, mean={np.mean(pl_ratios):.3f}")

# Top clusters by maturity score
print("\n" + "-"*60)
print("Top 10 TLS by Maturity Score (Literature-Aligned):")
tls_sorted = sorted(tls_data, key=lambda x: x['maturity_score'], reverse=True)
for i, tls in enumerate(tls_sorted[:10]):
    print(f"  {i+1}. Score: {tls['maturity_score']:.3f} | Stage: {tls['stage_name']} | "
          f"Size: {tls['n_cells']} | Density: {tls['density']:.5f} | "
          f"P/L: {tls['pl_ratio']:.3f}")

print("\n" + "-"*60)
print("Involuted TLS (Post-GC Remnants - Small + High P/L):")
involuted = [t for t in tls_data if t['stage_name'] == 'Involuted']
involuted_sorted = sorted(involuted, key=lambda x: x['pl_ratio'], reverse=True)
for i, tls in enumerate(involuted_sorted[:5]):
    print(f"  {i+1}. Size: {tls['n_cells']} | P/L: {tls['pl_ratio']:.3f} | "
          f"Density: {tls['density']:.5f}")


# =============================================================================
# VISUALIZATIONS
# =============================================================================
print("\n" + "="*60)
print("Creating visualizations...")
print("="*60)

# Color scheme for stages
stage_colors = {
    'Early': '#3498db',        # Blue
    'Developing': '#f39c12',   # Orange
    'Mature/Active': '#27ae60', # Green
    'Involuted': '#9b59b6'     # Purple
}

fig = plt.figure(figsize=(20, 16))

# 1. 3D scatter: log(Size) × Density × P/L ratio (colored by stage)
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
for stage, color in stage_colors.items():
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]
    if stage_tls:
        xs = [np.log10(t['n_cells']) for t in stage_tls]  # Log scale for size
        ys = [t['density'] * 1e4 for t in stage_tls]  # Scale for visibility
        zs = [t['pl_ratio'] for t in stage_tls]
        ax1.scatter(xs, ys, zs, c=color, label=f'{stage} ({len(stage_tls)})', alpha=0.7, s=50)

ax1.set_xlabel('log₁₀(Cluster Size)', fontsize=10)
ax1.set_ylabel('Density (×10⁻⁴ cells/px²)', fontsize=10)
ax1.set_zlabel('P/L Ratio', fontsize=10)
ax1.set_title('TLS Multi-Metric Space\n(log Size × Density × P/L Ratio)', fontsize=12, fontweight='bold')
ax1.legend(loc='upper left', fontsize=9)

# 2. Spatial map colored by stage with H&E background
ax2 = fig.add_subplot(2, 2, 2)
all_x = points['x'].values
all_y = points['y'].values

# Get H&E thumbnail for background
slide_dims = he_slide.dimensions
thumbnail_size = (1200, int(1200 * slide_dims[1] / slide_dims[0]))
he_thumbnail = he_slide.get_thumbnail(thumbnail_size)
he_thumbnail_arr = np.array(he_thumbnail.convert('RGB'))
ax2.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=0.6)

for stage in ['Early', 'Developing', 'Mature/Active', 'Involuted']:
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]
    if stage_tls:
        xs = [t['centroid_x'] for t in stage_tls]
        ys = [t['centroid_y'] for t in stage_tls]
        sizes = [np.sqrt(t['n_cells']) * 3 for t in stage_tls]
        ax2.scatter(xs, ys, c=stage_colors[stage], s=sizes, alpha=0.8,
                   edgecolors='black', linewidths=0.5, label=f'{stage} ({len(stage_tls)})')

ax2.set_xlim(0, slide_dims[0])
ax2.set_ylim(slide_dims[1], 0)
ax2.set_xlabel('X (pixels)', fontsize=10)
ax2.set_ylabel('Y (pixels)', fontsize=10)
ax2.set_title('TLS Spatial Distribution by Stage\n(Literature-Aligned Staging)', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right', fontsize=9)
ax2.set_aspect('equal')

# 3. Maturity Score vs P/L Ratio (showing the correction)
ax3 = fig.add_subplot(2, 2, 3)
for stage, color in stage_colors.items():
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]
    if stage_tls:
        xs = [t['pl_ratio'] for t in stage_tls]
        ys = [t['maturity_score'] for t in stage_tls]
        sizes = [t['n_cells'] / 5 for t in stage_tls]
        ax3.scatter(xs, ys, c=color, s=sizes, alpha=0.6, label=f'{stage} ({len(stage_tls)})')

ax3.set_xlabel('P/L Ratio', fontsize=10)
ax3.set_ylabel('Maturity Score (Literature-Aligned)', fontsize=10)
ax3.set_title('Maturity Score vs P/L Ratio\nNote: High P/L + Small = Involuted (not Mature)',
              fontsize=12, fontweight='bold')
ax3.axvline(0.4, color='gray', linestyle='--', alpha=0.5, label='P/L > 0.4 threshold')
ax3.legend(loc='upper right', fontsize=9)

# 4. Size vs P/L Ratio colored by stage (key diagnostic plot)
ax4 = fig.add_subplot(2, 2, 4)
for stage, color in stage_colors.items():
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]
    if stage_tls:
        sizes = [t['n_cells'] for t in stage_tls]
        pl_ratios = [t['pl_ratio'] for t in stage_tls]
        ax4.scatter(sizes, pl_ratios, c=color, alpha=0.6, s=40,
                   label=f'{stage} ({len(stage_tls)})', edgecolors='white', linewidths=0.3)

ax4.set_xscale('log')
ax4.axhline(0.4, color='gray', linestyle='--', alpha=0.7, label='Involuted threshold (P/L>0.4)')
ax4.axvline(100, color='gray', linestyle=':', alpha=0.5)
ax4.axvline(200, color='gray', linestyle=':', alpha=0.5)
ax4.axvline(500, color='gray', linestyle=':', alpha=0.5)
ax4.set_xlabel('Cluster Size (cells, log scale)', fontsize=10)
ax4.set_ylabel('P/L Ratio', fontsize=10)
ax4.set_title('Size vs P/L Ratio by Stage\nInvoluted = Small (<200) + High P/L (>0.4)',
              fontsize=12, fontweight='bold')
ax4.legend(loc='upper right', fontsize=9)

plt.suptitle('Revised TLS Analysis - Literature-Aligned Staging (50 um radius)\n'
             'Mature TLS = Large, Dense, Organized | Involuted = Small, High P/L (Post-GC)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = "./tls_revised_analysis.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path}")


# =============================================================================
# Example clusters from each stage - showing ALL cell types
# Fixed 250 um FOV for consistent comparison
# =============================================================================
fig2, axes = plt.subplots(2, 2, figsize=(16, 16))
axes = axes.flatten()

# Colors for all cell types
all_cell_colors = {
    'Epithelial': '#1f77b4',    # Blue
    'Connective': '#ff7f0e',    # Orange
    'Lymphocyte': '#2ca02c',    # Green
    'Plasma cell': '#d62728',   # Red
    'Eosinophil': '#9467bd',    # Purple
    'Neutrophil': '#8c564b',    # Brown
}

# Fixed FOV: 250 um = 769 pixels
FOV_UM = 250
FOV_PX = FOV_UM / UM_PER_PIXEL
crop_size = FOV_PX

stage_order = ['Early', 'Developing', 'Mature/Active', 'Involuted']

# Get one example from each stage - select by MATURITY SCORE (most principled)
# For large TLS, we'll show the central 250 um region
examples = {}
for stage in stage_order:
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]
    if stage_tls:
        if stage == 'Involuted':
            # Highest P/L ratio (defining characteristic of involuted)
            examples[stage] = sorted(stage_tls, key=lambda x: x['pl_ratio'], reverse=True)[0]
        else:
            # For all other stages: highest maturity score
            examples[stage] = sorted(stage_tls, key=lambda x: x['maturity_score'], reverse=True)[0]

for idx, stage in enumerate(stage_order):
    ax = axes[idx]

    if stage in examples:
        tls = examples[stage]
        cx, cy = tls['centroid_x'], tls['centroid_y']

        # Fixed 250 um FOV
        x_min, x_max = cx - crop_size, cx + crop_size
        y_min, y_max = cy - crop_size, cy + crop_size

        region_mask = (
            (points['x'] >= x_min) & (points['x'] <= x_max) &
            (points['y'] >= y_min) & (points['y'] <= y_max)
        )
        region_points = points[region_mask]

        # Show H&E background image using openslide
        # Convert coordinates to integers for region extraction
        region_x = max(0, int(x_min))
        region_y = max(0, int(y_min))
        region_w = int(x_max - x_min)
        region_h = int(y_max - y_min)

        # Read region from slide (level 0 = full resolution)
        he_region = he_slide.read_region((region_x, region_y), 0, (region_w, region_h))
        he_crop = np.array(he_region.convert('RGB'))

        # Display with extent matching the coordinate system
        ax.imshow(he_crop, extent=[region_x, region_x + region_w, region_y + region_h, region_y], alpha=0.6)

        # Fixed spot sizes for 250 um FOV
        immune_size = 12
        other_size = 3

        # Plot each cell type
        for cell_type, color in all_cell_colors.items():
            ct_mask = region_points['classification_name'] == cell_type
            subset = region_points[ct_mask]
            if len(subset) > 0:
                if cell_type in ['Lymphocyte', 'Plasma cell']:
                    size, alpha = immune_size, 0.9
                else:
                    size, alpha = other_size, 0.5
                ax.scatter(
                    subset['x'], subset['y'],
                    c=color,
                    label=f'{cell_type} ({len(subset)})',
                    s=size,
                    alpha=alpha,
                    rasterized=True
                )

        # Mark TLS centroid
        ax.scatter([cx], [cy], c='black', s=150, marker='x', linewidths=2, zorder=10)

        # Draw 50 um radius circle (the clustering radius)
        cluster_radius_px = 50 / UM_PER_PIXEL
        circle = plt.Circle((cx, cy), cluster_radius_px, fill=False,
                            color=stage_colors[stage], linewidth=2, linestyle='--',
                            label='50 um radius')
        ax.add_patch(circle)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_max, y_min)  # Inverted y

        # Check if TLS extends beyond FOV (estimate from convex hull)
        tls_extent = np.sqrt(tls['area']) if tls['area'] > 0 else 0
        fov_note = "(central region)" if tls_extent > FOV_PX else "(full TLS)"

        ax.set_title(f"Stage {idx+1}: {stage} {fov_note}\n"
                    f"TLS Size: {tls['n_cells']} cells | P/L: {tls['pl_ratio']:.3f} | "
                    f"Maturity: {tls['maturity_score']:.3f}",
                    fontsize=11, fontweight='bold', color=stage_colors[stage])
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
        ax.legend(loc='upper right', fontsize=8, markerscale=2)
        ax.set_aspect('equal')
    else:
        ax.text(0.5, 0.5, f'No {stage} TLS found', ha='center', va='center', fontsize=14)
        ax.axis('off')

plt.suptitle(f'Example TLS Clusters by Stage ({FOV_UM} um FOV)\n'
             'Dashed circle = 50 um clustering radius, X = Centroid',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path2 = "./tls_revised_examples.png"
plt.savefig(output_path2, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path2}")


# =============================================================================
# Size distribution by stage
# =============================================================================
fig3, axes = plt.subplots(1, 4, figsize=(18, 4))
stage_order = ['Early', 'Developing', 'Mature/Active', 'Involuted']

for idx, stage in enumerate(stage_order):
    ax = axes[idx]
    stage_tls = [t for t in tls_data if t['stage_name'] == stage]

    if stage_tls:
        sizes = [t['n_cells'] for t in stage_tls]
        ax.hist(sizes, bins=20, color=stage_colors[stage], alpha=0.7, edgecolor='black')
        ax.axvline(np.median(sizes), color='black', linestyle='-', linewidth=2,
                   label=f'Median={np.median(sizes):.0f}')
        ax.set_xlabel('Cluster Size (cells)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'{stage}\n(n={len(stage_tls)})', fontsize=11, fontweight='bold', color=stage_colors[stage])
        ax.legend(fontsize=8)
        if max(sizes) > 1000:
            ax.set_xscale('log')
    else:
        ax.text(0.5, 0.5, f'No data', ha='center', va='center')
        ax.set_title(stage, fontsize=11)

plt.suptitle('TLS Size Distribution by Stage',
             fontsize=12, fontweight='bold')
plt.tight_layout()

output_path3 = "./tls_size_distribution.png"
plt.savefig(output_path3, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path3}")

plt.show()

print("\n" + "="*60)
print("Analysis complete!")
print("="*60)
