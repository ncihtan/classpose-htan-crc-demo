"""
TLS Radius Analysis - Comparing 50, 100, 200 um clustering radii

Pixel scale: 0.325 um/pixel
- 50 um  = 154 pixels
- 100 um = 308 pixels
- 200 um = 615 pixels
"""

import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
import networkx as nx
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Pixel scale
UM_PER_PIXEL = 0.325

# Radii to test (in microns)
RADII_UM = [50, 100, 200]
RADII_PX = [int(r / UM_PER_PIXEL) for r in RADII_UM]

print(f"Pixel scale: {UM_PER_PIXEL} um/pixel")
print(f"Testing radii: {RADII_UM} um = {RADII_PX} pixels")

# Load data
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

print("\nLoading cell centroids...")
points = sdata.points['cell_centroids'].compute()

# Filter to immune cells
immune_mask = points['classification_name'].isin(['Lymphocyte', 'Plasma cell'])
immune_cells = points[immune_mask].reset_index(drop=True)
print(f"Immune cells: {len(immune_cells):,}")

coords = immune_cells[['x', 'y']].values
cell_types = immune_cells['classification_name'].values


def analyze_tls_at_radius(coords, cell_types, radius_px, radius_um):
    """Run TLS analysis at a specific radius."""

    # Build spatial graph
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=radius_px, output_type='ndarray')

    G = nx.Graph()
    G.add_nodes_from(range(len(coords)))
    for i, (ct, x, y) in enumerate(zip(cell_types, coords[:, 0], coords[:, 1])):
        G.nodes[i]['cell_type'] = ct
        G.nodes[i]['x'] = x
        G.nodes[i]['y'] = y
    G.add_edges_from(pairs)

    # Find connected components
    components = list(nx.connected_components(G))

    # Analyze TLS clusters
    tls_data = []
    for comp_id, comp in enumerate(components):
        if len(comp) < 15:
            continue

        comp_list = list(comp)
        comp_coords = np.array([[G.nodes[n]['x'], G.nodes[n]['y']] for n in comp_list])
        comp_types = [G.nodes[n]['cell_type'] for n in comp_list]
        type_counts = Counter(comp_types)

        n_cells = len(comp_list)
        n_lymph = type_counts.get('Lymphocyte', 0)
        n_plasma = type_counts.get('Plasma cell', 0)

        if n_lymph >= 10 and n_plasma >= 3:
            centroid = comp_coords.mean(axis=0)

            # Area via convex hull
            if n_cells >= 3:
                try:
                    hull = ConvexHull(comp_coords)
                    area = hull.volume
                except:
                    area = 1.0
            else:
                area = 1.0

            # Convert area to um²
            area_um2 = area * (UM_PER_PIXEL ** 2)

            # Density
            density = n_cells / area_um2 if area_um2 > 0 else 0

            # P/L ratio
            pl_ratio = n_plasma / n_lymph if n_lymph > 0 else 0

            # Plasma localization
            lymph_mask = np.array([ct == 'Lymphocyte' for ct in comp_types])
            plasma_mask = np.array([ct == 'Plasma cell' for ct in comp_types])

            if lymph_mask.sum() > 0 and plasma_mask.sum() > 0:
                lymph_coords = comp_coords[lymph_mask]
                plasma_coords = comp_coords[plasma_mask]
                lymph_dists = np.linalg.norm(lymph_coords - centroid, axis=1)
                plasma_dists = np.linalg.norm(plasma_coords - centroid, axis=1)
                plasma_dist_ratio = plasma_dists.mean() / lymph_dists.mean() if lymph_dists.mean() > 0 else 1.0
            else:
                plasma_dist_ratio = 1.0

            # Staging
            if n_cells < 200 and pl_ratio > 0.4:
                stage = 'Involuted'
            elif n_cells >= 500:
                stage = 'Mature/Active'
            elif n_cells >= 100:
                stage = 'Developing'
            else:
                stage = 'Early'

            tls_data.append({
                'n_cells': n_cells,
                'n_lymph': n_lymph,
                'n_plasma': n_plasma,
                'pl_ratio': pl_ratio,
                'area_um2': area_um2,
                'density': density,
                'plasma_dist_ratio': plasma_dist_ratio,
                'centroid_x': centroid[0],
                'centroid_y': centroid[1],
                'stage': stage,
                'nodes': comp
            })

    return {
        'radius_um': radius_um,
        'radius_px': radius_px,
        'n_components': len(components),
        'n_tls': len(tls_data),
        'tls_data': tls_data
    }


# Run analysis at each radius
print("\n" + "="*70)
results = {}
for radius_um, radius_px in zip(RADII_UM, RADII_PX):
    print(f"\nAnalyzing at {radius_um} um ({radius_px} px)...")
    results[radius_um] = analyze_tls_at_radius(coords, cell_types, radius_px, radius_um)

    r = results[radius_um]
    print(f"  Components: {r['n_components']:,}")
    print(f"  TLS clusters: {r['n_tls']}")


# Summary comparison
print("\n" + "="*70)
print("RADIUS COMPARISON SUMMARY")
print("="*70)

print(f"\n{'Metric':<30} {'50 um':>12} {'100 um':>12} {'200 um':>12}")
print("-"*70)

# Basic counts
print(f"{'Total components':<30} {results[50]['n_components']:>12,} {results[100]['n_components']:>12,} {results[200]['n_components']:>12,}")
print(f"{'TLS clusters (≥15 cells)':<30} {results[50]['n_tls']:>12} {results[100]['n_tls']:>12} {results[200]['n_tls']:>12}")

# Stage counts
for stage in ['Early', 'Developing', 'Mature/Active', 'Involuted']:
    counts = []
    for r in RADII_UM:
        count = sum(1 for t in results[r]['tls_data'] if t['stage'] == stage)
        counts.append(count)
    print(f"{'  ' + stage:<30} {counts[0]:>12} {counts[1]:>12} {counts[2]:>12}")

# Size statistics
print()
for stat_name, stat_func in [('Mean size (cells)', np.mean), ('Median size', np.median), ('Max size', np.max)]:
    values = []
    for r in RADII_UM:
        if results[r]['tls_data']:
            val = stat_func([t['n_cells'] for t in results[r]['tls_data']])
            values.append(f"{val:.0f}")
        else:
            values.append("N/A")
    print(f"{stat_name:<30} {values[0]:>12} {values[1]:>12} {values[2]:>12}")

# P/L ratio
print()
for stat_name, stat_func in [('Mean P/L ratio', np.mean), ('Median P/L ratio', np.median)]:
    values = []
    for r in RADII_UM:
        if results[r]['tls_data']:
            val = stat_func([t['pl_ratio'] for t in results[r]['tls_data']])
            values.append(f"{val:.3f}")
        else:
            values.append("N/A")
    print(f"{stat_name:<30} {values[0]:>12} {values[1]:>12} {values[2]:>12}")

# Area statistics
print()
for stat_name, stat_func in [('Mean area (um²)', np.mean), ('Median area (um²)', np.median)]:
    values = []
    for r in RADII_UM:
        if results[r]['tls_data']:
            val = stat_func([t['area_um2'] for t in results[r]['tls_data']])
            values.append(f"{val:.0f}")
        else:
            values.append("N/A")
    print(f"{stat_name:<30} {values[0]:>12} {values[1]:>12} {values[2]:>12}")


# =============================================================================
# VISUALIZATIONS
# =============================================================================
print("\n" + "="*70)
print("Creating visualizations...")

stage_colors = {
    'Early': '#3498db',
    'Developing': '#f39c12',
    'Mature/Active': '#27ae60',
    'Involuted': '#9b59b6'
}

fig = plt.figure(figsize=(20, 16))

# 1. Bar chart: TLS counts by stage for each radius
ax1 = fig.add_subplot(2, 2, 1)
x = np.arange(len(['Early', 'Developing', 'Mature/Active', 'Involuted']))
width = 0.25
stages = ['Early', 'Developing', 'Mature/Active', 'Involuted']

for i, (radius, offset) in enumerate(zip(RADII_UM, [-width, 0, width])):
    counts = [sum(1 for t in results[radius]['tls_data'] if t['stage'] == s) for s in stages]
    bars = ax1.bar(x + offset, counts, width, label=f'{radius} um', alpha=0.8)

ax1.set_xlabel('TLS Stage', fontsize=11)
ax1.set_ylabel('Number of TLS', fontsize=11)
ax1.set_title('TLS Stage Distribution by Clustering Radius', fontsize=12, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(stages)
ax1.legend(title='Radius')
ax1.grid(axis='y', alpha=0.3)

# 2. Size distributions
ax2 = fig.add_subplot(2, 2, 2)
for radius, color in zip(RADII_UM, ['#e74c3c', '#3498db', '#2ecc71']):
    sizes = [t['n_cells'] for t in results[radius]['tls_data']]
    if sizes:
        ax2.hist(sizes, bins=50, alpha=0.5, label=f'{radius} um (n={len(sizes)})', color=color, density=True)

ax2.set_xlabel('Cluster Size (cells)', fontsize=11)
ax2.set_ylabel('Density', fontsize=11)
ax2.set_title('TLS Size Distribution by Radius', fontsize=12, fontweight='bold')
ax2.set_xscale('log')
ax2.legend()
ax2.grid(alpha=0.3)

# 3. Spatial maps side by side
ax3 = fig.add_subplot(2, 2, 3)
all_x = points['x'].values
all_y = points['y'].values

# Plot background
ax3.scatter(all_x, all_y, c='lightgray', s=0.01, alpha=0.1, rasterized=True)

# Plot TLS for each radius with different markers
markers = ['o', 's', '^']
for radius, marker in zip(RADII_UM, markers):
    for stage, color in stage_colors.items():
        stage_tls = [t for t in results[radius]['tls_data'] if t['stage'] == stage]
        if stage_tls:
            xs = [t['centroid_x'] for t in stage_tls]
            ys = [t['centroid_y'] for t in stage_tls]
            sizes = [np.sqrt(t['n_cells']) * 2 for t in stage_tls]
            ax3.scatter(xs, ys, c=color, s=sizes, alpha=0.4, marker=marker,
                       edgecolors='black', linewidths=0.2)

ax3.set_xlim(all_x.min(), all_x.max())
ax3.set_ylim(all_y.max(), all_y.min())
ax3.set_xlabel('X (pixels)', fontsize=10)
ax3.set_ylabel('Y (pixels)', fontsize=10)
ax3.set_title('TLS Spatial Distribution\n(○=50um, □=100um, △=200um)', fontsize=12, fontweight='bold')
ax3.set_aspect('equal')

# Add legend for stages
for stage, color in stage_colors.items():
    ax3.scatter([], [], c=color, s=100, label=stage)
ax3.legend(loc='upper right', fontsize=9)

# 4. P/L ratio vs Size for each radius
ax4 = fig.add_subplot(2, 2, 4)
for radius, color, marker in zip(RADII_UM, ['#e74c3c', '#3498db', '#2ecc71'], ['o', 's', '^']):
    tls = results[radius]['tls_data']
    if tls:
        sizes = [t['n_cells'] for t in tls]
        pl_ratios = [t['pl_ratio'] for t in tls]
        ax4.scatter(sizes, pl_ratios, c=color, alpha=0.5, s=30, marker=marker, label=f'{radius} um')

ax4.set_xlabel('Cluster Size (cells)', fontsize=11)
ax4.set_ylabel('P/L Ratio', fontsize=11)
ax4.set_title('P/L Ratio vs Size by Radius', fontsize=12, fontweight='bold')
ax4.set_xscale('log')
ax4.axhline(0.4, color='gray', linestyle='--', alpha=0.5, label='Involuted threshold')
ax4.axvline(100, color='gray', linestyle=':', alpha=0.5)
ax4.axvline(500, color='gray', linestyle=':', alpha=0.5)
ax4.legend()
ax4.grid(alpha=0.3)

plt.suptitle(f'TLS Radius Analysis: 50 vs 100 vs 200 um\nPixel scale: {UM_PER_PIXEL} um/pixel',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = "./tls_radius_comparison.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path}")


# =============================================================================
# Detailed comparison of largest TLS at each radius
# =============================================================================
fig2, axes = plt.subplots(1, 3, figsize=(18, 6))

all_cell_colors = {
    'Epithelial': '#1f77b4',
    'Connective': '#ff7f0e',
    'Lymphocyte': '#2ca02c',
    'Plasma cell': '#d62728',
    'Eosinophil': '#9467bd',
    'Neutrophil': '#8c564b',
}

for idx, (radius, ax) in enumerate(zip(RADII_UM, axes)):
    tls_list = results[radius]['tls_data']
    if tls_list:
        # Get largest TLS
        largest = max(tls_list, key=lambda x: x['n_cells'])
        cx, cy = largest['centroid_x'], largest['centroid_y']

        # Crop size based on area
        crop_size = max(500, np.sqrt(largest['area_um2'] / (UM_PER_PIXEL**2)) + 300)

        x_min, x_max = cx - crop_size, cx + crop_size
        y_min, y_max = cy - crop_size, cy + crop_size

        region_mask = (
            (points['x'] >= x_min) & (points['x'] <= x_max) &
            (points['y'] >= y_min) & (points['y'] <= y_max)
        )
        region_points = points[region_mask]

        # Plot all cell types
        for cell_type, color in all_cell_colors.items():
            ct_mask = region_points['classification_name'] == cell_type
            subset = region_points[ct_mask]
            if len(subset) > 0:
                size = 15 if cell_type in ['Lymphocyte', 'Plasma cell'] else 4
                alpha = 0.9 if cell_type in ['Lymphocyte', 'Plasma cell'] else 0.5
                ax.scatter(subset['x'], subset['y'], c=color, s=size, alpha=alpha,
                          label=f'{cell_type} ({len(subset)})', rasterized=True)

        # Mark centroid
        ax.scatter([cx], [cy], c='black', s=200, marker='x', linewidths=3, zorder=10)

        # Draw radius circle
        circle = plt.Circle((cx, cy), radius / UM_PER_PIXEL, fill=False,
                            color='black', linewidth=2, linestyle='--')
        ax.add_patch(circle)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_max, y_min)
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
        ax.set_title(f'Largest TLS at {radius} um radius\n'
                    f'Size: {largest["n_cells"]} cells | P/L: {largest["pl_ratio"]:.3f}\n'
                    f'Stage: {largest["stage"]}',
                    fontsize=11, fontweight='bold')
        ax.legend(loc='upper right', fontsize=7, markerscale=1.5)
        ax.set_aspect('equal')

plt.suptitle('Largest TLS at Each Clustering Radius\nDashed circle = clustering radius',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path2 = "./tls_radius_largest.png"
plt.savefig(output_path2, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path2}")


# =============================================================================
# Metrics table
# =============================================================================
print("\n" + "="*70)
print("DETAILED METRICS BY RADIUS")
print("="*70)

for radius in RADII_UM:
    r = results[radius]
    tls = r['tls_data']

    print(f"\n{'='*30}")
    print(f"RADIUS: {radius} um ({r['radius_px']} pixels)")
    print(f"{'='*30}")

    if tls:
        sizes = [t['n_cells'] for t in tls]
        pl_ratios = [t['pl_ratio'] for t in tls]
        areas = [t['area_um2'] for t in tls]
        densities = [t['density'] for t in tls]

        print(f"TLS count: {len(tls)}")
        print(f"\nSize (cells):")
        print(f"  Range: {min(sizes)} - {max(sizes)}")
        print(f"  Mean: {np.mean(sizes):.1f}, Median: {np.median(sizes):.0f}")
        print(f"  25th/75th percentile: {np.percentile(sizes, 25):.0f} / {np.percentile(sizes, 75):.0f}")

        print(f"\nArea (um²):")
        print(f"  Range: {min(areas):.0f} - {max(areas):.0f}")
        print(f"  Mean: {np.mean(areas):.0f}, Median: {np.median(areas):.0f}")

        print(f"\nP/L Ratio:")
        print(f"  Range: {min(pl_ratios):.3f} - {max(pl_ratios):.3f}")
        print(f"  Mean: {np.mean(pl_ratios):.3f}, Median: {np.median(pl_ratios):.3f}")

        print(f"\nStage distribution:")
        stage_counts = Counter([t['stage'] for t in tls])
        for stage in ['Early', 'Developing', 'Mature/Active', 'Involuted']:
            count = stage_counts.get(stage, 0)
            pct = 100 * count / len(tls)
            print(f"  {stage}: {count} ({pct:.1f}%)")
    else:
        print("No TLS found")

plt.show()

print("\n" + "="*70)
print("Radius analysis complete!")
print("="*70)
