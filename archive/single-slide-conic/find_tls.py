import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
import pandas as pd

# Load the SpatialData object
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

# Get points with cell type info
print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()
print(f"Total cells: {len(points):,}")

# Extract lymphocytes and plasma cells
lymphocytes = points[points['classification_name'] == 'Lymphocyte'][['x', 'y']].values
plasma_cells = points[points['classification_name'] == 'Plasma cell'][['x', 'y']].values
print(f"Lymphocytes: {len(lymphocytes):,}")
print(f"Plasma cells: {len(plasma_cells):,}")

# TLS detection parameters
neighbor_radius = 100  # pixels - radius to count neighbors
min_lymphocytes = 20   # minimum lymphocytes in radius to consider
min_plasma = 5         # minimum plasma cells in radius to consider

print(f"\nSearching for TLS candidates (radius={neighbor_radius}px, min_lymph={min_lymphocytes}, min_plasma={min_plasma})...")

# Build KD-trees for fast neighbor queries
lymph_tree = cKDTree(lymphocytes)
plasma_tree = cKDTree(plasma_cells)

# For each plasma cell, count nearby lymphocytes and plasma cells
# TLS should have both cell types in proximity
plasma_lymph_counts = lymph_tree.query_ball_point(plasma_cells, neighbor_radius, return_length=True)
plasma_plasma_counts = plasma_tree.query_ball_point(plasma_cells, neighbor_radius, return_length=True)

# Find plasma cells that are in TLS-like environments
tls_mask = (plasma_lymph_counts >= min_lymphocytes) & (plasma_plasma_counts >= min_plasma)
tls_plasma_cells = plasma_cells[tls_mask]
print(f"Plasma cells in TLS-like regions: {len(tls_plasma_cells):,}")

# Cluster the TLS plasma cells to find distinct TLS regions
if len(tls_plasma_cells) > 0:
    from scipy.cluster.hierarchy import fcluster, linkage

    # Cluster TLS plasma cells
    if len(tls_plasma_cells) > 1:
        Z = linkage(tls_plasma_cells, method='average')
        cluster_labels = fcluster(Z, t=500, criterion='distance')  # 500px cluster distance
    else:
        cluster_labels = np.array([1])

    # Get cluster centers and sizes
    tls_regions = []
    for cluster_id in np.unique(cluster_labels):
        cluster_mask = cluster_labels == cluster_id
        cluster_points = tls_plasma_cells[cluster_mask]
        center_x = cluster_points[:, 0].mean()
        center_y = cluster_points[:, 1].mean()
        n_plasma_in_cluster = len(cluster_points)

        # Count lymphocytes near this cluster center
        lymph_near = lymph_tree.query_ball_point([center_x, center_y], 300, return_length=True)

        if n_plasma_in_cluster >= 5 and lymph_near >= 30:
            tls_regions.append({
                'center_x': center_x,
                'center_y': center_y,
                'n_plasma': n_plasma_in_cluster,
                'n_lymph_nearby': lymph_near,
                'radius': 300
            })

    print(f"\nFound {len(tls_regions)} potential TLS regions")

    # Sort by combined score (lymphocytes + plasma cells)
    tls_regions = sorted(tls_regions, key=lambda x: x['n_plasma'] + x['n_lymph_nearby'], reverse=True)

    # Print top TLS candidates
    print("\nTop TLS candidates:")
    for i, tls in enumerate(tls_regions[:10]):
        print(f"  TLS {i+1}: center=({tls['center_x']:.0f}, {tls['center_y']:.0f}), "
              f"plasma={tls['n_plasma']}, lymph={tls['n_lymph_nearby']}")

# Visualize TLS regions
cell_type_colors = {
    'Epithelial': '#1f77b4',
    'Connective': '#ff7f0e',
    'Lymphocyte': '#2ca02c',
    'Plasma cell': '#d62728',
    'Eosinophil': '#9467bd',
    'Neutrophil': '#8c564b',
}

# Create figure with top 4 TLS regions
n_show = min(4, len(tls_regions))
if n_show > 0:
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()

    for idx in range(4):
        ax = axes[idx]

        if idx < n_show:
            tls = tls_regions[idx]
            cx, cy = tls['center_x'], tls['center_y']
            crop_size = 800  # pixels

            x_min, x_max = cx - crop_size, cx + crop_size
            y_min, y_max = cy - crop_size, cy + crop_size

            # Filter points to region
            mask = (
                (points['x'] >= x_min) & (points['x'] <= x_max) &
                (points['y'] >= y_min) & (points['y'] <= y_max)
            )
            region_points = points[mask]

            # Plot each cell type
            for cell_type, color in cell_type_colors.items():
                ct_mask = region_points['classification_name'] == cell_type
                subset = region_points[ct_mask]
                if len(subset) > 0:
                    ax.scatter(
                        subset['x'], subset['y'],
                        c=color,
                        label=f'{cell_type} ({len(subset)})',
                        s=8 if cell_type in ['Lymphocyte', 'Plasma cell'] else 4,
                        alpha=0.8,
                        rasterized=True
                    )

            # Mark TLS center
            ax.scatter([cx], [cy], c='black', s=200, marker='x', linewidths=3, zorder=10)
            circle = plt.Circle((cx, cy), 300, fill=False, color='black', linewidth=2, linestyle='--')
            ax.add_patch(circle)

            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_max, y_min)
            ax.set_xlabel('X (pixels)')
            ax.set_ylabel('Y (pixels)')
            ax.set_title(f"TLS Candidate {idx+1}\nPlasma: {tls['n_plasma']}, Lymphocytes: {tls['n_lymph_nearby']}",
                        fontsize=12, fontweight='bold')
            ax.legend(loc='upper right', fontsize=8, markerscale=3)
            ax.set_aspect('equal')
        else:
            ax.axis('off')

    plt.suptitle('Tertiary Lymphoid Structure (TLS) Candidates\nPlasma cells + Lymphocyte clusters',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    output_path = "./tls_candidates.png"
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"\nPlot saved to {output_path}")
    plt.show()
else:
    print("No TLS regions found with current parameters")
