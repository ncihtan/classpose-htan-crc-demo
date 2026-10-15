import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
import networkx as nx
from collections import Counter

# Load the SpatialData object
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

# Get points with cell type info
print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()
print(f"Total cells: {len(points):,}")

# Filter to immune cells only (lymphocytes + plasma cells) for TLS analysis
immune_mask = points['classification_name'].isin(['Lymphocyte', 'Plasma cell'])
immune_cells = points[immune_mask].reset_index(drop=True)
print(f"Immune cells (Lymphocytes + Plasma): {len(immune_cells):,}")

# Extract coordinates and cell types
coords = immune_cells[['x', 'y']].values
cell_types = immune_cells['classification_name'].values

# Build spatial graph using radius neighbors
print("\nBuilding spatial graph...")
radius = 50  # pixels - connect cells within this distance
tree = cKDTree(coords)
pairs = tree.query_pairs(r=radius, output_type='ndarray')
print(f"Found {len(pairs):,} edges (radius={radius}px)")

# Create graph
G = nx.Graph()
G.add_nodes_from(range(len(immune_cells)))

# Add node attributes
for i, (ct, x, y) in enumerate(zip(cell_types, coords[:, 0], coords[:, 1])):
    G.nodes[i]['cell_type'] = ct
    G.nodes[i]['x'] = x
    G.nodes[i]['y'] = y

# Add edges
G.add_edges_from(pairs)
print(f"Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

# Find connected components (immune cell clusters)
print("\nFinding connected components...")
components = list(nx.connected_components(G))
print(f"Found {len(components):,} connected components")

# Analyze each component for TLS characteristics
tls_candidates = []
for comp_id, comp in enumerate(components):
    if len(comp) < 10:  # Skip tiny clusters
        continue

    # Count cell types in component
    comp_types = [G.nodes[n]['cell_type'] for n in comp]
    type_counts = Counter(comp_types)

    n_lymph = type_counts.get('Lymphocyte', 0)
    n_plasma = type_counts.get('Plasma cell', 0)

    # TLS criteria: need both cell types present
    if n_lymph >= 5 and n_plasma >= 3:
        # Calculate component centroid
        comp_coords = np.array([[G.nodes[n]['x'], G.nodes[n]['y']] for n in comp])
        centroid_x = comp_coords[:, 0].mean()
        centroid_y = comp_coords[:, 1].mean()

        # Calculate component density (cells per area)
        if len(comp) > 1:
            spread = np.sqrt(np.var(comp_coords[:, 0]) + np.var(comp_coords[:, 1]))
        else:
            spread = 1

        # TLS score based on:
        # - Total immune cells
        # - Balance of lymphocytes and plasma cells
        # - Density
        balance = min(n_lymph, n_plasma * 3) / max(n_lymph, n_plasma * 3)  # plasma cells are rarer
        tls_score = (n_lymph + n_plasma * 2) * balance  # weight plasma cells more

        tls_candidates.append({
            'component_id': comp_id,
            'nodes': comp,
            'n_cells': len(comp),
            'n_lymph': n_lymph,
            'n_plasma': n_plasma,
            'centroid_x': centroid_x,
            'centroid_y': centroid_y,
            'spread': spread,
            'tls_score': tls_score
        })

# Sort by TLS score
tls_candidates = sorted(tls_candidates, key=lambda x: x['tls_score'], reverse=True)
print(f"\nFound {len(tls_candidates)} TLS candidate clusters")

print("\nTop 10 TLS candidates:")
for i, tls in enumerate(tls_candidates[:10]):
    print(f"  {i+1}. Cells: {tls['n_cells']}, Lymph: {tls['n_lymph']}, Plasma: {tls['n_plasma']}, "
          f"Score: {tls['tls_score']:.1f}, Center: ({tls['centroid_x']:.0f}, {tls['centroid_y']:.0f})")

# Visualization
print("\nCreating visualizations...")

# Plot top 4 TLS clusters
fig, axes = plt.subplots(2, 2, figsize=(16, 16))
axes = axes.flatten()

colors = {'Lymphocyte': '#2ca02c', 'Plasma cell': '#d62728'}

for idx in range(4):
    ax = axes[idx]

    if idx < len(tls_candidates):
        tls = tls_candidates[idx]
        cx, cy = tls['centroid_x'], tls['centroid_y']

        # Get nodes in this component
        comp_nodes = list(tls['nodes'])

        # Create subgraph
        subgraph = G.subgraph(comp_nodes)

        # Get positions and colors
        pos = {n: (G.nodes[n]['x'], G.nodes[n]['y']) for n in comp_nodes}
        node_colors = [colors[G.nodes[n]['cell_type']] for n in comp_nodes]

        # Draw the graph
        nx.draw_networkx_edges(subgraph, pos, ax=ax, alpha=0.3, edge_color='gray', width=0.5)
        nx.draw_networkx_nodes(subgraph, pos, ax=ax, node_color=node_colors,
                               node_size=30, alpha=0.8)

        # Set limits with padding
        xs = [pos[n][0] for n in comp_nodes]
        ys = [pos[n][1] for n in comp_nodes]
        pad = 100
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(max(ys) + pad, min(ys) - pad)  # inverted y

        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
        ax.set_title(f"TLS Cluster {idx+1}\n{tls['n_cells']} cells: {tls['n_lymph']} Lymph (green) + {tls['n_plasma']} Plasma (red)\nScore: {tls['tls_score']:.1f}",
                    fontsize=11, fontweight='bold')
        ax.set_aspect('equal')
    else:
        ax.axis('off')

plt.suptitle('TLS Detection via Graph Analysis\nConnected components of immune cells (radius=50px)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
output_path = "./tls_graph_clusters.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path}")

# Create overview plot showing all TLS locations
fig, ax = plt.subplots(figsize=(14, 14))

# Plot all cells faintly
all_x = points['x'].values
all_y = points['y'].values
ax.scatter(all_x, all_y, c='lightgray', s=0.01, alpha=0.2, rasterized=True)

# Highlight TLS clusters
for i, tls in enumerate(tls_candidates[:50]):  # Top 50
    cx, cy = tls['centroid_x'], tls['centroid_y']
    size = np.sqrt(tls['n_cells']) * 20
    alpha = min(0.8, 0.3 + tls['tls_score'] / 100)

    circle = plt.Circle((cx, cy), tls['spread'] + 50, fill=False,
                        color='red', linewidth=2, alpha=alpha)
    ax.add_patch(circle)

    if i < 10:
        ax.annotate(f"{i+1}", (cx, cy), fontsize=10, fontweight='bold',
                   ha='center', va='center', color='red')

ax.set_xlim(all_x.min(), all_x.max())
ax.set_ylim(all_y.max(), all_y.min())
ax.set_xlabel('X (pixels)')
ax.set_ylabel('Y (pixels)')
ax.set_title('TLS Locations (Graph-based detection)\nTop 50 immune cell clusters (numbered 1-10 are highest scoring)',
            fontsize=14, fontweight='bold')
ax.set_aspect('equal')

plt.tight_layout()
output_path2 = "./tls_graph_overview.png"
plt.savefig(output_path2, dpi=200, bbox_inches='tight')
print(f"Saved: {output_path2}")

plt.show()
