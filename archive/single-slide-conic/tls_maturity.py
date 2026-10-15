import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
import networkx as nx
from collections import Counter

# Load the SpatialData object
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()

# Filter to immune cells
immune_mask = points['classification_name'].isin(['Lymphocyte', 'Plasma cell'])
immune_cells = points[immune_mask].reset_index(drop=True)
print(f"Immune cells: {len(immune_cells):,}")

# Build spatial graph
coords = immune_cells[['x', 'y']].values
cell_types = immune_cells['classification_name'].values

radius = 50
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

# Analyze TLS with maturity scoring
tls_data = []
for comp_id, comp in enumerate(components):
    if len(comp) < 15:  # Need minimum size for reliable ratio
        continue

    comp_types = [G.nodes[n]['cell_type'] for n in comp]
    type_counts = Counter(comp_types)
    n_lymph = type_counts.get('Lymphocyte', 0)
    n_plasma = type_counts.get('Plasma cell', 0)

    # Need both cell types
    if n_lymph >= 10 and n_plasma >= 5:
        comp_coords = np.array([[G.nodes[n]['x'], G.nodes[n]['y']] for n in comp])
        centroid_x = comp_coords[:, 0].mean()
        centroid_y = comp_coords[:, 1].mean()
        spread = np.sqrt(np.var(comp_coords[:, 0]) + np.var(comp_coords[:, 1]))

        # Maturity ratio: plasma / lymphocyte
        # Higher ratio = more mature TLS (more B cells differentiated to plasma)
        maturity_ratio = n_plasma / n_lymph

        tls_data.append({
            'component_id': comp_id,
            'nodes': comp,
            'n_cells': len(comp),
            'n_lymph': n_lymph,
            'n_plasma': n_plasma,
            'maturity_ratio': maturity_ratio,
            'centroid_x': centroid_x,
            'centroid_y': centroid_y,
            'spread': spread
        })

print(f"\nFound {len(tls_data)} TLS clusters")

# Calculate maturity statistics
ratios = [t['maturity_ratio'] for t in tls_data]
print(f"\nMaturity ratio (Plasma/Lymphocyte) statistics:")
print(f"  Min: {min(ratios):.3f}")
print(f"  Max: {max(ratios):.3f}")
print(f"  Mean: {np.mean(ratios):.3f}")
print(f"  Median: {np.median(ratios):.3f}")

# Categorize TLS by maturity
for tls in tls_data:
    r = tls['maturity_ratio']
    if r < 0.15:
        tls['maturity_stage'] = 'Early'
    elif r < 0.25:
        tls['maturity_stage'] = 'Intermediate'
    else:
        tls['maturity_stage'] = 'Mature'

# Count by stage
stage_counts = Counter([t['maturity_stage'] for t in tls_data])
print(f"\nTLS by maturity stage:")
for stage in ['Early', 'Intermediate', 'Mature']:
    print(f"  {stage}: {stage_counts[stage]}")

# Sort by maturity ratio
tls_by_maturity = sorted(tls_data, key=lambda x: x['maturity_ratio'], reverse=True)

print("\nTop 10 MOST MATURE TLS (highest plasma/lymphocyte ratio):")
for i, tls in enumerate(tls_by_maturity[:10]):
    print(f"  {i+1}. Ratio: {tls['maturity_ratio']:.3f} | "
          f"Cells: {tls['n_cells']} (L:{tls['n_lymph']}, P:{tls['n_plasma']}) | "
          f"Loc: ({tls['centroid_x']:.0f}, {tls['centroid_y']:.0f})")

print("\nTop 10 LEAST MATURE TLS (lowest plasma/lymphocyte ratio):")
tls_by_immaturity = sorted(tls_data, key=lambda x: x['maturity_ratio'])
for i, tls in enumerate(tls_by_immaturity[:10]):
    print(f"  {i+1}. Ratio: {tls['maturity_ratio']:.3f} | "
          f"Cells: {tls['n_cells']} (L:{tls['n_lymph']}, P:{tls['n_plasma']}) | "
          f"Loc: ({tls['centroid_x']:.0f}, {tls['centroid_y']:.0f})")

# Visualization
fig = plt.figure(figsize=(18, 14))

# 1. Maturity ratio distribution
ax1 = fig.add_subplot(2, 2, 1)
ax1.hist(ratios, bins=30, edgecolor='black', alpha=0.7, color='steelblue')
ax1.axvline(0.15, color='orange', linestyle='--', linewidth=2, label='Early/Intermediate threshold')
ax1.axvline(0.25, color='red', linestyle='--', linewidth=2, label='Intermediate/Mature threshold')
ax1.set_xlabel('Maturity Ratio (Plasma/Lymphocyte)', fontsize=12)
ax1.set_ylabel('Number of TLS', fontsize=12)
ax1.set_title('TLS Maturity Distribution', fontsize=14, fontweight='bold')
ax1.legend()

# 2. Scatter: Size vs Maturity
ax2 = fig.add_subplot(2, 2, 2)
sizes = [t['n_cells'] for t in tls_data]
colors_stage = {'Early': 'green', 'Intermediate': 'orange', 'Mature': 'red'}
c = [colors_stage[t['maturity_stage']] for t in tls_data]
ax2.scatter(sizes, ratios, c=c, alpha=0.6, s=50)
ax2.set_xlabel('TLS Size (total cells)', fontsize=12)
ax2.set_ylabel('Maturity Ratio (Plasma/Lymphocyte)', fontsize=12)
ax2.set_title('TLS Size vs Maturity', fontsize=14, fontweight='bold')
ax2.axhline(0.15, color='orange', linestyle='--', alpha=0.5)
ax2.axhline(0.25, color='red', linestyle='--', alpha=0.5)
# Legend
for stage, col in colors_stage.items():
    ax2.scatter([], [], c=col, label=stage, s=100)
ax2.legend(title='Maturity Stage')

# 3. Spatial map colored by maturity
ax3 = fig.add_subplot(2, 2, 3)
all_x = points['x'].values
all_y = points['y'].values
ax3.scatter(all_x, all_y, c='lightgray', s=0.01, alpha=0.2, rasterized=True)

# Plot TLS colored by maturity
for tls in tls_data:
    cx, cy = tls['centroid_x'], tls['centroid_y']
    color = colors_stage[tls['maturity_stage']]
    size = np.sqrt(tls['n_cells']) * 3
    ax3.scatter(cx, cy, c=color, s=size, alpha=0.6, edgecolors='black', linewidths=0.5)

ax3.set_xlim(all_x.min(), all_x.max())
ax3.set_ylim(all_y.max(), all_y.min())
ax3.set_xlabel('X (pixels)')
ax3.set_ylabel('Y (pixels)')
ax3.set_title('TLS Spatial Distribution by Maturity Stage', fontsize=14, fontweight='bold')
for stage, col in colors_stage.items():
    ax3.scatter([], [], c=col, label=f'{stage} ({stage_counts[stage]})', s=100)
ax3.legend(loc='upper right')
ax3.set_aspect('equal')

# 4. Compare mature vs early TLS structure
ax4 = fig.add_subplot(2, 2, 4)

# Get one mature and one early TLS of similar size for comparison
mature_tls = [t for t in tls_data if t['maturity_stage'] == 'Mature' and t['n_cells'] > 100]
early_tls = [t for t in tls_data if t['maturity_stage'] == 'Early' and t['n_cells'] > 100]

if mature_tls and early_tls:
    # Pick examples
    mature_ex = sorted(mature_tls, key=lambda x: x['maturity_ratio'], reverse=True)[0]
    early_ex = sorted(early_tls, key=lambda x: x['n_cells'], reverse=True)[0]

    # Plot both side by side in the same axes with offset
    offset = 3000  # pixels offset between the two

    colors_cell = {'Lymphocyte': '#2ca02c', 'Plasma cell': '#d62728'}

    for tls, x_off, label in [(early_ex, 0, 'Early'), (mature_ex, offset, 'Mature')]:
        comp_nodes = list(tls['nodes'])
        subgraph = G.subgraph(comp_nodes)
        pos = {n: (G.nodes[n]['x'] - tls['centroid_x'] + x_off,
                   G.nodes[n]['y'] - tls['centroid_y']) for n in comp_nodes}
        node_colors = [colors_cell[G.nodes[n]['cell_type']] for n in comp_nodes]

        nx.draw_networkx_edges(subgraph, pos, ax=ax4, alpha=0.2, edge_color='gray', width=0.3)
        nx.draw_networkx_nodes(subgraph, pos, ax=ax4, node_color=node_colors, node_size=15, alpha=0.8)

        # Label
        ax4.text(x_off, -tls['spread'] - 100, f"{label}\nRatio: {tls['maturity_ratio']:.2f}\n({tls['n_lymph']}L + {tls['n_plasma']}P)",
                ha='center', fontsize=10, fontweight='bold')

    ax4.set_xlim(-1500, offset + 1500)
    ax4.set_aspect('equal')
    ax4.axis('off')
    ax4.set_title('Early vs Mature TLS Structure\nGreen=Lymphocyte, Red=Plasma cell', fontsize=14, fontweight='bold')

plt.suptitle('TLS Maturity Analysis\nPlasma/Lymphocyte ratio as proxy for B-cell differentiation',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = "./tls_maturity.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\nSaved: {output_path}")

plt.show()
