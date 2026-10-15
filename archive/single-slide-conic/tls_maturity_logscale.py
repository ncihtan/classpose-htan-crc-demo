import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
import networkx as nx
from collections import Counter
from scipy import stats

# Load the SpatialData object
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()

# Filter to immune cells
immune_mask = points['classification_name'].isin(['Lymphocyte', 'Plasma cell'])
immune_cells = points[immune_mask].reset_index(drop=True)

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

# Analyze TLS
tls_data = []
for comp_id, comp in enumerate(components):
    if len(comp) < 15:
        continue

    comp_types = [G.nodes[n]['cell_type'] for n in comp]
    type_counts = Counter(comp_types)
    n_lymph = type_counts.get('Lymphocyte', 0)
    n_plasma = type_counts.get('Plasma cell', 0)

    if n_lymph >= 10 and n_plasma >= 5:
        comp_coords = np.array([[G.nodes[n]['x'], G.nodes[n]['y']] for n in comp])
        maturity_ratio = n_plasma / n_lymph

        # Categorize
        if maturity_ratio < 0.15:
            stage = 'Early'
        elif maturity_ratio < 0.25:
            stage = 'Intermediate'
        else:
            stage = 'Mature'

        tls_data.append({
            'n_cells': len(comp),
            'n_lymph': n_lymph,
            'n_plasma': n_plasma,
            'maturity_ratio': maturity_ratio,
            'maturity_stage': stage
        })

print(f"Found {len(tls_data)} TLS clusters")

# Extract data for plotting
sizes = np.array([t['n_cells'] for t in tls_data])
ratios = np.array([t['maturity_ratio'] for t in tls_data])
stages = [t['maturity_stage'] for t in tls_data]

# Calculate correlation
log_sizes = np.log10(sizes)
corr, pval = stats.spearmanr(sizes, ratios)
print(f"\nSpearman correlation (size vs maturity): r={corr:.3f}, p={pval:.2e}")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

colors_stage = {'Early': '#2ca02c', 'Intermediate': '#ff7f0e', 'Mature': '#d62728'}

# Plot 1: Log scale on x-axis
ax1 = axes[0]
for stage in ['Early', 'Intermediate', 'Mature']:
    mask = np.array([s == stage for s in stages])
    ax1.scatter(sizes[mask], ratios[mask], c=colors_stage[stage],
                label=stage, alpha=0.6, s=60, edgecolors='white', linewidths=0.5)

ax1.set_xscale('log')
ax1.set_xlabel('TLS Size (total cells) - log scale', fontsize=12)
ax1.set_ylabel('Maturity Ratio (Plasma/Lymphocyte)', fontsize=12)
ax1.set_title('TLS Size vs Maturity (log scale)', fontsize=14, fontweight='bold')
ax1.axhline(0.15, color='orange', linestyle='--', alpha=0.5, label='Early/Intermediate')
ax1.axhline(0.25, color='red', linestyle='--', alpha=0.5, label='Intermediate/Mature')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)

# Add correlation annotation
ax1.text(0.05, 0.95, f'Spearman r = {corr:.3f}\np = {pval:.2e}',
         transform=ax1.transAxes, fontsize=11, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 2: Log-log scale
ax2 = axes[1]
for stage in ['Early', 'Intermediate', 'Mature']:
    mask = np.array([s == stage for s in stages])
    ax2.scatter(sizes[mask], ratios[mask], c=colors_stage[stage],
                label=stage, alpha=0.6, s=60, edgecolors='white', linewidths=0.5)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('TLS Size (total cells) - log scale', fontsize=12)
ax2.set_ylabel('Maturity Ratio (Plasma/Lymphocyte) - log scale', fontsize=12)
ax2.set_title('TLS Size vs Maturity (log-log scale)', fontsize=14, fontweight='bold')
ax2.axhline(0.15, color='orange', linestyle='--', alpha=0.5)
ax2.axhline(0.25, color='red', linestyle='--', alpha=0.5)
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)

# Fit trend line on log-log
slope, intercept, r_value, p_value, std_err = stats.linregress(log_sizes, np.log10(ratios))
x_fit = np.logspace(np.log10(sizes.min()), np.log10(sizes.max()), 100)
y_fit = 10**(intercept + slope * np.log10(x_fit))
ax2.plot(x_fit, y_fit, 'k--', linewidth=2, alpha=0.7, label=f'Trend (slope={slope:.2f})')
ax2.legend(loc='upper right')

ax2.text(0.05, 0.95, f'Log-log slope = {slope:.3f}\nR² = {r_value**2:.3f}',
         transform=ax2.transAxes, fontsize=11, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('TLS Size vs Maturity Relationship\nSmaller TLS tend to be more mature (higher plasma/lymphocyte ratio)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = "./tls_size_maturity_log.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\nSaved: {output_path}")

plt.show()
