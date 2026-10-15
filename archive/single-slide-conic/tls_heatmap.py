import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

# Load the SpatialData object
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
sdata = sd.read_zarr(zarr_path)

# Get points with cell type info
print("Loading cell centroids...")
points = sdata.points['cell_centroids'].compute()
print(f"Total cells: {len(points):,}")

# Extract coordinates
all_x = points['x'].values
all_y = points['y'].values

# Extract lymphocytes and plasma cells
lymphocytes = points[points['classification_name'] == 'Lymphocyte'][['x', 'y']].values
plasma_cells = points[points['classification_name'] == 'Plasma cell'][['x', 'y']].values
print(f"Lymphocytes: {len(lymphocytes):,}")
print(f"Plasma cells: {len(plasma_cells):,}")

# Define grid for heatmap
x_min, x_max = all_x.min(), all_x.max()
y_min, y_max = all_y.min(), all_y.max()

# Grid resolution
grid_size = 200  # pixels per bin
n_bins_x = int(np.ceil((x_max - x_min) / grid_size))
n_bins_y = int(np.ceil((y_max - y_min) / grid_size))
print(f"\nGrid: {n_bins_x} x {n_bins_y} bins ({grid_size}px resolution)")

# Create 2D histograms for each cell type
lymph_hist, x_edges, y_edges = np.histogram2d(
    lymphocytes[:, 0], lymphocytes[:, 1],
    bins=[n_bins_x, n_bins_y],
    range=[[x_min, x_max], [y_min, y_max]]
)

plasma_hist, _, _ = np.histogram2d(
    plasma_cells[:, 0], plasma_cells[:, 1],
    bins=[n_bins_x, n_bins_y],
    range=[[x_min, x_max], [y_min, y_max]]
)

# Also create total cell density for normalization/masking
all_hist, _, _ = np.histogram2d(
    all_x, all_y,
    bins=[n_bins_x, n_bins_y],
    range=[[x_min, x_max], [y_min, y_max]]
)

# Apply Gaussian smoothing (rolling window effect)
sigma = 3  # bins
lymph_smooth = gaussian_filter(lymph_hist.astype(float), sigma=sigma)
plasma_smooth = gaussian_filter(plasma_hist.astype(float), sigma=sigma)
all_smooth = gaussian_filter(all_hist.astype(float), sigma=sigma)

# Compute TLS score: product of normalized densities
# This highlights regions where BOTH cell types are present
# Normalize by total cell density to get proportions
epsilon = 1e-10
lymph_prop = lymph_smooth / (all_smooth + epsilon)
plasma_prop = plasma_smooth / (all_smooth + epsilon)

# TLS score = geometric mean of proportions (highlights co-localization)
tls_score = np.sqrt(lymph_prop * plasma_prop)

# Mask regions with very few cells
tissue_mask = all_smooth > 5
tls_score_masked = np.ma.masked_where(~tissue_mask, tls_score)

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(18, 16))

# Plot extent
extent = [x_min, x_max, y_max, y_min]  # y inverted for image coords

# 1. Lymphocyte density
ax1 = axes[0, 0]
im1 = ax1.imshow(lymph_smooth.T, extent=extent, cmap='Greens', aspect='equal', origin='upper')
ax1.set_title('Lymphocyte Density', fontsize=14, fontweight='bold')
ax1.set_xlabel('X (pixels)')
ax1.set_ylabel('Y (pixels)')
plt.colorbar(im1, ax=ax1, label='Cell count (smoothed)')

# 2. Plasma cell density
ax2 = axes[0, 1]
im2 = ax2.imshow(plasma_smooth.T, extent=extent, cmap='Reds', aspect='equal', origin='upper')
ax2.set_title('Plasma Cell Density', fontsize=14, fontweight='bold')
ax2.set_xlabel('X (pixels)')
ax2.set_ylabel('Y (pixels)')
plt.colorbar(im2, ax=ax2, label='Cell count (smoothed)')

# 3. TLS Score Heatmap
ax3 = axes[1, 0]
im3 = ax3.imshow(tls_score_masked.T, extent=extent, cmap='hot', aspect='equal', origin='upper',
                  vmin=0, vmax=np.percentile(tls_score[tissue_mask], 99))
ax3.set_title('TLS Score (Plasma × Lymphocyte co-localization)', fontsize=14, fontweight='bold')
ax3.set_xlabel('X (pixels)')
ax3.set_ylabel('Y (pixels)')
plt.colorbar(im3, ax=ax3, label='TLS Score')

# 4. Overlay on cell types
ax4 = axes[1, 1]
# Plot all cells faintly
ax4.scatter(all_x, all_y, c='lightgray', s=0.01, alpha=0.3, rasterized=True)
# Overlay TLS heatmap
im4 = ax4.imshow(tls_score_masked.T, extent=extent, cmap='hot', aspect='equal', origin='upper',
                  alpha=0.7, vmin=0, vmax=np.percentile(tls_score[tissue_mask], 99))
ax4.set_title('TLS Heatmap Overlay on Tissue', fontsize=14, fontweight='bold')
ax4.set_xlabel('X (pixels)')
ax4.set_ylabel('Y (pixels)')
ax4.set_xlim(x_min, x_max)
ax4.set_ylim(y_max, y_min)
plt.colorbar(im4, ax=ax4, label='TLS Score')

plt.suptitle('Tertiary Lymphoid Structure (TLS) Spatial Heatmap\nRolling window analysis of Plasma + Lymphocyte co-localization',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()

output_path = "./tls_heatmap.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\nPlot saved to {output_path}")

# Find and report top TLS hotspots
print("\nTop TLS hotspots (grid coordinates):")
tls_flat = tls_score.flatten()
top_indices = np.argsort(tls_flat)[-10:][::-1]
for i, idx in enumerate(top_indices):
    grid_x, grid_y = np.unravel_index(idx, tls_score.shape)
    real_x = x_min + (grid_x + 0.5) * grid_size
    real_y = y_min + (grid_y + 0.5) * grid_size
    score = tls_flat[idx]
    print(f"  {i+1}. ({real_x:.0f}, {real_y:.0f}) - TLS score: {score:.4f}")

plt.show()
