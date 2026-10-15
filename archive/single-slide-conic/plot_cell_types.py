import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
import openslide

# Base paths
BASE_DIR = "."
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"

# Data paths for both models
conic_zarr_path = f"{BASE_DIR}/{SAMPLE_ID}/{SAMPLE_ID}_spatialdata.zarr"
consep_zarr_path = f"{BASE_DIR}/dual_model_output/consep/{SAMPLE_ID}_consep_spatialdata.zarr"
he_slide_path = f"{BASE_DIR}/{SAMPLE_ID}.tif"

# Cell type color palettes
conic_colors = {
    "Epithelial": "#1f77b4",  # blue
    "Connective": "#ff7f0e",  # orange
    "Lymphocyte": "#2ca02c",  # green
    "Plasma cell": "#d62728",  # red
    "Eosinophil": "#9467bd",  # purple
    "Neutrophil": "#8c564b",  # brown
}

consep_colors = {
    "Healthy epithelial": "#3498db",  # light blue
    "Malignant epithelial": "#e74c3c",  # red
    "Inflammatory": "#2ecc71",  # green
    "Stroma": "#f39c12",  # orange
    "Muscle": "#9b59b6",  # purple
    "Other": "#95a5a6",  # gray
}

# Load H&E slide
print("Loading H&E slide...")
he_slide = openslide.OpenSlide(he_slide_path)
slide_dims = he_slide.dimensions
print(f"Slide dimensions: {slide_dims}")

# Load both SpatialData objects
print("\nLoading CoNIC data...")
sdata_conic = sd.read_zarr(conic_zarr_path)
points_conic = sdata_conic.points["cell_centroids"].compute()
print(f"CoNIC total cells: {len(points_conic):,}")
print(f"CoNIC cell classes: {points_conic['classification_name'].unique().tolist()}")

print("\nLoading ConSep data...")
sdata_consep = sd.read_zarr(consep_zarr_path)
points_consep = sdata_consep.points["cell_centroids"].compute()
print(f"ConSep total cells: {len(points_consep):,}")
print(f"ConSep cell classes: {points_consep['classification_name'].unique().tolist()}")

# Create side-by-side plot
print("\nCreating dual plot...")
fig, axes = plt.subplots(1, 2, figsize=(24, 12))

# Get H&E thumbnail
thumbnail_size = (1200, int(1200 * slide_dims[1] / slide_dims[0]))
he_thumbnail = he_slide.get_thumbnail(thumbnail_size)
he_thumbnail_arr = np.array(he_thumbnail.convert('RGB'))

# Plot CoNIC (left panel)
ax_conic = axes[0]
ax_conic.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=0.6)

print("\nPlotting CoNIC cell types:")
for cell_type, color in conic_colors.items():
    mask = points_conic["classification_name"] == cell_type
    subset = points_conic[mask]
    n_cells = len(subset)
    print(f"  {cell_type}: {n_cells:,} cells")
    ax_conic.scatter(
        subset["x"],
        subset["y"],
        c=color,
        label=f"{cell_type} ({n_cells:,})",
        s=0.001,
        alpha=0.6,
        rasterized=True,
    )

ax_conic.set_xlim(0, slide_dims[0])
ax_conic.set_ylim(slide_dims[1], 0)
ax_conic.set_xlabel("X (pixels)")
ax_conic.set_ylabel("Y (pixels)")
ax_conic.set_title(f"CoNIC Model - {len(points_conic):,} cells\n(6 cell types: immune-focused)", fontsize=14, fontweight='bold')
ax_conic.legend(markerscale=30, loc="upper right", fontsize=9)
ax_conic.set_aspect("equal")

# Plot ConSep (right panel)
ax_consep = axes[1]
ax_consep.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=0.6)

print("\nPlotting ConSep cell types:")
for cell_type, color in consep_colors.items():
    mask = points_consep["classification_name"] == cell_type
    subset = points_consep[mask]
    n_cells = len(subset)
    print(f"  {cell_type}: {n_cells:,} cells")
    ax_consep.scatter(
        subset["x"],
        subset["y"],
        c=color,
        label=f"{cell_type} ({n_cells:,})",
        s=0.001,
        alpha=0.6,
        rasterized=True,
    )

ax_consep.set_xlim(0, slide_dims[0])
ax_consep.set_ylim(slide_dims[1], 0)
ax_consep.set_xlabel("X (pixels)")
ax_consep.set_ylabel("Y (pixels)")
ax_consep.set_title(f"ConSep Model - {len(points_consep):,} cells\n(6 cell types: epithelial-focused)", fontsize=14, fontweight='bold')
ax_consep.legend(markerscale=30, loc="upper right", fontsize=9)
ax_consep.set_aspect("equal")

plt.suptitle(f"Cell Type Comparison - Sample {SAMPLE_ID.split('_')[-1]}", fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
output_path = f"{BASE_DIR}/cell_types_spatial_plot.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"\nPlot saved to {output_path}")
plt.show()
