import spatialdata as sd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import numpy as np
import openslide

# Pixel scale
UM_PER_PIXEL = 0.325
FOV_UM = 500  # 500 um field of view
FOV_PX = int(FOV_UM / UM_PER_PIXEL)  # ~1538 pixels

print(f"FOV: {FOV_UM} um = {FOV_PX} pixels")

# Base paths
BASE_DIR = "."
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"

# Data paths for both models
conic_zarr_path = f"{BASE_DIR}/{SAMPLE_ID}/{SAMPLE_ID}_spatialdata.zarr"
consep_zarr_path = f"{BASE_DIR}/dual_model_output/consep/{SAMPLE_ID}_consep_spatialdata.zarr"
he_slide_path = f"{BASE_DIR}/{SAMPLE_ID}.tif"

# Cell type color palettes
conic_colors = {
    'Epithelial': '#1f77b4',    # blue
    'Connective': '#ff7f0e',    # orange
    'Lymphocyte': '#2ca02c',    # green
    'Plasma cell': '#d62728',   # red
    'Eosinophil': '#9467bd',    # purple
    'Neutrophil': '#8c564b',    # brown
}

consep_colors = {
    'Healthy epithelial': '#3498db',  # light blue
    'Malignant epithelial': '#e74c3c',  # red
    'Inflammatory': '#2ecc71',  # green
    'Stroma': '#f39c12',  # orange
    'Muscle': '#9b59b6',  # purple
    'Other': '#95a5a6',  # gray
}

# Load H&E slide
print("Loading H&E slide...")
he_slide = openslide.OpenSlide(he_slide_path)
slide_dims = he_slide.dimensions
print(f"Slide dimensions: {slide_dims}")

# Load both SpatialData objects
print("\nLoading CoNIC data...")
sdata_conic = sd.read_zarr(conic_zarr_path)
contours_conic = sdata_conic.shapes['cell_contours']
print(f"CoNIC total contours: {len(contours_conic):,}")

print("\nLoading ConSep data...")
sdata_consep = sd.read_zarr(consep_zarr_path)
contours_consep = sdata_consep.shapes['cell_contours']
print(f"ConSep total contours: {len(contours_consep):,}")

# Define regions of interest (center_x, center_y, description)
regions = [
    (58000, 18000, "Immune infiltration zone"),
    (10000, 40000, "Epithelial dominant: Glandular region"),
    (65000, 40000, "Stromal region"),
    (40000, 30000, "Tumor-stroma interface"),
]

region_colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']


def plot_contours_for_model(ax, contours, cell_type_colors, x_min, x_max, y_min, y_max, title):
    """Plot cell contours for a given model within the specified region."""
    # Filter contours to region using centroid coordinates
    centroids = contours.geometry.centroid
    mask = (
        (centroids.x >= x_min) & (centroids.x <= x_max) &
        (centroids.y >= y_min) & (centroids.y <= y_max)
    )
    contours_in_region = contours[mask]

    print(f"\n  {title}: {len(contours_in_region):,} contours")

    # Plot contours by cell type
    for cell_type, color in cell_type_colors.items():
        ct_contours = contours_in_region[contours_in_region['classification_name'] == cell_type]
        n_cells = len(ct_contours)

        if n_cells > 0:
            pct = 100 * n_cells / len(contours_in_region) if len(contours_in_region) > 0 else 0
            print(f"    {cell_type}: {n_cells:,} ({pct:.1f}%)")

            # Create polygon patches for each contour
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
                    patches,
                    facecolor='none',
                    edgecolor=color,
                    alpha=1.0,
                    linewidth=0.5
                )
                ax.add_collection(collection)
                ax.scatter([], [], c=color, label=f'{cell_type} ({pct:.0f}%)', s=50)

    return len(contours_in_region)


# Create figure: 4 rows (one per region), 3 columns (overview, conic, consep)
fig = plt.figure(figsize=(22, 28))

# Create overview panel at top spanning full width
ax_overview = fig.add_subplot(5, 3, (1, 3))

# Get H&E thumbnail for overview
thumbnail_size = (1800, int(1800 * slide_dims[1] / slide_dims[0]))
he_thumbnail = he_slide.get_thumbnail(thumbnail_size)
he_thumbnail_arr = np.array(he_thumbnail.convert('RGB'))
ax_overview.imshow(he_thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0])

# Draw boxes for each crop region
for idx, ((cx, cy, title), color) in enumerate(zip(regions, region_colors)):
    half_fov = FOV_PX // 2
    x_min, x_max = cx - half_fov, cx + half_fov
    y_min, y_max = cy - half_fov, cy + half_fov

    rect = plt.Rectangle((x_min, y_min), FOV_PX, FOV_PX,
                         fill=False, edgecolor=color, linewidth=3)
    ax_overview.add_patch(rect)
    ax_overview.text(cx, y_min - 500, f"{idx+1}", fontsize=14, fontweight='bold',
                    color=color, ha='center', va='bottom')

ax_overview.set_xlim(0, slide_dims[0])
ax_overview.set_ylim(slide_dims[1], 0)
ax_overview.set_title('Overview - Crop Locations', fontsize=14, fontweight='bold')
ax_overview.axis('off')

# Create crop panels: 4 regions x 2 models
for idx, ((cx, cy, title), region_color) in enumerate(zip(regions, region_colors)):
    print(f"\n{'='*60}")
    print(f"Region {idx+1}: {title}")
    print(f"  Center: ({cx}, {cy}), FOV: {FOV_UM} um")

    # Calculate crop bounds
    half_fov = FOV_PX // 2
    x_min, x_max = cx - half_fov, cx + half_fov
    y_min, y_max = cy - half_fov, cy + half_fov

    # Read H&E crop once
    he_region = he_slide.read_region((x_min, y_min), 0, (FOV_PX, FOV_PX))
    he_crop = np.array(he_region.convert('RGB'))

    # CoNIC panel (left)
    ax_conic = fig.add_subplot(5, 2, 3 + idx * 2)
    ax_conic.imshow(he_crop, extent=[x_min, x_max, y_max, y_min], alpha=1.0)
    n_conic = plot_contours_for_model(ax_conic, contours_conic, conic_colors,
                                       x_min, x_max, y_min, y_max, "CoNIC")
    ax_conic.set_xlim(x_min, x_max)
    ax_conic.set_ylim(y_max, y_min)
    ax_conic.set_xlabel('X (pixels)')
    ax_conic.set_ylabel('Y (pixels)')
    ax_conic.set_title(f"{idx+1}. {title}\nCoNIC ({n_conic:,} cells)", fontsize=11, fontweight='bold', color=region_color)
    ax_conic.legend(loc='upper right', fontsize=7)
    ax_conic.set_aspect('equal')

    # ConSep panel (right)
    ax_consep = fig.add_subplot(5, 2, 4 + idx * 2)
    ax_consep.imshow(he_crop, extent=[x_min, x_max, y_max, y_min], alpha=1.0)
    n_consep = plot_contours_for_model(ax_consep, contours_consep, consep_colors,
                                        x_min, x_max, y_min, y_max, "ConSep")
    ax_consep.set_xlim(x_min, x_max)
    ax_consep.set_ylim(y_max, y_min)
    ax_consep.set_xlabel('X (pixels)')
    ax_consep.set_ylabel('Y (pixels)')
    ax_consep.set_title(f"{idx+1}. {title}\nConSep ({n_consep:,} cells)", fontsize=11, fontweight='bold', color=region_color)
    ax_consep.legend(loc='upper right', fontsize=7)
    ax_consep.set_aspect('equal')

plt.suptitle(f'Cell Type Comparison - CoNIC vs ConSep\nSample {SAMPLE_ID.split("_")[-1]} ({FOV_UM} um FOV)',
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
output_path = f"{BASE_DIR}/cell_types_crops.png"
plt.savefig(output_path, dpi=200, bbox_inches='tight')
print(f"\nPlot saved to {output_path}")
plt.show()
