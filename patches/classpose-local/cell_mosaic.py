#!/usr/bin/env python
"""Create mosaic of example cells by type with segmentation boundaries."""

import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openslide
from matplotlib.patches import Polygon as MplPolygon
from PIL import Image

# Paths
slide_path = Path("test_data/CMU-1-Small-Region.svs")
geojson_path = Path("output/CMU-1-Small-Region_cell_contours.geojson")
output_path = Path("output/cell_mosaic.png")

# Parameters
CELLS_PER_TYPE = 5
ROI_SIZE = 50  # fixed square size for all ROIs (pixels)
MPP = 0.499  # microns per pixel from slide metadata

# Load slide
slide = openslide.OpenSlide(str(slide_path))
print(f"Slide dimensions: {slide.dimensions}")

# Load GeoJSON
with open(geojson_path) as f:
    geojson = json.load(f)

# Group cells by type with their coordinates
cells_by_type = {}
for feature in geojson['features']:
    cell_type = feature['properties']['classification']['name']
    color = feature['properties']['classification']['color']
    coords = np.array(feature['geometry']['coordinates'][0])
    centroid_x = feature['properties']['measurements'][2]['value']
    centroid_y = feature['properties']['measurements'][3]['value']
    area = feature['properties']['measurements'][0]['value']

    if cell_type not in cells_by_type:
        cells_by_type[cell_type] = {
            'color': tuple(c/255 for c in color),
            'cells': []
        }

    cells_by_type[cell_type]['cells'].append({
        'coords': coords,
        'centroid': (centroid_x, centroid_y),
        'area': area,
        'bbox': (coords[:, 0].min(), coords[:, 1].min(),
                 coords[:, 0].max(), coords[:, 1].max())
    })

print(f"Cell types: {list(cells_by_type.keys())}")

# Select 5 random cells per type (prefer medium-sized cells)
random.seed(42)
selected_cells = {}
for cell_type, data in cells_by_type.items():
    cells = data['cells']
    # Sort by area and pick from middle range for better examples
    cells_sorted = sorted(cells, key=lambda c: c['area'])
    mid_start = len(cells_sorted) // 4
    mid_end = 3 * len(cells_sorted) // 4
    mid_cells = cells_sorted[mid_start:mid_end] if mid_end > mid_start else cells_sorted

    n_select = min(CELLS_PER_TYPE, len(mid_cells))
    selected = random.sample(mid_cells, n_select)
    selected_cells[cell_type] = {
        'color': data['color'],
        'cells': selected
    }
    print(f"{cell_type}: selected {n_select} cells")

# Create mosaic
n_types = len(selected_cells)
fig, axes = plt.subplots(n_types, CELLS_PER_TYPE, figsize=(CELLS_PER_TYPE * 2, n_types * 2))

if n_types == 1:
    axes = [axes]

for row_idx, (cell_type, data) in enumerate(selected_cells.items()):
    color = data['color']

    for col_idx in range(CELLS_PER_TYPE):
        ax = axes[row_idx][col_idx] if n_types > 1 else axes[col_idx]

        if col_idx < len(data['cells']):
            cell = data['cells'][col_idx]
            coords = cell['coords']
            centroid = cell['centroid']

            # Calculate square ROI centered on cell centroid
            half_size = ROI_SIZE // 2
            x_center = int(centroid[0])
            y_center = int(centroid[1])

            x_min = x_center - half_size
            y_min = y_center - half_size

            # Ensure within slide bounds
            x_min = max(0, min(x_min, slide.dimensions[0] - ROI_SIZE))
            y_min = max(0, min(y_min, slide.dimensions[1] - ROI_SIZE))

            # Read square region from slide
            region = slide.read_region((x_min, y_min), 0, (ROI_SIZE, ROI_SIZE))
            region = region.convert('RGB')
            region_arr = np.array(region)

            # Plot image
            ax.imshow(region_arr)

            # Adjust coordinates to ROI space and draw polygon
            local_coords = coords - np.array([x_min, y_min])
            polygon = MplPolygon(local_coords,
                                fill=False,
                                edgecolor=color,
                                linewidth=2.5)
            ax.add_patch(polygon)

            # Add cell type label on first column
            if col_idx == 0:
                ax.set_ylabel(cell_type, fontsize=12, fontweight='bold',
                            color='black', rotation=0, ha='right', va='center')

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

# Add scale bar to bottom-right tile
scale_bar_um = 10  # 10 micron scale bar
scale_bar_px = scale_bar_um / MPP  # convert to pixels

# Get the last axis (bottom-right)
last_ax = axes[-1][-1] if n_types > 1 else axes[-1]

# Draw scale bar
bar_y = ROI_SIZE - 5  # near bottom
bar_x_start = ROI_SIZE - scale_bar_px - 5  # near right edge
last_ax.plot([bar_x_start, bar_x_start + scale_bar_px], [bar_y, bar_y],
             color='black', linewidth=2)
last_ax.text(bar_x_start + scale_bar_px/2, bar_y - 3, f'{scale_bar_um} µm',
             ha='center', va='bottom', fontsize=8, color='black')

# Add title
fig.suptitle('ClassPose Cell Detection Examples\n(5 cells per type with segmentation boundaries)',
             fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\nSaved mosaic to: {output_path}")

plt.show()
