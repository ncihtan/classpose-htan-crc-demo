#!/usr/bin/env python
"""Visualize ClassPose cell detections overlaid on slide thumbnail."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openslide
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

# Paths
slide_path = Path("test_data/CMU-1-Small-Region.svs")
geojson_path = Path("output/CMU-1-Small-Region_cell_contours.geojson")
output_path = Path("output/cell_visualization.png")

# Load slide and get thumbnail
slide = openslide.OpenSlide(str(slide_path))
dims = slide.dimensions  # (width, height)
print(f"Slide dimensions: {dims}")

# Get thumbnail at reasonable resolution
thumbnail_size = (800, int(800 * dims[1] / dims[0]))
thumbnail = slide.get_thumbnail(thumbnail_size)
thumbnail_arr = np.array(thumbnail)
print(f"Thumbnail size: {thumbnail_arr.shape}")

# Calculate scale factor
scale_x = thumbnail_size[0] / dims[0]
scale_y = thumbnail_size[1] / dims[1]

# Load GeoJSON
with open(geojson_path) as f:
    geojson = json.load(f)

print(f"Number of cells: {len(geojson['features'])}")

# Group cells by type
cells_by_type = {}
for feature in geojson['features']:
    cell_type = feature['properties']['classification']['name']
    color = feature['properties']['classification']['color']
    coords = feature['geometry']['coordinates'][0]

    if cell_type not in cells_by_type:
        cells_by_type[cell_type] = {
            'color': tuple(c/255 for c in color),
            'coords': []
        }
    cells_by_type[cell_type]['coords'].append(coords)

print(f"Cell types found: {list(cells_by_type.keys())}")
for ct, data in cells_by_type.items():
    print(f"  {ct}: {len(data['coords'])} cells")

# Create figure
fig, ax = plt.subplots(figsize=(12, 12 * dims[1] / dims[0]))

# Show thumbnail
ax.imshow(thumbnail_arr, extent=[0, dims[0], dims[1], 0])

# Plot cells by type
for cell_type, data in cells_by_type.items():
    patches = []
    for coords in data['coords']:
        poly = MplPolygon(coords, closed=True)
        patches.append(poly)

    collection = PatchCollection(
        patches,
        facecolor=(*data['color'], 0.5),  # Semi-transparent fill
        edgecolor=data['color'],
        linewidth=0.5
    )
    ax.add_collection(collection)

# Create legend
legend_handles = []
for cell_type, data in cells_by_type.items():
    handle = plt.Rectangle((0, 0), 1, 1,
                           facecolor=(*data['color'], 0.5),
                           edgecolor=data['color'],
                           label=f"{cell_type} ({len(data['coords'])})")
    legend_handles.append(handle)

ax.legend(handles=legend_handles, loc='upper right', fontsize=10)

ax.set_xlim(0, dims[0])
ax.set_ylim(dims[1], 0)  # Flip y-axis
ax.set_aspect('equal')
ax.set_title('ClassPose Cell Detection Results\nCMU-1-Small-Region.svs', fontsize=14)
ax.axis('off')

plt.tight_layout()
plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\nSaved visualization to: {output_path}")

plt.show()
