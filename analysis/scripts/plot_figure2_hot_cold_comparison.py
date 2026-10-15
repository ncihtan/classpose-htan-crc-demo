#!/usr/bin/env python3
"""
Figure 2: Immune-Hot vs Immune-Cold Sample Comparison

Visualizes both extreme samples from the cohort:
- HTA7_917 (immune-hot): 38.3% inflammatory
- HTA7_932 (immune-cold): 6.0% inflammatory

Layout:
- Row 1: Full H&E overviews (no contours) with bounding boxes
- Row 2: Full H&E overviews with filled cell contours
- Row 3: Zoom to dense malignant region
- Row 4: Zoom to dense inflammatory region
- Row 5: Zoom to dense stromal region
"""

import spatialdata as sd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle
import numpy as np
import openslide
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import UM_PER_PIXEL, CONSEP_COLORS, CONSEP_ORDER, add_scale_bar, get_project_root

# =============================================================================
# Configuration
# =============================================================================

# Sample definitions
SAMPLES = {
    "hot": {
        "sample_id": "dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF",
        "participant": "HTA7_917",
        "label": "Immune-Hot (38% Inflammatory)",
    },
    "cold": {
        "sample_id": "dg_4DFC_98B41CE1-C47F-4C37-BB1D-F4B49E1C4619",
        "participant": "HTA7_932",
        "label": "Immune-Cold (6% Inflammatory)",
    },
}

DATA_DIR = SAMPLE_DATA_DIR

# Field of view for zoom regions
FOV_UM = 300
FOV_PX = int(FOV_UM / UM_PER_PIXEL)


def load_spatialdata_contours(zarr_path: Path):
    """Load cell contours from SpatialData zarr."""
    print(f"  Loading contours from {zarr_path.name}...")
    sdata = sd.read_zarr(zarr_path)
    contours = sdata.shapes["cell_contours"]
    print(f"    Loaded {len(contours):,} contours")
    return contours


def find_dense_regions(contours_gdf, slide_dims, n_bins=25):
    """
    Find regions with highest density of each cell type.

    Returns dict of {cell_type: (center_x, center_y)}
    """
    print("  Finding dense regions for each cell type...")

    # Get centroids
    centroids = contours_gdf.geometry.centroid
    contours_gdf = contours_gdf.copy()
    contours_gdf["cx"] = centroids.x.values
    contours_gdf["cy"] = centroids.y.values

    dense_regions = {}

    for cell_type in ["Malignant epithelial", "Inflammatory", "Stroma"]:
        ct_df = contours_gdf[contours_gdf["classification_name"] == cell_type]

        if len(ct_df) < 100:
            print(f"    {cell_type}: insufficient cells ({len(ct_df)})")
            continue

        # Create 2D histogram to find densest region
        x_bins = np.linspace(0, slide_dims[0], n_bins + 1)
        y_bins = np.linspace(0, slide_dims[1], n_bins + 1)

        hist, x_edges, y_edges = np.histogram2d(
            ct_df["cx"], ct_df["cy"], bins=[x_bins, y_bins]
        )

        # Find bin with maximum count
        max_idx = np.unravel_index(hist.argmax(), hist.shape)
        center_x = (x_edges[max_idx[0]] + x_edges[max_idx[0] + 1]) / 2
        center_y = (y_edges[max_idx[1]] + y_edges[max_idx[1] + 1]) / 2

        # Make sure region is within slide bounds with margin
        margin = FOV_PX // 2 + 100
        center_x = np.clip(center_x, margin, slide_dims[0] - margin)
        center_y = np.clip(center_y, margin, slide_dims[1] - margin)

        dense_regions[cell_type] = (int(center_x), int(center_y))
        print(f"    {cell_type}: ({center_x:.0f}, {center_y:.0f}) - {hist.max():.0f} cells in bin")

    return dense_regions


def plot_filled_contours(ax, contours_gdf, x_min, x_max, y_min, y_max, alpha=0.6):
    """Plot filled cell contours in a region."""
    # Filter contours to region
    centroids = contours_gdf.geometry.centroid
    mask = (
        (centroids.x >= x_min) & (centroids.x <= x_max) &
        (centroids.y >= y_min) & (centroids.y <= y_max)
    )
    region_contours = contours_gdf[mask]

    cell_counts = {}

    for cell_type in CONSEP_ORDER:
        color = CONSEP_COLORS[cell_type]
        ct_contours = region_contours[region_contours["classification_name"] == cell_type]
        cell_counts[cell_type] = len(ct_contours)

        if len(ct_contours) == 0:
            continue

        patches = []
        for geom in ct_contours.geometry:
            if geom.geom_type == "Polygon":
                coords = np.array(geom.exterior.coords)
                patches.append(Polygon(coords, closed=True))
            elif geom.geom_type == "MultiPolygon":
                for poly in geom.geoms:
                    coords = np.array(poly.exterior.coords)
                    patches.append(Polygon(coords, closed=True))

        if patches:
            collection = PatchCollection(
                patches,
                facecolor=color,
                edgecolor=color,
                alpha=alpha,
                linewidth=0.2,
            )
            ax.add_collection(collection)

    return cell_counts


def create_figure2_hot_cold():
    """Create Figure 2 comparing immune-hot and immune-cold samples."""
    print("=" * 70)
    print("Figure 2: Immune-Hot vs Immune-Cold Comparison")
    print("=" * 70)

    # Check data availability
    for key, info in SAMPLES.items():
        sample_dir = DATA_DIR / info["sample_id"]
        tif_path = sample_dir / f"{info['sample_id']}.tif"
        zarr_path = sample_dir / f"{info['sample_id']}_consep_spatialdata.zarr"

        if not tif_path.exists():
            raise FileNotFoundError(f"H&E not found for {info['participant']}: {tif_path}")
        if not zarr_path.exists():
            raise FileNotFoundError(f"SpatialData zarr not found for {info['participant']}: {zarr_path}")

    # Load data for both samples
    sample_data = {}
    for key, info in SAMPLES.items():
        print(f"\nLoading {info['participant']} ({key})...")
        sample_dir = DATA_DIR / info["sample_id"]

        # Load H&E
        tif_path = sample_dir / f"{info['sample_id']}.tif"
        slide = openslide.OpenSlide(str(tif_path))
        print(f"  H&E dimensions: {slide.dimensions}")

        # Load contours from SpatialData zarr
        zarr_path = sample_dir / f"{info['sample_id']}_consep_spatialdata.zarr"
        contours = load_spatialdata_contours(zarr_path)

        # Find dense regions
        dense_regions = find_dense_regions(contours, slide.dimensions)

        sample_data[key] = {
            "slide": slide,
            "contours": contours,
            "dense_regions": dense_regions,
            "dims": slide.dimensions,
            **info,
        }

    # Create figure layout: 5 rows x 2 columns
    fig = plt.figure(figsize=(16, 24))
    gs = fig.add_gridspec(5, 3, width_ratios=[1, 1, 0.15], height_ratios=[1.5, 1.5, 1, 1, 1],
                          hspace=0.12, wspace=0.1)

    # Row labels
    row_labels = ["H&E Only", "H&E + Contours", "Malignant Region", "Inflammatory Region", "Stromal Region"]

    # Process each row
    for row_idx, row_name in enumerate(row_labels):
        print(f"\nProcessing row {row_idx}: {row_name}...")

        for col_idx, (key, data) in enumerate(sample_data.items()):
            ax = fig.add_subplot(gs[row_idx, col_idx])

            slide = data["slide"]
            contours = data["contours"]
            dims = data["dims"]

            if row_idx == 0:
                # Row 0: H&E only (no contours) with bounding boxes
                thumb_width = 1500
                thumb_height = int(thumb_width * dims[1] / dims[0])
                thumbnail = slide.get_thumbnail((thumb_width, thumb_height))
                thumb_arr = np.array(thumbnail.convert("RGB"))

                ax.imshow(thumb_arr, extent=[0, dims[0], dims[1], 0])

                # Draw prominent bounding boxes for zoom regions
                for ct, (cx, cy) in data["dense_regions"].items():
                    color = CONSEP_COLORS[ct]
                    half = FOV_PX // 2
                    # White outline for contrast
                    rect_outline = Rectangle(
                        (cx - half, cy - half), FOV_PX, FOV_PX,
                        fill=False, edgecolor="white", linewidth=5
                    )
                    ax.add_patch(rect_outline)
                    # Colored box
                    rect = Rectangle(
                        (cx - half, cy - half), FOV_PX, FOV_PX,
                        fill=False, edgecolor=color, linewidth=3
                    )
                    ax.add_patch(rect)

                ax.set_xlim(0, dims[0])
                ax.set_ylim(dims[1], 0)
                add_scale_bar(ax, 5000, UM_PER_PIXEL, position="lower right", color="black", fontsize=8)

                # Title
                if col_idx == 0:
                    ax.set_title(f"A. {data['label']}", fontsize=12, fontweight="bold", loc="left")
                else:
                    ax.set_title(f"B. {data['label']}", fontsize=12, fontweight="bold", loc="left")

            elif row_idx == 1:
                # Row 1: H&E with cell contours overlay
                thumb_width = 1500
                thumb_height = int(thumb_width * dims[1] / dims[0])
                thumbnail = slide.get_thumbnail((thumb_width, thumb_height))
                thumb_arr = np.array(thumbnail.convert("RGB"))

                ax.imshow(thumb_arr, extent=[0, dims[0], dims[1], 0], alpha=0.6)

                # Plot filled contours
                print(f"    Plotting overview contours for {data['participant']}...")
                cell_counts = plot_filled_contours(
                    ax, contours, 0, dims[0], 0, dims[1], alpha=0.5
                )

                # Draw prominent bounding boxes for zoom regions
                for ct, (cx, cy) in data["dense_regions"].items():
                    color = CONSEP_COLORS[ct]
                    half = FOV_PX // 2
                    # White outline for contrast
                    rect_outline = Rectangle(
                        (cx - half, cy - half), FOV_PX, FOV_PX,
                        fill=False, edgecolor="white", linewidth=5
                    )
                    ax.add_patch(rect_outline)
                    # Colored box
                    rect = Rectangle(
                        (cx - half, cy - half), FOV_PX, FOV_PX,
                        fill=False, edgecolor=color, linewidth=3
                    )
                    ax.add_patch(rect)

                ax.set_xlim(0, dims[0])
                ax.set_ylim(dims[1], 0)
                add_scale_bar(ax, 5000, UM_PER_PIXEL, position="lower right", color="black", fontsize=8)

                # Title
                if col_idx == 0:
                    ax.set_title("C. ConSeP Cell Segmentation", fontsize=12, fontweight="bold", loc="left")
                else:
                    ax.set_title("D. ConSeP Cell Segmentation", fontsize=12, fontweight="bold", loc="left")

            else:
                # Zoom rows (rows 2-4)
                cell_type_map = {2: "Malignant epithelial", 3: "Inflammatory", 4: "Stroma"}
                cell_type = cell_type_map[row_idx]

                if cell_type not in data["dense_regions"]:
                    ax.text(0.5, 0.5, f"No dense {cell_type}\nregion found",
                            ha="center", va="center", transform=ax.transAxes)
                    ax.axis("off")
                    continue

                cx, cy = data["dense_regions"][cell_type]
                half = FOV_PX // 2
                x_min, x_max = cx - half, cx + half
                y_min, y_max = cy - half, cy + half

                # Get H&E crop
                he_region = slide.read_region((x_min, y_min), 0, (FOV_PX, FOV_PX))
                he_crop = np.array(he_region.convert("RGB"))

                # Display H&E
                ax.imshow(he_crop, extent=[x_min, x_max, y_max, y_min], alpha=0.5)

                # Plot filled contours
                print(f"    Plotting {cell_type} zoom for {data['participant']}...")
                cell_counts = plot_filled_contours(
                    ax, contours, x_min, x_max, y_min, y_max, alpha=0.7
                )

                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_max, y_min)
                add_scale_bar(ax, 50, UM_PER_PIXEL, position="lower right", color="white", fontsize=7)

                # Panel label (E-J for zoom rows)
                panel_letters = {2: ("E", "F"), 3: ("G", "H"), 4: ("I", "J")}
                letter = panel_letters[row_idx][col_idx]

                # Add cell count annotation
                total = sum(cell_counts.values())
                target_count = cell_counts.get(cell_type, 0)
                target_pct = 100 * target_count / total if total > 0 else 0

                ax.annotate(
                    f"{cell_type}: {target_pct:.0f}%\n(n={total:,})",
                    xy=(0.02, 0.98), xycoords="axes fraction",
                    ha="left", va="top", fontsize=7,
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.7, edgecolor="none")
                )

                if col_idx == 0:
                    ax.set_title(f"{letter}. {cell_type} Region", fontsize=10, fontweight="bold", loc="left")
                else:
                    ax.set_title(f"{letter}.", fontsize=10, fontweight="bold", loc="left")

            ax.axis("off")

    # Legend panel
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis("off")
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)

    y_pos = 0.95
    ax_legend.text(0.1, y_pos, "Cell Types", fontsize=10, fontweight="bold", transform=ax_legend.transAxes)
    y_pos -= 0.04

    for ct in CONSEP_ORDER:
        color = CONSEP_COLORS[ct]
        ax_legend.add_patch(plt.Rectangle(
            (0.1, y_pos - 0.015), 0.15, 0.025,
            facecolor=color, edgecolor="none",
            transform=ax_legend.transAxes
        ))
        ax_legend.text(0.3, y_pos, ct, fontsize=8, va="center", transform=ax_legend.transAxes)
        y_pos -= 0.04

    # Main title
    fig.suptitle(
        "Figure 2: Cell Segmentation - Immune-Hot vs Immune-Cold",
        fontsize=14, fontweight="bold", y=0.98
    )

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure2_hot_cold_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    plt.close()

    # Close slides
    for data in sample_data.values():
        data["slide"].close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure2_hot_cold()
