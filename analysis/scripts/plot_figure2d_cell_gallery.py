#!/usr/bin/env python3
"""
Figure 2D: Cell-Level Classification Examples

Gallery of individual cells showing each ConSeP cell type.
Shows what each classification looks like at high magnification.
"""

import spatialdata as sd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import openslide
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import (
    UM_PER_PIXEL, CONSEP_COLORS, CONSEP_ORDER, add_scale_bar, get_project_root
)

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
PARTICIPANT_ID = "HTA7_954"

# Input paths (SAMPLE_DATA_DIR from config.py)
HE_TIF_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}.tif"
CONSEP_ZARR_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}_consep_spatialdata.zarr"

# Cell gallery configuration
CELLS_PER_TYPE = 4
CELL_CROP_UM = 40  # Size of each cell crop in micrometers
CELL_CROP_PX = int(CELL_CROP_UM / UM_PER_PIXEL)


def find_representative_cells(contours, cell_type: str, n_cells: int = 4):
    """
    Find representative cells of a given type.

    Selects cells from different regions and of typical sizes.
    """
    ct_contours = contours[contours["classification_name"] == cell_type].copy()

    if len(ct_contours) == 0:
        return []

    # Get centroids
    centroids = ct_contours.geometry.centroid
    ct_contours["cx"] = centroids.x.values
    ct_contours["cy"] = centroids.y.values

    # Calculate cell areas for filtering
    ct_contours["area"] = ct_contours.geometry.area

    # Filter to typical sizes (middle 80% of area distribution)
    q10 = ct_contours["area"].quantile(0.1)
    q90 = ct_contours["area"].quantile(0.9)
    typical_cells = ct_contours[
        (ct_contours["area"] >= q10) & (ct_contours["area"] <= q90)
    ]

    if len(typical_cells) < n_cells:
        typical_cells = ct_contours

    # Sample cells spread across the slide
    # Divide into spatial bins and sample from each
    if len(typical_cells) >= n_cells:
        x_bins = np.linspace(typical_cells["cx"].min(), typical_cells["cx"].max(), n_cells + 1)
        selected = []

        for i in range(n_cells):
            bin_cells = typical_cells[
                (typical_cells["cx"] >= x_bins[i]) &
                (typical_cells["cx"] < x_bins[i + 1])
            ]
            if len(bin_cells) > 0:
                # Select cell closest to median area in this bin
                median_area = bin_cells["area"].median()
                best_idx = (bin_cells["area"] - median_area).abs().idxmin()
                selected.append(bin_cells.loc[best_idx])

        # If we didn't get enough, sample randomly from remaining
        if len(selected) < n_cells:
            remaining = typical_cells[~typical_cells.index.isin([s.name for s in selected])]
            n_needed = n_cells - len(selected)
            if len(remaining) >= n_needed:
                additional = remaining.sample(n_needed)
                selected.extend([row for _, row in additional.iterrows()])

        return selected[:n_cells]

    return [row for _, row in typical_cells.head(n_cells).iterrows()]


def plot_cell_with_contour(ax, slide, cell_row, color, crop_px):
    """Plot a single cell with its contour overlay."""
    cx, cy = cell_row["cx"], cell_row["cy"]
    geom = cell_row.geometry

    # Calculate crop bounds
    half = crop_px // 2
    x_min = int(cx - half)
    y_min = int(cy - half)

    # Read H&E crop
    try:
        he_region = slide.read_region((x_min, y_min), 0, (crop_px, crop_px))
        he_crop = np.array(he_region.convert("RGB"))
    except Exception:
        # If region is outside slide bounds, return False
        return False

    # Display H&E
    ax.imshow(he_crop, extent=[x_min, x_min + crop_px, y_min + crop_px, y_min])

    # Draw contour
    patches = []
    if geom.geom_type == "Polygon":
        coords = np.array(geom.exterior.coords)
        patches.append(Polygon(coords, closed=True))
    elif geom.geom_type == "MultiPolygon":
        for poly in geom.geoms:
            coords = np.array(poly.exterior.coords)
            patches.append(Polygon(coords, closed=True))

    if patches:
        collection = PatchCollection(
            patches, facecolor="none", edgecolor=color,
            linewidth=1.5, alpha=1.0
        )
        ax.add_collection(collection)

    # Set limits
    ax.set_xlim(x_min, x_min + crop_px)
    ax.set_ylim(y_min + crop_px, y_min)
    ax.axis("off")

    return True


def create_figure2d():
    """Create Figure 2D: Cell gallery."""
    print("=" * 60)
    print("Figure 2D: Cell Gallery")
    print("=" * 60)

    # Check data availability
    if not HE_TIF_PATH.exists():
        raise FileNotFoundError(f"H&E tif not found: {HE_TIF_PATH}")
    if not CONSEP_ZARR_PATH.exists():
        raise FileNotFoundError(f"ConSeP zarr not found: {CONSEP_ZARR_PATH}")

    # Load data
    print("Loading H&E slide...")
    slide = openslide.OpenSlide(str(HE_TIF_PATH))

    print("Loading ConSeP SpatialData...")
    sdata = sd.read_zarr(CONSEP_ZARR_PATH)
    contours = sdata.shapes["cell_contours"]
    print(f"  Total contours: {len(contours):,}")

    # Find representative cells for each type
    print("\nSelecting representative cells...")
    selected_cells = {}
    for cell_type in CONSEP_ORDER:
        cells = find_representative_cells(contours, cell_type, CELLS_PER_TYPE)
        selected_cells[cell_type] = cells
        print(f"  {cell_type}: {len(cells)} cells selected")

    # Create figure: grid of cell types (rows) x examples (columns)
    n_types = len(CONSEP_ORDER)
    n_cols = CELLS_PER_TYPE + 1  # Extra column for labels

    fig, axes = plt.subplots(
        n_types, n_cols,
        figsize=(2 * n_cols, 2 * n_types),
        gridspec_kw={"width_ratios": [0.5] + [1] * CELLS_PER_TYPE}
    )

    # Plot cells
    for row_idx, cell_type in enumerate(CONSEP_ORDER):
        color = CONSEP_COLORS[cell_type]
        cells = selected_cells[cell_type]

        # Label column
        ax_label = axes[row_idx, 0]
        ax_label.axis("off")
        ax_label.text(
            0.9, 0.5, cell_type,
            ha="right", va="center",
            fontsize=9, fontweight="bold",
            color=color,
            transform=ax_label.transAxes
        )

        # Cell columns
        for col_idx, cell in enumerate(cells):
            ax = axes[row_idx, col_idx + 1]
            success = plot_cell_with_contour(ax, slide, cell, color, CELL_CROP_PX)
            if not success:
                ax.axis("off")
                ax.text(0.5, 0.5, "N/A", ha="center", va="center", transform=ax.transAxes)

        # Fill remaining columns if fewer cells
        for col_idx in range(len(cells), CELLS_PER_TYPE):
            ax = axes[row_idx, col_idx + 1]
            ax.axis("off")

    # Add scale bar to bottom-right cell
    ax_scale = axes[-1, -1]
    add_scale_bar(ax_scale, 10, UM_PER_PIXEL, position="lower right", color="white", fontsize=7)

    # Title
    fig.suptitle(
        "D. Cell-Level Classification Examples",
        fontsize=14, fontweight="bold", y=0.98
    )

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure2d_cell_gallery.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "figure2d_cell_gallery.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()
    slide.close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure2d()
