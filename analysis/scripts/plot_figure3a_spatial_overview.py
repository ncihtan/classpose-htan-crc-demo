#!/usr/bin/env python3
"""
Figure 3A: Spatial Cell Type Distribution Map

Shows full slide with cell centroids colored by ConSeP cell type.
Reveals tissue architecture through spatial distribution of cell types.
"""

import spatialdata as sd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import openslide
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import (
    UM_PER_PIXEL,
    CONSEP_COLORS,
    CONSEP_ORDER,
    add_scale_bar,
    get_project_root,
)

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
PARTICIPANT_ID = "HTA7_954"

# Input paths (SAMPLE_DATA_DIR from config.py)
HE_TIF_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}.tif"
CONSEP_ZARR_PATH = (
    SAMPLE_DATA_DIR
    / SAMPLE_ID
    / f"{SAMPLE_ID}_consep_spatialdata.zarr"
)

# Visualization parameters
POINT_SIZE = 0.1  # Small points for dense visualization
POINT_ALPHA = 0.7


def create_figure3a():
    """Create Figure 3A: Spatial cell type distribution."""
    print("=" * 60)
    print("Figure 3A: Spatial Distribution Overview")
    print("=" * 60)

    # Check data availability
    if not HE_TIF_PATH.exists():
        raise FileNotFoundError(f"H&E tif not found: {HE_TIF_PATH}")
    if not CONSEP_ZARR_PATH.exists():
        raise FileNotFoundError(f"ConSeP zarr not found: {CONSEP_ZARR_PATH}")

    # Load data
    print("Loading H&E slide for dimensions...")
    slide = openslide.OpenSlide(str(HE_TIF_PATH))
    slide_dims = slide.dimensions
    print(f"  Slide dimensions: {slide_dims[0]:,} x {slide_dims[1]:,} px")

    # Get thumbnail for background
    thumbnail_width = 2000
    thumbnail_height = int(thumbnail_width * slide_dims[1] / slide_dims[0])
    thumbnail = slide.get_thumbnail((thumbnail_width, thumbnail_height))
    thumbnail_arr = np.array(thumbnail.convert("RGB"))

    print("Loading ConSeP SpatialData...")
    sdata = sd.read_zarr(CONSEP_ZARR_PATH)

    # Get cell centroids
    print("Loading cell centroids...")
    points = sdata.points["cell_centroids"].compute()
    print(f"  Total cells: {len(points):,}")

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Display H&E thumbnail as faint background
    ax.imshow(thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=0.4)

    # Plot centroids by cell type (in order from bottom to top layer)
    # Plot in reverse order so inflammatory (most interesting) is on top
    print("Plotting cell centroids by type...")
    cell_counts = {}

    for cell_type in reversed(CONSEP_ORDER):
        color = CONSEP_COLORS[cell_type]
        mask = points["classification_name"] == cell_type
        subset = points[mask]
        cell_counts[cell_type] = len(subset)

        if len(subset) > 0:
            ax.scatter(
                subset["x"],
                subset["y"],
                c=color,
                s=POINT_SIZE,
                alpha=POINT_ALPHA,
                rasterized=True,  # Important for large point counts
                label=cell_type,
            )
            print(f"  {cell_type}: {len(subset):,} cells")

    # Set limits
    ax.set_xlim(0, slide_dims[0])
    ax.set_ylim(slide_dims[1], 0)  # Inverted for image coordinates
    ax.axis("off")

    # Add scale bar
    add_scale_bar(
        ax, 5000, UM_PER_PIXEL, position="lower right", color="black", fontsize=10
    )

    # Create legend
    legend_elements = [
        Patch(facecolor=CONSEP_COLORS[ct], label=f"{ct} ({cell_counts[ct]:,})")
        for ct in CONSEP_ORDER
        if cell_counts.get(ct, 0) > 0
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=9,
        title="Cell Type",
        title_fontsize=10,
        framealpha=0.9,
    )

    # Title
    ax.set_title(
        f"A. Spatial Cell Type Distribution ({PARTICIPANT_ID})",
        fontsize=14,
        fontweight="bold",
        loc="left",
        pad=10,
    )

    # Add annotation with total
    total_cells = sum(cell_counts.values())
    ax.annotate(
        f"Total: {total_cells:,} cells\nConSeP classification",
        xy=(0.02, 0.02),
        xycoords="axes fraction",
        ha="left",
        va="bottom",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9
        ),
    )

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure3a_spatial_overview.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "figure3a_spatial_overview.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()
    slide.close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure3a()
