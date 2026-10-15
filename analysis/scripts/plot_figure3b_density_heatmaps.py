#!/usr/bin/env python3
"""
Figure 3B: Cell Density Heatmaps

Shows kernel density estimation for each cell type, revealing where
different cell populations concentrate in the tissue.
"""

import spatialdata as sd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy import ndimage
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

# Density estimation parameters
GRID_SIZE = 500  # Number of bins for density grid
SIGMA = 10  # Gaussian smoothing sigma (in grid units)


def compute_density_grid(x_coords, y_coords, x_range, y_range, grid_size, sigma):
    """
    Compute kernel density estimation on a grid.

    Parameters
    ----------
    x_coords, y_coords : array-like
        Coordinates of points
    x_range, y_range : tuple
        (min, max) ranges for x and y
    grid_size : int
        Number of bins in each dimension
    sigma : float
        Gaussian smoothing sigma

    Returns
    -------
    density : ndarray
        2D density grid
    extent : tuple
        (x_min, x_max, y_min, y_max) for plotting
    """
    # Create 2D histogram
    x_bins = np.linspace(x_range[0], x_range[1], grid_size + 1)
    y_bins = np.linspace(y_range[0], y_range[1], grid_size + 1)

    hist, _, _ = np.histogram2d(x_coords, y_coords, bins=[x_bins, y_bins])

    # Apply Gaussian smoothing
    density = ndimage.gaussian_filter(hist.T, sigma=sigma)

    extent = (x_range[0], x_range[1], y_range[1], y_range[0])  # Note: y inverted

    return density, extent


def create_colormap_from_color(color, name="custom"):
    """Create a white-to-color colormap."""
    # Parse hex color
    r = int(color[1:3], 16) / 255
    g = int(color[3:5], 16) / 255
    b = int(color[5:7], 16) / 255

    colors = [(1, 1, 1, 0), (r, g, b, 1)]  # White transparent to color opaque
    return LinearSegmentedColormap.from_list(name, colors)


def create_figure3b():
    """Create Figure 3B: Cell density heatmaps."""
    print("=" * 60)
    print("Figure 3B: Cell Density Heatmaps")
    print("=" * 60)

    # Check data availability
    if not HE_TIF_PATH.exists():
        raise FileNotFoundError(f"H&E tif not found: {HE_TIF_PATH}")
    if not CONSEP_ZARR_PATH.exists():
        raise FileNotFoundError(f"ConSeP zarr not found: {CONSEP_ZARR_PATH}")

    # Load data
    print("Loading H&E slide for background...")
    slide = openslide.OpenSlide(str(HE_TIF_PATH))
    slide_dims = slide.dimensions
    print(f"  Slide dimensions: {slide_dims[0]:,} x {slide_dims[1]:,} px")

    # Get thumbnail
    thumbnail_width = 1000
    thumbnail_height = int(thumbnail_width * slide_dims[1] / slide_dims[0])
    thumbnail = slide.get_thumbnail((thumbnail_width, thumbnail_height))
    thumbnail_arr = np.array(thumbnail.convert("RGB"))

    print("Loading ConSeP SpatialData...")
    sdata = sd.read_zarr(CONSEP_ZARR_PATH)
    points = sdata.points["cell_centroids"].compute()
    print(f"  Total cells: {len(points):,}")

    # Define coordinate ranges
    x_range = (0, slide_dims[0])
    y_range = (0, slide_dims[1])

    # Create figure: 2x3 grid of cell types
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    print("\nComputing density maps...")

    for idx, cell_type in enumerate(CONSEP_ORDER):
        ax = axes[idx]
        color = CONSEP_COLORS[cell_type]

        # Get points for this cell type
        mask = points["classification_name"] == cell_type
        subset = points[mask]
        n_cells = len(subset)
        print(f"  {cell_type}: {n_cells:,} cells")

        # Display H&E background
        ax.imshow(thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0], alpha=0.5)

        if n_cells > 100:  # Only compute density if enough points
            # Compute density
            density, extent = compute_density_grid(
                subset["x"].values, subset["y"].values,
                x_range, y_range, GRID_SIZE, SIGMA
            )

            # Normalize density
            if density.max() > 0:
                density = density / density.max()

            # Create colormap
            cmap = create_colormap_from_color(color)

            # Overlay density
            ax.imshow(
                density, extent=extent, origin="upper",
                cmap=cmap, alpha=0.7, vmin=0, vmax=1
            )

        # Formatting
        ax.set_xlim(0, slide_dims[0])
        ax.set_ylim(slide_dims[1], 0)
        ax.axis("off")

        # Title with count
        ax.set_title(
            f"{cell_type}\n(n = {n_cells:,})",
            fontsize=10, fontweight="bold", color=color
        )

    # Add scale bar to last panel
    add_scale_bar(axes[-1], 5000, UM_PER_PIXEL, position="lower right", color="black", fontsize=8)

    # Main title
    fig.suptitle(
        f"B. Cell Density Maps by Type ({PARTICIPANT_ID})",
        fontsize=14, fontweight="bold", y=0.98
    )

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure3b_density_heatmaps.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "figure3b_density_heatmaps.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()
    slide.close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure3b()
