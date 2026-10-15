#!/usr/bin/env python3
"""
Figure 2C: Epithelial/Glandular Region with Cell Contours

Shows zoomed view of glandular epithelial region with ConSeP cell contours.
Highlights distinction between malignant and healthy epithelial cells.
"""

import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
import openslide
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import (
    UM_PER_PIXEL, CONSEP_COLORS, CONSEP_ORDER,
    add_scale_bar, plot_contours, create_legend, get_project_root
)

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
PARTICIPANT_ID = "HTA7_954"

# Input paths (SAMPLE_DATA_DIR from config.py)
HE_TIF_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}.tif"
CONSEP_ZARR_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}_consep_spatialdata.zarr"

# Region of interest: Glandular epithelial region
REGION_CENTER = (15000, 38000)
REGION_NAME = "Glandular Epithelial"
FOV_UM = 500
FOV_PX = int(FOV_UM / UM_PER_PIXEL)


def create_figure2c():
    """Create Figure 2C: Epithelial region with cell contours."""
    print("=" * 60)
    print("Figure 2C: Epithelial Region")
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

    # Calculate region bounds
    cx, cy = REGION_CENTER
    half_fov = FOV_PX // 2
    x_min, x_max = cx - half_fov, cx + half_fov
    y_min, y_max = cy - half_fov, cy + half_fov

    print(f"\nRegion: {REGION_NAME}")
    print(f"  Center: ({cx}, {cy})")
    print(f"  Bounds: x=[{x_min}, {x_max}], y=[{y_min}, {y_max}]")
    print(f"  FOV: {FOV_UM} µm ({FOV_PX} px)")

    # Read H&E crop
    he_region = slide.read_region((x_min, y_min), 0, (FOV_PX, FOV_PX))
    he_crop = np.array(he_region.convert("RGB"))

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Display H&E background
    ax.imshow(he_crop, extent=[x_min, x_max, y_max, y_min])

    # Plot cell contours
    print("\nPlotting cell contours...")
    cell_counts, total_cells = plot_contours(
        ax, contours, CONSEP_COLORS,
        x_min, x_max, y_min, y_max,
        linewidth=0.5, alpha=0.9
    )

    # Print cell type breakdown
    print(f"  Total cells in region: {total_cells:,}")
    for ct in CONSEP_ORDER:
        count = cell_counts.get(ct, 0)
        pct = 100 * count / total_cells if total_cells > 0 else 0
        if count > 0:
            print(f"    {ct}: {count:,} ({pct:.1f}%)")

    # Calculate epithelial ratio for annotation
    malignant = cell_counts.get("Malignant epithelial", 0)
    healthy = cell_counts.get("Healthy epithelial", 0)
    total_epithelial = malignant + healthy
    if total_epithelial > 0:
        malignant_pct = 100 * malignant / total_epithelial
        print(f"\n  Epithelial breakdown: {malignant_pct:.1f}% malignant, {100-malignant_pct:.1f}% healthy")

    # Set limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_max, y_min)  # Inverted for image coordinates
    ax.axis("off")

    # Add scale bar
    add_scale_bar(ax, 100, UM_PER_PIXEL, position="lower right", color="white", fontsize=9)

    # Add legend
    create_legend(
        ax, CONSEP_COLORS, cell_counts,
        title=f"n = {total_cells:,}",
        fontsize=8, loc="upper right"
    )

    # Title
    ax.set_title(
        f"C. {REGION_NAME} Region",
        fontsize=14, fontweight="bold", loc="left", pad=10
    )

    # Add region info
    ax.annotate(
        f"ConSeP model\n{PARTICIPANT_ID}",
        xy=(0.02, 0.98), xycoords="axes fraction",
        ha="left", va="top", fontsize=8,
        color="white",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.7, edgecolor="none")
    )

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure2c_region_epithelial.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "figure2c_region_epithelial.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()
    slide.close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure2c()
