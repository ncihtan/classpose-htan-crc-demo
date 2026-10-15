#!/usr/bin/env python3
"""
Figure 2A: Slide Overview with Region Annotations

Shows full H&E thumbnail with annotated regions of interest for detailed views.
Uses HTA7_954 (dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A) as representative sample.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import openslide
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import UM_PER_PIXEL, add_scale_bar, get_project_root

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
PARTICIPANT_ID = "HTA7_954"

# Input paths (SAMPLE_DATA_DIR from config.py)
HE_TIF_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}.tif"

# Field of view for region crops (in micrometers)
FOV_UM = 500
FOV_PX = int(FOV_UM / UM_PER_PIXEL)

# Regions of interest: (center_x, center_y, name, color)
# Selected to show different tissue types
REGIONS = [
    (58000, 18000, "Immune-rich", "#e41a1c"),       # Red - matches Inflammatory
    (15000, 38000, "Glandular epithelial", "#377eb8"),  # Blue - matches Healthy epithelial
    (40000, 30000, "Tumor-stroma interface", "#4daf4a"),  # Green - matches Stroma
]


def create_figure2a():
    """Create Figure 2A: Slide overview with region annotations."""
    print("=" * 60)
    print("Figure 2A: Slide Overview")
    print("=" * 60)

    # Check that H&E file exists
    if not HE_TIF_PATH.exists():
        raise FileNotFoundError(
            f"H&E tif not found: {HE_TIF_PATH}\n"
            "Please ensure the prototype sample data is available."
        )

    # Load H&E slide
    print(f"Loading H&E slide: {SAMPLE_ID}...")
    slide = openslide.OpenSlide(str(HE_TIF_PATH))
    slide_dims = slide.dimensions
    print(f"  Slide dimensions: {slide_dims[0]:,} x {slide_dims[1]:,} pixels")
    print(f"  Physical size: {slide_dims[0] * UM_PER_PIXEL / 1000:.1f} x {slide_dims[1] * UM_PER_PIXEL / 1000:.1f} mm")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Get thumbnail
    thumbnail_width = 2000
    thumbnail_height = int(thumbnail_width * slide_dims[1] / slide_dims[0])
    thumbnail = slide.get_thumbnail((thumbnail_width, thumbnail_height))
    thumbnail_arr = np.array(thumbnail.convert("RGB"))

    # Display thumbnail
    ax.imshow(thumbnail_arr, extent=[0, slide_dims[0], slide_dims[1], 0])

    # Draw region boxes
    for idx, (cx, cy, name, color) in enumerate(REGIONS):
        half_fov = FOV_PX // 2
        x_min = cx - half_fov
        y_min = cy - half_fov

        # Draw rectangle
        rect = Rectangle(
            (x_min, y_min), FOV_PX, FOV_PX,
            fill=False, edgecolor=color, linewidth=3
        )
        ax.add_patch(rect)

        # Add label
        ax.text(
            cx, y_min - 800,
            f"{chr(66 + idx)}",  # B, C, D, ...
            fontsize=14, fontweight="bold",
            color=color, ha="center", va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor=color,
                alpha=0.9
            )
        )

    # Formatting
    ax.set_xlim(0, slide_dims[0])
    ax.set_ylim(slide_dims[1], 0)
    ax.axis("off")

    # Add scale bar (5 mm)
    add_scale_bar(ax, 5000, UM_PER_PIXEL, position="lower right", color="black", fontsize=10)

    # Title with sample info
    ax.set_title(
        f"A. Representative H&E Image ({PARTICIPANT_ID})",
        fontsize=14, fontweight="bold", loc="left", pad=10
    )

    # Add sample info annotation
    n_regions = len(REGIONS)
    region_text = ", ".join([f"{chr(66+i)}: {r[2]}" for i, r in enumerate(REGIONS)])
    ax.annotate(
        f"Regions: {region_text}\nFOV: {FOV_UM} µm",
        xy=(0.99, 0.02), xycoords="axes fraction",
        ha="right", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9)
    )

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure2a_slide_overview.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_path}")

    output_pdf = output_dir / "figure2a_slide_overview.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()

    # Close slide
    slide.close()

    print("\nDone!")


if __name__ == "__main__":
    create_figure2a()
