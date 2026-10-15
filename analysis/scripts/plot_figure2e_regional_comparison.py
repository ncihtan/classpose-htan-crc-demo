#!/usr/bin/env python3
"""
Figure 2E: Regional Cell Composition Comparison

Shows how cell type composition varies between different tissue regions.
Compares the immune-rich and epithelial regions from panels B and C.
"""

import spatialdata as sd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from utils.spatial_utils import (
    UM_PER_PIXEL, CONSEP_COLORS, CONSEP_ORDER, get_project_root
)

# =============================================================================
# Configuration
# =============================================================================
SAMPLE_ID = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
PARTICIPANT_ID = "HTA7_954"

# Input paths (SAMPLE_DATA_DIR from config.py)
CONSEP_ZARR_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}_consep_spatialdata.zarr"

# Regions to compare (matching panels B and C)
FOV_UM = 500
FOV_PX = int(FOV_UM / UM_PER_PIXEL)

REGIONS = {
    "Immune-rich": (58000, 18000),
    "Glandular": (15000, 38000),
    "Tumor-stroma": (40000, 30000),
}


def get_region_composition(contours, region_name: str, center: tuple, fov_px: int) -> dict:
    """Get cell type composition for a region."""
    cx, cy = center
    half = fov_px // 2
    x_min, x_max = cx - half, cx + half
    y_min, y_max = cy - half, cy + half

    # Filter contours to region
    centroids = contours.geometry.centroid
    mask = (
        (centroids.x >= x_min) & (centroids.x <= x_max) &
        (centroids.y >= y_min) & (centroids.y <= y_max)
    )
    region_contours = contours[mask]

    # Count by cell type
    counts = {}
    for cell_type in CONSEP_ORDER:
        ct_count = len(region_contours[region_contours["classification_name"] == cell_type])
        counts[cell_type] = ct_count

    return counts


def create_figure2e():
    """Create Figure 2E: Regional composition comparison."""
    print("=" * 60)
    print("Figure 2E: Regional Cell Composition")
    print("=" * 60)

    # Check data availability
    if not CONSEP_ZARR_PATH.exists():
        raise FileNotFoundError(f"ConSeP zarr not found: {CONSEP_ZARR_PATH}")

    # Load data
    print("Loading ConSeP SpatialData...")
    sdata = sd.read_zarr(CONSEP_ZARR_PATH)
    contours = sdata.shapes["cell_contours"]
    print(f"  Total contours: {len(contours):,}")

    # Get composition for each region
    print("\nCalculating regional compositions...")
    region_data = {}
    for region_name, center in REGIONS.items():
        counts = get_region_composition(contours, region_name, center, FOV_PX)
        total = sum(counts.values())
        print(f"  {region_name}: {total:,} cells")
        region_data[region_name] = {
            "counts": counts,
            "total": total,
            "percentages": {ct: 100 * counts[ct] / total if total > 0 else 0 for ct in counts}
        }

    # Also get whole-slide composition
    print("\nCalculating whole-slide composition...")
    slide_counts = {}
    for cell_type in CONSEP_ORDER:
        slide_counts[cell_type] = len(contours[contours["classification_name"] == cell_type])
    slide_total = sum(slide_counts.values())
    slide_pcts = {ct: 100 * slide_counts[ct] / slide_total if slide_total > 0 else 0 for ct in slide_counts}
    print(f"  Whole slide: {slide_total:,} cells")

    # Create figure with grouped bar chart
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2, 1]})

    # Panel 1: Grouped bar chart comparing regions
    ax1 = axes[0]
    x = np.arange(len(CONSEP_ORDER))
    width = 0.22
    n_regions = len(REGIONS)

    bars_list = []
    for idx, (region_name, data) in enumerate(region_data.items()):
        offset = (idx - (n_regions - 1) / 2) * width
        pcts = [data["percentages"].get(ct, 0) for ct in CONSEP_ORDER]
        bars = ax1.bar(
            x + offset, pcts, width,
            label=region_name,
            alpha=0.8,
            edgecolor="white",
            linewidth=0.5
        )
        bars_list.append(bars)

    ax1.set_xlabel("Cell Type", fontsize=11)
    ax1.set_ylabel("Proportion (%)", fontsize=11)
    ax1.set_title("Regional Cell Type Distribution", fontsize=12, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels([ct.replace(" ", "\n") for ct in CONSEP_ORDER], fontsize=9)
    ax1.legend(loc="upper right", fontsize=9, title="Region")
    ax1.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax1.set_axisbelow(True)
    ax1.set_ylim(0, max(
        max(region_data[r]["percentages"].values()) for r in region_data
    ) * 1.15)

    # Panel 2: Whole-slide pie chart
    ax2 = axes[1]

    # Filter to cell types with >1%
    pie_data = [(ct, slide_pcts[ct]) for ct in CONSEP_ORDER if slide_pcts[ct] > 1]
    other_pct = sum(slide_pcts[ct] for ct in CONSEP_ORDER if slide_pcts[ct] <= 1)
    if other_pct > 0:
        pie_data.append(("Other (<1%)", other_pct))

    labels = [ct for ct, _ in pie_data]
    sizes = [pct for _, pct in pie_data]
    colors = [CONSEP_COLORS.get(ct, "#999999") for ct, _ in pie_data]

    wedges, texts, autotexts = ax2.pie(
        sizes, labels=None, colors=colors,
        autopct=lambda p: f"{p:.0f}%" if p > 5 else "",
        startangle=90, counterclock=False,
        wedgeprops=dict(edgecolor="white", linewidth=1)
    )

    # Add legend instead of labels (cleaner)
    ax2.legend(
        wedges, [f"{label} ({pct:.1f}%)" for label, pct in pie_data],
        loc="center left", bbox_to_anchor=(0.9, 0.5),
        fontsize=8
    )

    ax2.set_title(f"Whole Slide Composition\n(n = {slide_total:,})", fontsize=12, fontweight="bold")

    # Main title
    fig.suptitle(
        f"E. Cell Type Composition Analysis ({PARTICIPANT_ID})",
        fontsize=14, fontweight="bold", y=1.02
    )

    plt.tight_layout()

    # Save
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "figure2e_regional_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "figure2e_regional_comparison.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()

    # Print summary table
    print("\n" + "=" * 60)
    print("Regional Composition Summary")
    print("=" * 60)
    print(f"{'Cell Type':<22} {'Immune-rich':>12} {'Glandular':>12} {'Tumor-stroma':>12} {'Whole Slide':>12}")
    print("-" * 70)
    for ct in CONSEP_ORDER:
        r1 = region_data["Immune-rich"]["percentages"].get(ct, 0)
        r2 = region_data["Glandular"]["percentages"].get(ct, 0)
        r3 = region_data["Tumor-stroma"]["percentages"].get(ct, 0)
        ws = slide_pcts.get(ct, 0)
        print(f"{ct:<22} {r1:>11.1f}% {r2:>11.1f}% {r3:>11.1f}% {ws:>11.1f}%")

    print("\nDone!")


if __name__ == "__main__":
    create_figure2e()
