#!/usr/bin/env python3
"""
Plot region assignments overlaid on greyscale H&E images.

Creates visualizations showing tumor core, margin, stroma, and other regions
as a semi-transparent categorical heatmap over the H&E image.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation
import openslide

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root


# Region colors (with good contrast against greyscale)
REGION_COLORS = {
    0: "#00000000",  # Other - transparent
    1: "#d62728",    # Tumor core - red
    2: "#ff7f0e",    # Tumor margin - orange
    3: "#2ca02c",    # Stroma - green
}

REGION_LABELS = {
    0: "Other",
    1: "Tumor core",
    2: "Tumor margin",
    3: "Stroma",
}

# ConSeP cell types
CONSEP_TYPES = [
    "Healthy epithelial",
    "Malignant epithelial",
    "Stroma",
    "Inflammatory",
    "Muscle",
    "Other",
]

# Representative samples
PHENOTYPE_SAMPLES = {
    "Hot": "dg_4DFC_B0A48536-0171-4411-A200-88344086F3D4",
    "Altered-immunosuppressed": "dg_4DFC_96A5875F-826A-4917-A7DE-18E489024DF9",
    "Altered-excluded": "dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF",
    "Immune-desert": "dg_4DFC_8C77314D-9FAC-4917-83CE-52A3B70FE014",
}


def load_sample_points(sample_id: str) -> pd.DataFrame:
    """Load cell centroids from parquet."""
    parquet_path = (
        SAMPLE_DATA_DIR / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
        / "points" / "cell_centroids" / "points.parquet"
    )
    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


def load_he_thumbnail(sample_id: str, max_size: int = 2000) -> tuple[np.ndarray, float]:
    """
    Load H&E image as greyscale thumbnail.

    Returns (image_array, downsample_factor)
    """
    tif_path = SAMPLE_DATA_DIR / sample_id / f"{sample_id}.tif"

    if not tif_path.exists():
        return None, 1.0

    slide = openslide.OpenSlide(str(tif_path))

    # Calculate downsample to fit in max_size
    w, h = slide.dimensions
    downsample = max(w, h) / max_size

    # Get thumbnail
    thumb_size = (int(w / downsample), int(h / downsample))
    thumbnail = slide.get_thumbnail(thumb_size)

    # Convert to greyscale
    img_array = np.array(thumbnail.convert("L"))

    slide.close()
    return img_array, downsample


def compute_region_grid(
    points_df: pd.DataFrame,
    n_bins: int = 50,
    min_cells_per_bin: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute region assignments for each grid bin.

    Returns (region_grid, x_bins, y_bins)
    Region codes: 0=other, 1=tumor_core, 2=tumor_margin, 3=stroma
    """
    # Get coordinate ranges
    x_min, x_max = points_df["x"].min(), points_df["x"].max()
    y_min, y_max = points_df["y"].min(), points_df["y"].max()

    x_bins = np.linspace(x_min, x_max, n_bins + 1)
    y_bins = np.linspace(y_min, y_max, n_bins + 1)

    # Count cells by type per bin
    histograms = {}
    for cell_type in CONSEP_TYPES:
        mask = points_df["classification_name"] == cell_type
        hist, _, _ = np.histogram2d(
            points_df.loc[mask, "x"].values,
            points_df.loc[mask, "y"].values,
            bins=[x_bins, y_bins],
        )
        histograms[cell_type] = hist

    # Compute totals and percentages
    total_per_bin = sum(histograms.values())
    total_safe = total_per_bin + 1e-6

    malignant_pct = histograms["Malignant epithelial"] / total_safe
    stroma_pct = histograms["Stroma"] / total_safe

    # Initialize region grid
    region_grid = np.zeros((n_bins, n_bins), dtype=int)

    # Tumor core: top 25% malignant density among bins with malignant cells
    malignant_bins = malignant_pct > 0.1
    if malignant_bins.sum() > 0:
        threshold = np.percentile(malignant_pct[malignant_bins], 75)
        core_mask = malignant_pct >= threshold
        region_grid[core_mask] = 1  # tumor core

        # Tumor margin: adjacent to core, moderate malignant
        core_dilated = binary_dilation(core_mask, iterations=1)
        margin_mask = core_dilated & ~core_mask & (malignant_pct > 0.05)
        region_grid[margin_mask] = 2  # tumor margin

    # Stroma: low malignant, high stroma, sufficient cells
    stroma_mask = (
        (malignant_pct < 0.10) &
        (stroma_pct > 0.30) &
        (total_per_bin > min_cells_per_bin) &
        (region_grid == 0)
    )
    region_grid[stroma_mask] = 3  # stroma

    return region_grid, x_bins, y_bins


def plot_region_overlay(
    sample_id: str,
    phenotype: str,
    ax: plt.Axes,
    n_bins: int = 50,
    alpha: float = 0.5,
):
    """Plot region assignments overlaid on greyscale H&E or cell scatter."""

    # Load cell points
    points_df = load_sample_points(sample_id)

    # Compute region grid
    region_grid, x_bins, y_bins = compute_region_grid(points_df, n_bins=n_bins)

    # Load H&E thumbnail
    he_img, downsample = load_he_thumbnail(sample_id)

    # Get coordinate extents
    x_extent = (x_bins[0], x_bins[-1])
    y_extent = (y_bins[0], y_bins[-1])

    if he_img is not None:
        # Plot greyscale H&E
        ax.imshow(
            he_img,
            cmap="gray",
            extent=[0, he_img.shape[1] * downsample, he_img.shape[0] * downsample, 0],
            aspect="equal",
        )
    else:
        # No H&E available - plot cells as grey background
        ax.scatter(
            points_df["x"],
            points_df["y"],
            c="#cccccc",
            s=0.05,
            alpha=0.5,
            linewidths=0,
            rasterized=True,
        )
        ax.set_facecolor("white")

    # Create categorical colormap
    cmap_colors = ["#ffffff00", "#d62728", "#ff7f0e", "#2ca02c"]  # transparent, red, orange, green
    cmap = ListedColormap(cmap_colors)

    # Plot region overlay
    # Need to mask out "other" (0) regions for transparency
    region_masked = np.ma.masked_where(region_grid == 0, region_grid)

    ax.imshow(
        region_masked.T,  # Transpose to match image orientation
        cmap=cmap,
        extent=[x_extent[0], x_extent[1], y_extent[1], y_extent[0]],
        alpha=alpha,
        aspect="equal",
        interpolation="nearest",
        vmin=0,
        vmax=3,
    )

    # Set axis limits to data extent
    ax.set_xlim(x_extent)
    ax.set_ylim(y_extent[1], y_extent[0])  # Inverted for image coords

    # Get sample metadata
    classification = pd.read_csv(get_project_root() / "data" / "immune_phenotype_classification.csv")
    sample_info = classification[classification["sample_id"] == sample_id].iloc[0]

    # Title
    title = f"{phenotype}\n{sample_info['participant_id']}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])

    # Count regions
    n_core = (region_grid == 1).sum()
    n_margin = (region_grid == 2).sum()
    n_stroma = (region_grid == 3).sum()

    return {"core": n_core, "margin": n_margin, "stroma": n_stroma}


def create_region_comparison():
    """Create 2x2 figure showing region assignments for each phenotype."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    axes = axes.flatten()

    phenotype_order = ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]

    for ax, phenotype in zip(axes, phenotype_order):
        sample_id = PHENOTYPE_SAMPLES[phenotype]
        print(f"Plotting {phenotype} ({sample_id})...")
        plot_region_overlay(sample_id, phenotype, ax, n_bins=75, alpha=0.5)

    # Add legend
    legend_elements = [
        Patch(facecolor="#d62728", alpha=0.7, label="Tumor core"),
        Patch(facecolor="#ff7f0e", alpha=0.7, label="Tumor margin"),
        Patch(facecolor="#2ca02c", alpha=0.7, label="Stroma"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        fontsize=11,
        frameon=True,
        bbox_to_anchor=(0.5, 0.02),
    )

    plt.suptitle(
        "Tumor Microenvironment Region Assignment",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])

    output_path = get_project_root() / "figures" / "region_assignment_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nSaved to: {output_path}")


def create_individual_region_panels():
    """Create individual high-res region overlay panels."""

    output_dir = get_project_root() / "figures"

    for phenotype, sample_id in PHENOTYPE_SAMPLES.items():
        print(f"Plotting {phenotype}...")

        fig, ax = plt.subplots(figsize=(12, 12))
        counts = plot_region_overlay(sample_id, phenotype, ax, n_bins=100, alpha=0.5)

        # Add legend
        legend_elements = [
            Patch(facecolor="#d62728", alpha=0.7, label=f"Tumor core ({counts['core']} bins)"),
            Patch(facecolor="#ff7f0e", alpha=0.7, label=f"Tumor margin ({counts['margin']} bins)"),
            Patch(facecolor="#2ca02c", alpha=0.7, label=f"Stroma ({counts['stroma']} bins)"),
        ]
        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=10,
            frameon=True,
            facecolor="white",
            edgecolor="gray",
        )

        safe_name = phenotype.lower().replace("-", "_")
        output_path = output_dir / f"region_overlay_{safe_name}.png"
        plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"  Saved: {output_path}")


if __name__ == "__main__":
    print("Creating region assignment overlays...")
    print("=" * 50)

    # Create comparison figure
    create_region_comparison()

    # Create individual panels
    print("\nCreating individual panels...")
    create_individual_region_panels()

    print("\nDone!")
