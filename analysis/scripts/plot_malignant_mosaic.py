#!/usr/bin/env python3
"""
Mosaic grid of all 41 samples showing all cell types.

Each subplot shows the spatial distribution of all ConSeP cell types
for one sample, arranged in a grid sorted by immune phenotype.
Uses spatialdata to load data, matplotlib scatter for accurate rendering.
Points are kept extremely small to minimise overlap at this scale.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import spatialdata as sd

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root
from utils.spatial_utils import CONSEP_COLORS, CONSEP_ORDER

# =============================================================================
# Configuration
# =============================================================================
PROJECT_ROOT = get_project_root()

# Plot order: background types first, malignant + inflammatory on top
PLOT_ORDER = ["Other", "Muscle", "Healthy epithelial", "Stroma", "Malignant epithelial", "Inflammatory"]

# Phenotype ordering and colors for border/label
PHENOTYPE_ORDER = ["Hot", "Altered-excluded", "Altered-immunosuppressed", "Immune-desert"]
PHENOTYPE_COLORS = {
    "Hot": "#d62728",
    "Altered-excluded": "#ff7f0e",
    "Altered-immunosuppressed": "#2ca02c",
    "Immune-desert": "#1f77b4",
}


def load_sample_points(sample_id: str) -> pd.DataFrame:
    """Load cell centroids via spatialdata."""
    zarr_path = (
        SAMPLE_DATA_DIR
        / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
    )
    sdata = sd.read_zarr(zarr_path)
    points = sdata.points["cell_centroids"].compute()
    return points[["x", "y", "classification_name"]]


def create_mosaic():
    """Create a mosaic grid showing all cell types for all samples."""
    print("=" * 60)
    print("Cell Type Mosaic - All 41 Samples")
    print("=" * 60)

    # Load phenotype classification for ordering and labels
    classification = pd.read_csv(
        PROJECT_ROOT / "data" / "immune_phenotype_classification.csv"
    )

    # Sort by phenotype group, then by inflammatory % within group
    classification["phenotype_order"] = classification["phenotype"].map(
        {p: i for i, p in enumerate(PHENOTYPE_ORDER)}
    )
    classification = classification.sort_values(
        ["phenotype_order", "inflammatory_pct_overall"],
        ascending=[True, False],
    ).reset_index(drop=True)

    n_samples = len(classification)
    n_cols = 7
    n_rows = int(np.ceil(n_samples / n_cols))

    print(f"Grid: {n_rows} rows x {n_cols} cols for {n_samples} samples")

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 3, n_rows * 3),
        constrained_layout=True,
    )
    axes = axes.flatten()

    for idx, (_, row) in enumerate(classification.iterrows()):
        ax = axes[idx]
        sample_id = row["sample_id"]
        participant_id = row["participant_id"]
        phenotype = row["phenotype"]

        print(f"  [{idx + 1}/{n_samples}] {participant_id} ({phenotype})...", end="", flush=True)

        try:
            points_df = load_sample_points(sample_id)
            n_total = len(points_df)

            # Plot all cell types, background first, malignant/inflammatory on top
            for cell_type in PLOT_ORDER:
                mask = points_df["classification_name"] == cell_type
                subset = points_df[mask]
                if len(subset) == 0:
                    continue
                ax.scatter(
                    subset["x"],
                    subset["y"],
                    c=CONSEP_COLORS[cell_type],
                    s=0.005,
                    alpha=0.7,
                    linewidths=0,
                    rasterized=True,
                )

            ax.invert_yaxis()
            ax.set_aspect("equal")

            n_malignant = (points_df["classification_name"] == "Malignant epithelial").sum()
            pct = 100 * n_malignant / n_total if n_total > 0 else 0
            print(f" {n_total:,} cells ({pct:.0f}% malignant)")

        except Exception as e:
            print(f" ERROR: {e}")
            ax.text(
                0.5, 0.5, "No data",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=8, color="gray",
            )

        # Style
        ax.set_xticks([])
        ax.set_yticks([])

        # Colored border by phenotype
        border_color = PHENOTYPE_COLORS.get(phenotype, "black")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(border_color)
            spine.set_linewidth(2)

        # Title: participant ID
        ax.set_title(
            participant_id.replace("HTA7_", ""),
            fontsize=8,
            fontweight="bold",
            color=border_color,
            pad=2,
        )

    # Hide unused axes
    for idx in range(n_samples, len(axes)):
        axes[idx].set_visible(False)

    # Cell type legend
    cell_legend = [
        Patch(facecolor=CONSEP_COLORS[ct], edgecolor="none", label=ct)
        for ct in CONSEP_ORDER
    ]
    # Phenotype border legend
    phenotype_legend = [
        Patch(facecolor="white", edgecolor=PHENOTYPE_COLORS[p], linewidth=2, label=p)
        for p in PHENOTYPE_ORDER
    ]

    fig.legend(
        handles=cell_legend + phenotype_legend,
        loc="lower center",
        ncol=5,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, -0.01),
    )

    fig.suptitle(
        "Spatial Cell Type Distribution Across All Samples",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )

    # Save
    output_dir = PROJECT_ROOT / "figures"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "malignant_mosaic_all_samples.png"
    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {output_path}")

    output_pdf = output_dir / "malignant_mosaic_all_samples.pdf"
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_pdf}")

    plt.close()
    print("Done!")


if __name__ == "__main__":
    create_mosaic()
