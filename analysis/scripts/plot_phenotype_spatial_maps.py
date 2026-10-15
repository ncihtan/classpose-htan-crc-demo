#!/usr/bin/env python3
"""
Plot spatial cell type density maps for representative samples of each immune phenotype.

Creates filled contour plots showing the spatial distribution of:
- Malignant epithelial (tumor): dark blue
- Inflammatory (immune): magenta
- Stroma: mustard yellow
- Other tissues: greyscale
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root


# Cell type color schemes
CELL_COLORS = {
    "Malignant epithelial": {
        "cmap": LinearSegmentedColormap.from_list("tumor", ["white", "#08306b"]),  # dark blue
        "label": "Tumor",
    },
    "Inflammatory": {
        "cmap": LinearSegmentedColormap.from_list("immune", ["white", "#c51b8a"]),  # magenta
        "label": "Immune",
    },
    "Stroma": {
        "cmap": LinearSegmentedColormap.from_list("stroma", ["white", "#d4a017"]),  # mustard
        "label": "Stroma",
    },
    "Healthy epithelial": {
        "cmap": LinearSegmentedColormap.from_list("healthy", ["white", "#737373"]),  # grey
        "label": "Healthy Epi.",
    },
    "Muscle": {
        "cmap": LinearSegmentedColormap.from_list("muscle", ["white", "#525252"]),  # darker grey
        "label": "Muscle",
    },
    "Other": {
        "cmap": LinearSegmentedColormap.from_list("other", ["white", "#969696"]),  # light grey
        "label": "Other",
    },
}

# Representative samples for each phenotype
PHENOTYPE_SAMPLES = {
    "Hot": "dg_4DFC_B0A48536-0171-4411-A200-88344086F3D4",
    "Altered-immunosuppressed": "dg_4DFC_96A5875F-826A-4917-A7DE-18E489024DF9",
    "Altered-excluded": "dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF",
    "Immune-desert": "dg_4DFC_8C77314D-9FAC-4917-83CE-52A3B70FE014",
}


def load_sample_points(sample_id: str) -> pd.DataFrame:
    """Load cell centroids from parquet."""
    parquet_path = (
        SAMPLE_DATA_DIR
        / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
        / "points"
        / "cell_centroids"
        / "points.parquet"
    )
    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


def compute_density_grid(points_df: pd.DataFrame, cell_type: str, n_bins: int = 100):
    """Compute 2D density histogram for a cell type."""
    mask = points_df["classification_name"] == cell_type
    subset = points_df[mask]

    x_range = (points_df["x"].min(), points_df["x"].max())
    y_range = (points_df["y"].min(), points_df["y"].max())

    if len(subset) == 0:
        hist = np.zeros((n_bins, n_bins))
    else:
        hist, _, _ = np.histogram2d(
            subset["x"].values,
            subset["y"].values,
            bins=n_bins,
            range=[x_range, y_range],
        )

    return hist.T  # Transpose for correct orientation


def plot_single_sample(
    sample_id: str,
    phenotype: str,
    ax: plt.Axes,
    point_size: float = 0.1,
    subsample: int | None = None,
):
    """Plot individual cells colored by cell type."""
    points_df = load_sample_points(sample_id)

    # Get sample metadata
    classification = pd.read_csv(get_project_root() / "data" / "immune_phenotype_classification.csv")
    sample_info = classification[classification["sample_id"] == sample_id].iloc[0]

    # Cell type colors
    cell_colors = {
        "Malignant epithelial": "#08306b",  # dark blue
        "Inflammatory": "#c51b8a",           # magenta
        "Stroma": "#d4a017",                 # mustard
        "Healthy epithelial": "#737373",     # grey
        "Muscle": "#525252",                 # darker grey
        "Other": "#bdbdbd",                  # light grey
    }

    # Plot order: background first, then stroma, then tumor, then immune on top
    plot_order = ["Other", "Muscle", "Healthy epithelial", "Stroma", "Malignant epithelial", "Inflammatory"]

    for cell_type in plot_order:
        mask = points_df["classification_name"] == cell_type
        subset = points_df[mask]

        if len(subset) == 0:
            continue

        # Subsample if requested (for speed)
        if subsample and len(subset) > subsample:
            subset = subset.sample(n=subsample, random_state=42)

        ax.scatter(
            subset["x"],
            subset["y"],
            c=cell_colors.get(cell_type, "#999999"),
            s=point_size,
            alpha=0.8,
            linewidths=0,
            rasterized=True,  # For faster rendering
        )

    # Invert y-axis for image-like orientation
    ax.invert_yaxis()

    # Add title with phenotype and metrics
    title = f"{phenotype}\n{sample_info['participant_id']}: {sample_info['inflammatory_pct_overall']:.1f}% immune"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")

    # Clean up spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.5)


def create_phenotype_comparison():
    """Create a 2x2 figure comparing all phenotypes."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    axes = axes.flatten()

    phenotype_order = ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]

    for ax, phenotype in zip(axes, phenotype_order):
        sample_id = PHENOTYPE_SAMPLES[phenotype]
        print(f"Plotting {phenotype} ({sample_id})...")
        plot_single_sample(sample_id, phenotype, ax, point_size=0.05)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#08306b", label="Tumor"),
        Patch(facecolor="#c51b8a", label="Immune"),
        Patch(facecolor="#d4a017", label="Stroma"),
        Patch(facecolor="#737373", label="Other tissue"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=4,
        fontsize=11,
        frameon=True,
        bbox_to_anchor=(0.5, 0.02),
    )

    plt.suptitle(
        "Spatial Cell Type Distribution by Immune Phenotype",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])

    # Save
    output_path = get_project_root() / "figures" / "immune_phenotype_spatial_comparison.png"
    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nSaved to: {output_path}")


def create_individual_panels():
    """Create individual high-res panels for each phenotype."""
    output_dir = get_project_root() / "figures"

    for phenotype, sample_id in PHENOTYPE_SAMPLES.items():
        print(f"Plotting {phenotype}...")

        fig, ax = plt.subplots(figsize=(10, 10))
        plot_single_sample(sample_id, phenotype, ax, point_size=0.15)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor="#08306b", label="Tumor"),
            Patch(facecolor="#c51b8a", label="Immune"),
            Patch(facecolor="#d4a017", label="Stroma"),
            Patch(facecolor="#737373", label="Other tissue"),
        ]
        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=10,
            frameon=True,
            facecolor="white",
            edgecolor="gray",
        )

        # Save
        safe_name = phenotype.lower().replace("-", "_")
        output_path = output_dir / f"spatial_map_{safe_name}.png"
        plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"  Saved: {output_path}")


if __name__ == "__main__":
    print("Creating spatial cell type maps...")
    print("=" * 50)

    # Create comparison figure
    create_phenotype_comparison()

    # Create individual high-res panels
    print("\nCreating individual panels...")
    create_individual_panels()

    print("\nDone!")
