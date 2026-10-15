#!/usr/bin/env python3
"""
Plot KNN-based region assignments for hot and cold samples.

Shows individual cells colored by their microenvironment classification
based on neighborhood composition (tumor core, TSI, stroma, other).
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root

# KNN config (must match v2 script)
CONFIG = {
    "n_neighbors": 30,
    "tumor_core_malignant_frac": 0.60,
    "tsi_malignant_min": 0.15,
    "stroma_stroma_frac": 0.40,
    "stroma_malignant_max": 0.15,
}

REGION_COLORS = {
    0: "#bdbdbd",  # Other - grey
    1: "#d62728",  # Tumor core - red
    2: "#ff7f0e",  # TSI - orange
    3: "#2ca02c",  # Stroma - green
}

REGION_LABELS = {
    0: "Other",
    1: "Tumor core",
    2: "TSI (margin)",
    3: "Stroma",
}

# Representative samples
SAMPLES = {
    "Hot": "dg_4DFC_B0A48536-0171-4411-A200-88344086F3D4",  # HTA7_941
    "Cold": "dg_4DFC_8C77314D-9FAC-4917-83CE-52A3B70FE014",  # HTA7_972
}


def load_sample_data(sample_id: str) -> pd.DataFrame:
    """Load sample data from parquet."""
    parquet_path = (
        SAMPLE_DATA_DIR / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
        / "points" / "cell_centroids" / "points.parquet"
    )
    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


def classify_regions_knn(df: pd.DataFrame) -> np.ndarray:
    """Classify cell microenvironments using KNN."""
    coords = df[["x", "y"]].values
    cell_types = np.array(df["classification_name"].tolist())
    n_cells = len(df)

    # Build KNN graph
    k = CONFIG["n_neighbors"]
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1, workers=-1)
    neighbor_indices = indices[:, 1:]

    # Compute neighborhood composition
    neighbor_types = cell_types[neighbor_indices]
    malignant_frac = (neighbor_types == "Malignant epithelial").sum(axis=1) / k
    stroma_frac = (neighbor_types == "Stroma").sum(axis=1) / k

    # Classify regions
    regions = np.zeros(n_cells, dtype=int)

    # Tumor core: high malignant neighborhood
    tumor_core_mask = malignant_frac >= CONFIG["tumor_core_malignant_frac"]
    regions[tumor_core_mask] = 1

    # Stroma: high stroma, low malignant
    stroma_mask = (
        (stroma_frac >= CONFIG["stroma_stroma_frac"]) &
        (malignant_frac < CONFIG["stroma_malignant_max"]) &
        ~tumor_core_mask
    )
    regions[stroma_mask] = 3

    # TSI: interface region with mixed composition
    tsi_mask = (
        (malignant_frac >= CONFIG["tsi_malignant_min"]) &
        (malignant_frac < CONFIG["tumor_core_malignant_frac"]) &
        ~stroma_mask
    )
    regions[tsi_mask] = 2

    return regions


def plot_regions(sample_id: str, label: str, ax: plt.Axes, point_size: float = 0.05):
    """Plot cells colored by KNN-defined region."""
    print(f"  Loading {label} ({sample_id})...")
    df = load_sample_data(sample_id)

    print(f"  Computing KNN regions for {len(df):,} cells...")
    regions = classify_regions_knn(df)

    # Plot by region (other first, then stroma, TSI, tumor core)
    for region_code in [0, 3, 2, 1]:
        mask = regions == region_code
        if mask.sum() == 0:
            continue
        ax.scatter(
            df.loc[mask, "x"],
            df.loc[mask, "y"],
            c=REGION_COLORS[region_code],
            s=point_size,
            alpha=0.7,
            linewidths=0,
            rasterized=True,
        )

    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])

    # Count regions
    counts = {REGION_LABELS[i]: (regions == i).sum() for i in range(4)}
    title = f"{label}\n"
    title += f"Core: {counts['Tumor core']:,} | TSI: {counts['TSI (margin)']:,} | Stroma: {counts['Stroma']:,}"
    ax.set_title(title, fontsize=11, fontweight="bold")

    return counts


def main():
    print("Creating KNN region plots...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, (label, sample_id) in zip(axes, SAMPLES.items()):
        plot_regions(sample_id, label, ax)

    # Legend
    legend_elements = [
        Patch(facecolor=REGION_COLORS[1], label="Tumor core (>60% malignant neighbors)"),
        Patch(facecolor=REGION_COLORS[2], label="TSI/Margin (15-60% malignant)"),
        Patch(facecolor=REGION_COLORS[3], label="Stroma (>40% stroma, <15% malignant)"),
        Patch(facecolor=REGION_COLORS[0], label="Other"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=4,
        fontsize=10,
        frameon=True,
        bbox_to_anchor=(0.5, 0.02),
    )

    plt.suptitle(
        "KNN-Based Microenvironment Classification (K=30)",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])

    output_path = get_project_root() / "figures" / "knn_region_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
