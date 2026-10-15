#!/usr/bin/env python3
"""
Stacked bar chart showing immune cell distribution across compartments.

Each bar = one sample
Stacks = % of immune cells in tumor core, TSI, stroma, other
Sorted by phenotype then infiltration ratio.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root

# KNN config
CONFIG = {
    "n_neighbors": 30,
    "tumor_core_malignant_frac": 0.60,
    "tsi_malignant_min": 0.15,
    "stroma_stroma_frac": 0.40,
    "stroma_malignant_max": 0.15,
}

# Compartment colors (matching region plot)
COMPARTMENT_COLORS = {
    "tumor_core": "#08306b",  # dark blue
    "tsi": "#4292c6",         # medium blue
    "stroma": "#d4a017",      # mustard
    "other": "#bdbdbd",       # grey
}

PHENOTYPE_ORDER = ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]
PHENOTYPE_COLORS = {
    "Hot": "#b2182b",
    "Altered-immunosuppressed": "#f4a582",
    "Altered-excluded": "#92c5de",
    "Immune-desert": "#2166ac",
}


def load_sample_data(sample_id: str) -> pd.DataFrame:
    """Load sample data from parquet."""
    parquet_path = (
        SAMPLE_DATA_DIR / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
        / "points" / "cell_centroids" / "points.parquet"
    )
    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


def compute_immune_compartments(df: pd.DataFrame) -> dict:
    """Compute % of immune cells in each compartment."""
    coords = df[["x", "y"]].values
    cell_types = np.array(df["classification_name"].tolist())
    n_cells = len(df)

    # KNN classification
    k = CONFIG["n_neighbors"]
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1, workers=-1)
    neighbor_indices = indices[:, 1:]

    neighbor_types = cell_types[neighbor_indices]
    malignant_frac = (neighbor_types == "Malignant epithelial").sum(axis=1) / k
    stroma_frac = (neighbor_types == "Stroma").sum(axis=1) / k

    # Classify regions
    regions = np.zeros(n_cells, dtype=int)
    tumor_core_mask = malignant_frac >= CONFIG["tumor_core_malignant_frac"]
    regions[tumor_core_mask] = 1

    stroma_mask = (
        (stroma_frac >= CONFIG["stroma_stroma_frac"]) &
        (malignant_frac < CONFIG["stroma_malignant_max"]) &
        ~tumor_core_mask
    )
    regions[stroma_mask] = 3

    tsi_mask = (
        (malignant_frac >= CONFIG["tsi_malignant_min"]) &
        (malignant_frac < CONFIG["tumor_core_malignant_frac"]) &
        ~stroma_mask
    )
    regions[tsi_mask] = 2

    # Count immune cells per compartment
    is_immune = cell_types == "Inflammatory"
    total_immune = is_immune.sum()

    if total_immune == 0:
        return {"tumor_core": 0, "tsi": 0, "stroma": 0, "other": 0}

    return {
        "tumor_core": (is_immune & (regions == 1)).sum() / total_immune * 100,
        "tsi": (is_immune & (regions == 2)).sum() / total_immune * 100,
        "stroma": (is_immune & (regions == 3)).sum() / total_immune * 100,
        "other": (is_immune & (regions == 0)).sum() / total_immune * 100,
    }


def main():
    project_root = get_project_root()

    # Load v2 classification results
    metrics_df = pd.read_csv(project_root / "data" / "immune_phenotype_metrics_v2.csv")

    print("Computing immune cell compartment distribution...")

    results = []
    for _, row in tqdm(metrics_df.iterrows(), total=len(metrics_df), desc="Processing"):
        sample_id = row["sample_id"]
        try:
            df = load_sample_data(sample_id)
            compartments = compute_immune_compartments(df)
            compartments["sample_id"] = sample_id
            compartments["participant_id"] = row["participant_id"]
            compartments["phenotype"] = row["phenotype"]
            compartments["infiltration_ratio"] = row["infiltration_ratio"]
            compartments["density_overall"] = row["density_overall"]
            results.append(compartments)
        except Exception as e:
            print(f"  Error {sample_id}: {e}")

    results_df = pd.DataFrame(results)

    # Sort by phenotype order, then by infiltration ratio (descending)
    results_df["phenotype_order"] = results_df["phenotype"].map(
        {p: i for i, p in enumerate(PHENOTYPE_ORDER)}
    )
    results_df = results_df.sort_values(
        ["phenotype_order", "infiltration_ratio"],
        ascending=[True, False]
    ).reset_index(drop=True)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(16, 6))

    n_samples = len(results_df)
    x = np.arange(n_samples)
    bar_width = 0.8

    # Stack order: tumor_core (bottom), tsi, stroma, other (top)
    compartments = ["tumor_core", "tsi", "stroma", "other"]
    labels = ["Tumor core", "TSI (margin)", "Stroma", "Other"]

    bottom = np.zeros(n_samples)
    for compartment, label in zip(compartments, labels):
        values = results_df[compartment].values
        ax.bar(x, values, bar_width, bottom=bottom,
               color=COMPARTMENT_COLORS[compartment], label=label,
               edgecolor="white", linewidth=0.3)
        bottom += values

    # Add phenotype color indicators at bottom
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax.plot([i], [-3], marker="s", markersize=8,
                color=PHENOTYPE_COLORS[row["phenotype"]], clip_on=False)

    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(results_df["participant_id"], rotation=90, fontsize=8)

    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("% of Immune Cells", fontsize=12)
    ax.set_title("Distribution of Immune Cells Across Microenvironment Compartments",
                 fontsize=14, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.5, n_samples - 0.5)

    # Compartment legend
    ax.legend(loc="upper right", fontsize=10, title="Compartment")

    # Phenotype legend (manual)
    from matplotlib.patches import Patch
    phenotype_patches = [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in PHENOTYPE_ORDER]
    leg2 = ax.legend(handles=phenotype_patches, loc="upper left", fontsize=9,
                     title="Phenotype", framealpha=0.9)
    ax.add_artist(ax.legend(loc="upper right", fontsize=10, title="Compartment"))

    # Add vertical lines between phenotype groups
    current_phenotype = results_df["phenotype"].iloc[0]
    for i, phenotype in enumerate(results_df["phenotype"]):
        if phenotype != current_phenotype:
            ax.axvline(x=i - 0.5, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
            current_phenotype = phenotype

    plt.tight_layout()

    output_path = project_root / "figures" / "immune_compartment_stacked_bars.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"\nSaved to: {output_path}")

    # Also save the data
    results_df.to_csv(project_root / "data" / "immune_compartment_distribution.csv", index=False)
    print(f"Data saved to: {project_root / 'data' / 'immune_compartment_distribution.csv'}")

    # Print summary
    print("\nMean % immune cells by compartment and phenotype:")
    summary = results_df.groupby("phenotype")[["tumor_core", "tsi", "stroma", "other"]].mean()
    summary = summary.reindex(PHENOTYPE_ORDER)
    print(summary.round(1))


if __name__ == "__main__":
    main()
