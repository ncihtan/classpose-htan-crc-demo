#!/usr/bin/env python3
"""
Create paired panels showing cell types and KNN regions for each phenotype.

4x2 grid:
- Rows: Hot, Altered-immunosuppressed, Altered-excluded, Immune-desert
- Columns: Cell types | KNN regions
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

# KNN config
CONFIG = {
    "n_neighbors": 30,
    "tumor_core_malignant_frac": 0.60,
    "tsi_malignant_min": 0.15,
    "stroma_stroma_frac": 0.40,
    "stroma_malignant_max": 0.15,
}

# Cell type colors
CELL_COLORS = {
    "Malignant epithelial": "#08306b",  # dark blue (tumor)
    "Inflammatory": "#c51b8a",           # magenta (immune)
    "Stroma": "#d4a017",                 # mustard
    "Healthy epithelial": "#969696",     # grey
    "Muscle": "#737373",                 # darker grey
    "Other": "#bdbdbd",                  # light grey
}

# Region colors
REGION_COLORS = {
    0: "#bdbdbd",  # Other - grey
    1: "#08306b",  # Tumor core - dark blue
    2: "#4292c6",  # TSI - medium blue
    3: "#d4a017",  # Stroma - mustard
}

REGION_LABELS = {
    0: "Other",
    1: "Tumor core",
    2: "TSI (margin)",
    3: "Stroma",
}

# Representative samples for each phenotype (from v2 results - updated for tumor_associated_pct)
PHENOTYPE_SAMPLES = {
    "Hot": "dg_4DFC_6E7B033F-549D-4063-8F75-60D6418F616C",  # HTA7_908 - 9457/mm², 41.5% tumor-assoc
    "Altered-immunosuppressed": "dg_4DFC_98B41CE1-C47F-4C37-BB1D-F4B49E1C4619",  # HTA7_932 - 2817/mm², 33.9% tumor-assoc
    "Altered-excluded": "dg_4DFC_1A504F88-B55A-40DB-B090-CEE00B3105EF",  # HTA7_917 - 17728/mm², 4.8% tumor-assoc
    "Immune-desert": "dg_4DFC_ACAB6731-71AC-4351-8789-16141B385EC4",  # HTA7_933 - 1945/mm², 27.5% tumor-assoc
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

    k = CONFIG["n_neighbors"]
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1, workers=-1)
    neighbor_indices = indices[:, 1:]

    neighbor_types = cell_types[neighbor_indices]
    malignant_frac = (neighbor_types == "Malignant epithelial").sum(axis=1) / k
    stroma_frac = (neighbor_types == "Stroma").sum(axis=1) / k

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

    return regions


def plot_cell_types(df: pd.DataFrame, ax: plt.Axes, point_size: float = 0.02):
    """Plot cells colored by cell type."""
    # Plot order: other first, then stroma, tumor, immune on top
    plot_order = ["Other", "Muscle", "Healthy epithelial", "Stroma", "Malignant epithelial", "Inflammatory"]

    for cell_type in plot_order:
        mask = df["classification_name"] == cell_type
        if mask.sum() == 0:
            continue
        ax.scatter(
            df.loc[mask, "x"],
            df.loc[mask, "y"],
            c=CELL_COLORS.get(cell_type, "#999999"),
            s=point_size,
            alpha=0.8,
            linewidths=0,
            rasterized=True,
        )

    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_regions(df: pd.DataFrame, regions: np.ndarray, ax: plt.Axes, point_size: float = 0.02):
    """Plot cells colored by KNN-defined region."""
    coords = df[["x", "y"]].values

    # Plot order: other first, then stroma, TSI, tumor core
    for region_code in [0, 3, 2, 1]:
        mask = regions == region_code
        if mask.sum() == 0:
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
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


def get_sample_info(sample_id: str) -> dict:
    """Get sample metrics from v2 results."""
    project_root = get_project_root()
    metrics_path = project_root / "data" / "immune_phenotype_metrics_v2.csv"
    if metrics_path.exists():
        df = pd.read_csv(metrics_path)
        row = df[df["sample_id"] == sample_id]
        if len(row) > 0:
            return row.iloc[0].to_dict()
    return {}


def main():
    print("Creating paired phenotype panels (4x2 grid)...")

    fig, axes = plt.subplots(4, 2, figsize=(14, 24))

    phenotype_order = ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]

    for row, phenotype in enumerate(phenotype_order):
        sample_id = PHENOTYPE_SAMPLES[phenotype]
        print(f"Processing {phenotype}...")

        # Load data
        df = load_sample_data(sample_id)
        regions = classify_regions_knn(df)

        # Get metrics
        info = get_sample_info(sample_id)
        participant = info.get("participant_id", "")
        density = info.get("density_overall", 0)
        tumor_assoc = info.get("tumor_associated_pct", 0)

        # Left: cell types
        ax_cells = axes[row, 0]
        plot_cell_types(df, ax_cells)
        ax_cells.set_title(f"{phenotype}\n{participant}: {density:.0f}/mm²", fontsize=11, fontweight="bold")
        if row == 0:
            ax_cells.set_ylabel("Cell Types", fontsize=12, fontweight="bold")

        # Right: regions
        ax_regions = axes[row, 1]
        plot_regions(df, regions, ax_regions)

        # Count regions
        n_core = (regions == 1).sum()
        n_tsi = (regions == 2).sum()
        n_stroma = (regions == 3).sum()
        ax_regions.set_title(f"Tumor-associated: {tumor_assoc:.1f}%\nCore: {n_core:,} | TSI: {n_tsi:,} | Stroma: {n_stroma:,}",
                            fontsize=10)
        if row == 0:
            ax_regions.set_ylabel("KNN Regions", fontsize=12, fontweight="bold")

    # Cell type legend
    cell_legend = [
        Patch(facecolor="#08306b", label="Tumor"),
        Patch(facecolor="#c51b8a", label="Immune"),
        Patch(facecolor="#d4a017", label="Stroma"),
        Patch(facecolor="#969696", label="Other tissue"),
    ]

    # Region legend
    region_legend = [
        Patch(facecolor="#08306b", label="Tumor core"),
        Patch(facecolor="#4292c6", label="TSI (margin)"),
        Patch(facecolor="#d4a017", label="Stroma"),
        Patch(facecolor="#bdbdbd", label="Other"),
    ]

    # Add legends at bottom
    fig.legend(handles=cell_legend, loc="lower left", ncol=4, fontsize=10,
               frameon=True, bbox_to_anchor=(0.05, 0.01))
    fig.legend(handles=region_legend, loc="lower right", ncol=4, fontsize=10,
               frameon=True, bbox_to_anchor=(0.95, 0.01))

    plt.suptitle(
        "Immune Phenotype Classification: Cell Types vs KNN Microenvironment Regions",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.99])

    output_path = get_project_root() / "figures" / "phenotype_paired_panels.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()
