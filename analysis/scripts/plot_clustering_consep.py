#!/usr/bin/env python3
"""
Hierarchical clustering of samples based on proportional cell type composition.
ConSeP segmentation.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

from load_cell_densities import (
    get_project_root,
    load_cell_densities,
    pivot_for_stacked_bar,
)

# Color palette for ConSeP cell types
CONSEP_COLORS = {
    "Malignant epithelial": "#e41a1c",
    "Healthy epithelial": "#ff7f7f",
    "Stroma": "#377eb8",
    "Muscle": "#4daf4a",
    "Inflammatory": "#984ea3",
    "Other": "#999999",
}


def plot_clustering_consep():
    """Create clustered heatmap of ConSeP cell type proportions."""
    # Load and prepare data
    df = load_cell_densities("consep")
    pivoted = pivot_for_stacked_bar(df)

    # Convert to proportions (0-1)
    proportions = pivoted.div(pivoted.sum(axis=1), axis=0)

    # Create clustermap
    # Use seaborn clustermap which combines hierarchical clustering with heatmap
    col_colors = [CONSEP_COLORS.get(col, "#999999") for col in proportions.columns]

    g = sns.clustermap(
        proportions,
        method="ward",  # Ward's minimum variance method
        metric="euclidean",
        cmap="RdYlBu_r",
        figsize=(12, 14),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        col_colors=col_colors,
        col_cluster=True,  # Also cluster cell types
        row_cluster=True,  # Cluster samples
        linewidths=0.5,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        vmin=0,
        vmax=0.7,  # Cap at 70% for better color contrast
    )

    # Adjust labels
    g.ax_heatmap.set_xlabel("Cell Type", fontsize=12)
    g.ax_heatmap.set_ylabel("Sample (HTAN Participant ID)", fontsize=12)
    g.ax_heatmap.tick_params(axis="x", rotation=45)
    g.ax_heatmap.tick_params(axis="y", labelsize=9)

    # Add title
    g.fig.suptitle(
        "Hierarchical Clustering of Samples by Cell Type Proportions\nConSeP Segmentation",
        fontsize=14,
        y=1.02,
    )

    # Adjust colorbar label
    g.ax_cbar.set_ylabel("Proportion", fontsize=10)

    plt.tight_layout()

    # Save figure
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "clustering_consep.png"
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")

    plt.close()

    # Also print the cluster order for reference
    print("\n=== Sample Order (from clustering) ===")
    clustered_order = g.dendrogram_row.reordered_ind
    for i, idx in enumerate(clustered_order):
        sample = proportions.index[idx]
        print(f"{i+1:2d}. {sample}")


if __name__ == "__main__":
    plot_clustering_consep()
