#!/usr/bin/env python3
"""
Hierarchical clustering of samples using Centered Log-Ratio (CLR) transformation.
This addresses the compositional data problem where proportions must sum to 100%.

The CLR transformation maps compositional data from the simplex to unconstrained
Euclidean space, making standard clustering methods appropriate.
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

# Color palettes
CONIC_COLORS = {
    "Epithelial": "#e41a1c",
    "Connective": "#377eb8",
    "Lymphocyte": "#4daf4a",
    "Plasma cell": "#984ea3",
    "Neutrophil": "#ff7f00",
    "Eosinophil": "#ffff33",
}

CONSEP_COLORS = {
    "Malignant epithelial": "#e41a1c",
    "Healthy epithelial": "#ff7f7f",
    "Stroma": "#377eb8",
    "Muscle": "#4daf4a",
    "Inflammatory": "#984ea3",
    "Other": "#999999",
}


def clr_transform(df):
    """
    Apply Centered Log-Ratio (CLR) transformation to compositional data.

    CLR(x) = log(x / geometric_mean(x))

    This transforms data from the simplex (where values sum to 1) to
    unconstrained Euclidean space, making standard statistical methods valid.

    Note: Zeros are replaced with a small value (pseudocount) before transformation.
    """
    # Replace zeros with small pseudocount (half the minimum non-zero value)
    df_nonzero = df.replace(0, np.nan)
    min_val = df_nonzero.min().min()
    pseudocount = min_val / 2 if not np.isnan(min_val) else 1e-6
    df_pseudo = df.replace(0, pseudocount)

    # Calculate geometric mean per row
    log_data = np.log(df_pseudo)
    geometric_mean = np.exp(log_data.mean(axis=1))

    # CLR transformation: log(x / geometric_mean)
    clr_data = log_data.subtract(np.log(geometric_mean), axis=0)

    return clr_data


def plot_clustering_clr(model):
    """Create clustered heatmap using CLR-transformed data."""
    colors = CONIC_COLORS if model == "conic" else CONSEP_COLORS

    # Load and prepare data
    df = load_cell_densities(model)
    pivoted = pivot_for_stacked_bar(df)

    # Convert to proportions
    proportions = pivoted.div(pivoted.sum(axis=1), axis=0)

    # Apply CLR transformation
    clr_data = clr_transform(proportions)

    # Create clustermap with CLR-transformed data
    col_colors = [colors.get(col, "#999999") for col in clr_data.columns]

    g = sns.clustermap(
        clr_data,
        method="ward",
        metric="euclidean",  # Now appropriate after CLR transformation
        cmap="RdBu_r",  # Diverging colormap since CLR values center around 0
        figsize=(12, 14),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        col_colors=col_colors,
        col_cluster=True,
        row_cluster=True,
        linewidths=0.5,
        linecolor="white",
        xticklabels=True,
        yticklabels=True,
        center=0,  # Center colormap at 0 (geometric mean)
    )

    # Adjust labels
    g.ax_heatmap.set_xlabel("Cell Type", fontsize=12)
    g.ax_heatmap.set_ylabel("Sample (HTAN Participant ID)", fontsize=12)
    g.ax_heatmap.tick_params(axis="x", rotation=45)
    g.ax_heatmap.tick_params(axis="y", labelsize=9)

    model_name = "CoNIC" if model == "conic" else "ConSeP"
    g.fig.suptitle(
        f"Hierarchical Clustering with CLR Transformation\n{model_name} Segmentation",
        fontsize=14,
        y=1.02,
    )

    # Adjust colorbar label
    g.ax_cbar.set_ylabel("CLR value\n(log-ratio to geometric mean)", fontsize=9)

    plt.tight_layout()

    # Save figure
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f"clustering_clr_{model}.png"
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")

    plt.close()

    # Print interpretation
    print(f"\n=== CLR Clustering Results ({model_name}) ===")
    print("\nInterpretation of CLR values:")
    print("  - Positive values: cell type is MORE abundant than average (geometric mean)")
    print("  - Negative values: cell type is LESS abundant than average")
    print("  - Zero: cell type is at the geometric mean abundance")
    print("\nThis approach avoids spurious correlations from compositional data.")

    # Print cluster order
    print("\n=== Sample Order (from clustering) ===")
    clustered_order = g.dendrogram_row.reordered_ind
    for i, idx in enumerate(clustered_order):
        sample = clr_data.index[idx]
        print(f"{i+1:2d}. {sample}")

    return clr_data


def compare_clustering_methods():
    """
    Compare clustering results between raw proportions and CLR transformation.
    Shows why CLR is more appropriate for compositional data.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    for i, model in enumerate(["consep", "conic"]):
        colors = CONSEP_COLORS if model == "consep" else CONIC_COLORS
        model_name = "ConSeP" if model == "consep" else "CoNIC"

        # Load data
        df = load_cell_densities(model)
        pivoted = pivot_for_stacked_bar(df)
        proportions = pivoted.div(pivoted.sum(axis=1), axis=0)
        clr_data = clr_transform(proportions)

        # Calculate correlation matrices
        # Raw proportions - will show spurious negative correlations
        corr_raw = proportions.corr()
        # CLR-transformed - more meaningful correlations
        corr_clr = clr_data.corr()

        # Plot raw correlation
        ax1 = axes[i, 0]
        sns.heatmap(
            corr_raw,
            ax=ax1,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            square=True,
        )
        ax1.set_title(f"{model_name}: Raw Proportions\n(contains spurious correlations)", fontsize=11)
        ax1.tick_params(axis="x", rotation=45)
        ax1.tick_params(axis="y", rotation=0)

        # Plot CLR correlation
        ax2 = axes[i, 1]
        sns.heatmap(
            corr_clr,
            ax=ax2,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            square=True,
        )
        ax2.set_title(f"{model_name}: CLR-Transformed\n(compositionally appropriate)", fontsize=11)
        ax2.tick_params(axis="x", rotation=45)
        ax2.tick_params(axis="y", rotation=0)

    fig.suptitle(
        "Correlation Between Cell Types: Raw Proportions vs CLR Transformation\n"
        "CLR removes spurious negative correlations inherent in compositional data",
        fontsize=14,
        y=1.02,
    )

    plt.tight_layout()

    # Save
    output_dir = get_project_root() / "figures"
    output_path = output_dir / "correlation_comparison_raw_vs_clr.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {output_path}")
    plt.close()

    print("\n=== Why CLR Transformation Matters ===")
    print("In raw proportions, cell types show negative correlations simply because")
    print("if one type increases, others must decrease (they sum to 100%).")
    print("The CLR transformation removes this artifact, revealing true biological")
    print("relationships between cell types.")


if __name__ == "__main__":
    # Generate CLR-based clustering for both models
    print("=" * 60)
    print("Compositional Data Analysis with CLR Transformation")
    print("=" * 60)

    plot_clustering_clr("consep")
    print()
    plot_clustering_clr("conic")

    # Compare correlation matrices
    print("\n" + "=" * 60)
    print("Comparing Raw vs CLR Correlations")
    print("=" * 60)
    compare_clustering_methods()
