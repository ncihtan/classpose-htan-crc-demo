#!/usr/bin/env python3
"""
Generate 100% stacked bar chart of cell proportions by type for ConSeP segmentation.
Ordered by percentage of immune cells (descending).
"""

import matplotlib.pyplot as plt
from pathlib import Path

from load_cell_densities import (
    get_project_root,
    load_cell_densities,
    pivot_for_stacked_bar,
    CONSEP_IMMUNE_CLASSES,
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


def plot_stacked_100_bar_consep():
    """Create 100% stacked bar chart for ConSeP cell type proportions."""
    # Load and prepare data - sorted by immune percentage
    df = load_cell_densities("consep")
    pivoted = pivot_for_stacked_bar(
        df, sort_by_immune=True, immune_classes=CONSEP_IMMUNE_CLASSES
    )

    # Convert to percentages
    pivoted_pct = pivoted.div(pivoted.sum(axis=1), axis=0) * 100

    # Set up the figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Get colors in column order
    colors = [CONSEP_COLORS.get(col, "#999999") for col in pivoted_pct.columns]

    # Create stacked bar chart
    pivoted_pct.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8,
        edgecolor="white",
        linewidth=0.5,
    )

    # Formatting
    ax.set_xlabel("Sample (HTAN Participant ID) - ordered by % immune cells", fontsize=12)
    ax.set_ylabel("Proportion (%)", fontsize=12)
    ax.set_title(
        "Cell Type Proportions by Sample - ConSeP Segmentation (100% Stacked)",
        fontsize=14,
    )

    # Rotate x-axis labels
    plt.xticks(rotation=45, ha="right", fontsize=9)

    # Set y-axis limits
    ax.set_ylim(0, 100)

    # Legend
    ax.legend(
        title="Cell Type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=10,
    )

    # Add gridlines
    ax.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "stacked_100_bar_consep.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")

    plt.close()


if __name__ == "__main__":
    plot_stacked_100_bar_consep()
