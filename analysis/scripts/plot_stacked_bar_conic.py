#!/usr/bin/env python3
"""
Generate stacked bar chart of cell counts by type for CoNIC segmentation.
Ordered by percentage of immune cells (descending).
"""

import matplotlib.pyplot as plt
from pathlib import Path

from load_cell_densities import (
    get_project_root,
    load_cell_densities,
    pivot_for_stacked_bar,
    CONIC_IMMUNE_CLASSES,
)

# Color palette for CoNIC cell types
CONIC_COLORS = {
    "Epithelial": "#e41a1c",
    "Connective": "#377eb8",
    "Lymphocyte": "#4daf4a",
    "Plasma cell": "#984ea3",
    "Neutrophil": "#ff7f00",
    "Eosinophil": "#ffff33",
}


def plot_stacked_bar_conic():
    """Create stacked bar chart for CoNIC cell type counts."""
    # Load and prepare data - sorted by immune percentage
    df = load_cell_densities("conic")
    pivoted = pivot_for_stacked_bar(
        df, sort_by_immune=True, immune_classes=CONIC_IMMUNE_CLASSES
    )

    # Set up the figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Get colors in column order
    colors = [CONIC_COLORS.get(col, "#999999") for col in pivoted.columns]

    # Create stacked bar chart
    pivoted.plot(
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
    ax.set_ylabel("Cell Count", fontsize=12)
    ax.set_title("Cell Type Distribution by Sample - CoNIC Segmentation", fontsize=14)

    # Rotate x-axis labels
    plt.xticks(rotation=45, ha="right", fontsize=9)

    # Format y-axis with thousands separator
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x:,.0f}"))

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
    output_path = output_dir / "stacked_bar_conic.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")

    plt.close()


if __name__ == "__main__":
    plot_stacked_bar_conic()
