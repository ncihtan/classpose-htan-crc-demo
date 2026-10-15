#!/usr/bin/env python3
"""
Generate scatter plot comparing total cell counts from CoNIC vs ConSeP segmentation.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

from load_cell_densities import (
    get_project_root,
    load_cell_densities,
    get_total_cells_per_sample,
)


def plot_scatter_conic_vs_consep():
    """Create scatter plot comparing total cells from CoNIC vs ConSeP."""
    # Load data for both models
    conic_df = load_cell_densities("conic")
    consep_df = load_cell_densities("consep")

    # Get total cells per sample
    conic_totals = get_total_cells_per_sample(conic_df)
    consep_totals = get_total_cells_per_sample(consep_df)

    # Align samples (in case any are missing)
    common_samples = conic_totals.index.intersection(consep_totals.index)
    conic_totals = conic_totals.loc[common_samples]
    consep_totals = consep_totals.loc[common_samples]

    # Create comparison DataFrame
    comparison = pd.DataFrame(
        {
            "conic": conic_totals,
            "consep": consep_totals,
            "participant_id": common_samples,
        }
    )

    # Calculate correlation
    r, p_value = stats.pearsonr(comparison["conic"], comparison["consep"])

    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Create scatter plot
    ax.scatter(
        comparison["conic"],
        comparison["consep"],
        s=80,
        alpha=0.7,
        c="#1f77b4",
        edgecolors="white",
        linewidth=0.5,
    )

    # Add identity line (y=x)
    max_val = max(comparison["conic"].max(), comparison["consep"].max())
    min_val = min(comparison["conic"].min(), comparison["consep"].min())
    ax.plot(
        [min_val, max_val],
        [min_val, max_val],
        "k--",
        alpha=0.5,
        label="Identity line (y=x)",
    )

    # Add regression line
    slope, intercept, r_val, p_val, std_err = stats.linregress(
        comparison["conic"], comparison["consep"]
    )
    x_line = np.array([min_val, max_val])
    y_line = slope * x_line + intercept
    ax.plot(
        x_line,
        y_line,
        "r-",
        alpha=0.7,
        label=f"Linear fit (slope={slope:.2f})",
    )

    # Formatting
    ax.set_xlabel("Total Cells - CoNIC Segmentation", fontsize=12)
    ax.set_ylabel("Total Cells - ConSeP Segmentation", fontsize=12)
    ax.set_title(
        f"Comparison of Total Cell Counts: CoNIC vs ConSeP\n"
        f"(Pearson r = {r:.3f}, p < {p_value:.2e}, n = {len(common_samples)} samples)",
        fontsize=14,
    )

    # Format axes with thousands separator
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x:,.0f}"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x:,.0f}"))

    # Make axes equal
    ax.set_aspect("equal", adjustable="box")

    # Add gridlines
    ax.grid(True, linestyle="--", alpha=0.3)

    # Legend
    ax.legend(loc="lower right", fontsize=10)

    plt.tight_layout()

    # Save figure
    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "scatter_conic_vs_consep.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")

    plt.close()

    # Print summary statistics
    print(f"\n=== Summary Statistics ===")
    print(f"Number of samples: {len(common_samples)}")
    print(f"CoNIC total cells: {comparison['conic'].sum():,.0f}")
    print(f"ConSeP total cells: {comparison['consep'].sum():,.0f}")
    print(f"Pearson correlation: r = {r:.4f} (p = {p_value:.2e})")
    print(f"Linear fit: slope = {slope:.4f}, intercept = {intercept:,.0f}")


if __name__ == "__main__":
    plot_scatter_conic_vs_consep()
