#!/usr/bin/env python3
"""
Classify tumor immune phenotypes for HTAN CRC samples.

Uses spatial analysis of cell distributions from ConSeP SpatialData zarrs
to classify samples into immune phenotypes based on the Immunoscore system:
- Immune-desert: No inflammatory cells anywhere
- Altered-excluded: Inflammatory in stroma, blocked from tumor
- Altered-immunosuppressed: Low inflammatory throughout, no clustering
- Hot: Robust tumor infiltration, inflammatory clusters detected

Usage:
    python scripts/classify_immune_phenotypes.py
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation
from sklearn.cluster import DBSCAN
from tqdm import tqdm

# Filter warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root, load_cell_densities, get_sample_metadata

# ConSeP cell type order
CONSEP_ORDER = [
    "Healthy epithelial",
    "Malignant epithelial",
    "Stroma",
    "Inflammatory",
    "Muscle",
    "Other",
]

# Configuration
CONFIG = {
    # Grid parameters
    "n_bins": 50,
    "min_cells_per_bin": 10,
    # Region thresholds
    "core_percentile": 75,  # Top 25% malignant-containing bins
    "margin_min_malignant_pct": 0.05,
    "stroma_max_malignant_pct": 0.10,
    "stroma_min_stroma_pct": 0.30,
    # Cluster detection parameters (DBSCAN) - kept as supplementary metric
    "cluster_eps": 100,  # Max distance between samples in cluster
    "cluster_min_samples": 30,  # Min samples for core point
    "cluster_min_size": 30,  # Min cells for valid cluster
    # Phenotype classification - uses percentile-based thresholds
    # These are computed dynamically from cohort distribution
    "inflammatory_low_percentile": 33,  # Bottom third = low inflammatory
    "inflammatory_high_percentile": 67,  # Top third = high inflammatory
    "infiltration_low_percentile": 33,  # Bottom third = poor infiltration
}


def load_sample_spatial_data(zarr_path: Path) -> pd.DataFrame:
    """
    Load cell centroids from a SpatialData zarr.

    Uses direct parquet reading for speed (~45x faster than spatialdata).

    Parameters
    ----------
    zarr_path : Path
        Path to the .zarr directory

    Returns
    -------
    pd.DataFrame
        DataFrame with x, y, classification_name columns
    """
    # Direct parquet reading is much faster than going through spatialdata
    parquet_path = zarr_path / "points" / "cell_centroids" / "points.parquet"
    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


def compute_region_grid(
    points_df: pd.DataFrame,
    n_bins: int = 50,
    min_cells_per_bin: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Divide slide into grid and classify each bin as tumor_core/margin/stroma/other.

    Parameters
    ----------
    points_df : pd.DataFrame
        DataFrame with x, y, classification_name
    n_bins : int
        Number of bins in each dimension
    min_cells_per_bin : int
        Minimum cells for a bin to be classified

    Returns
    -------
    tuple
        (region_grid, x_bins, y_bins, histograms)
        - region_grid: 2D array with region labels (0=other, 1=core, 2=margin, 3=stroma)
        - x_bins, y_bins: bin edges
        - histograms: dict of 2D histograms per cell type
    """
    # Determine slide dimensions
    x_min, x_max = points_df["x"].min(), points_df["x"].max()
    y_min, y_max = points_df["y"].min(), points_df["y"].max()

    x_bins = np.linspace(x_min, x_max, n_bins + 1)
    y_bins = np.linspace(y_min, y_max, n_bins + 1)

    # Count cells by type per bin
    histograms = {}
    for cell_type in CONSEP_ORDER:
        mask = points_df["classification_name"] == cell_type
        hist, _, _ = np.histogram2d(
            points_df.loc[mask, "x"].values,
            points_df.loc[mask, "y"].values,
            bins=[x_bins, y_bins],
        )
        histograms[cell_type] = hist

    # Compute total cells per bin and percentages
    total_per_bin = sum(histograms.values())
    total_per_bin_safe = total_per_bin + 1e-6

    malignant_pct = histograms["Malignant epithelial"] / total_per_bin_safe
    stroma_pct = histograms["Stroma"] / total_per_bin_safe

    # Initialize region grid: 0 = other, 1 = core, 2 = margin, 3 = stroma
    region_grid = np.zeros_like(total_per_bin, dtype=int)

    # Find bins with significant malignant cells
    malignant_bins = malignant_pct > 0.1
    if malignant_bins.sum() > 0:
        # Tumor core: top 25% malignant density among malignant-containing bins
        malignant_threshold = np.percentile(
            malignant_pct[malignant_bins], CONFIG["core_percentile"]
        )
        core_mask = malignant_pct >= malignant_threshold
        region_grid[core_mask] = 1  # core

        # Tumor margin: adjacent to core, moderate malignant
        core_dilated = binary_dilation(core_mask, iterations=1)
        margin_mask = (
            core_dilated
            & ~core_mask
            & (malignant_pct > CONFIG["margin_min_malignant_pct"])
        )
        region_grid[margin_mask] = 2  # margin

    # Stroma-only: low malignant, high stroma, sufficient cells
    stroma_mask = (
        (malignant_pct < CONFIG["stroma_max_malignant_pct"])
        & (stroma_pct > CONFIG["stroma_min_stroma_pct"])
        & (total_per_bin > min_cells_per_bin)
    )
    # Don't overwrite core or margin
    stroma_mask = stroma_mask & (region_grid == 0)
    region_grid[stroma_mask] = 3  # stroma

    return region_grid, x_bins, y_bins, histograms


def compute_regional_metrics(
    points_df: pd.DataFrame,
    region_grid: np.ndarray,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
) -> dict:
    """
    Compute inflammatory cell metrics for each region.

    Parameters
    ----------
    points_df : pd.DataFrame
        DataFrame with x, y, classification_name
    region_grid : np.ndarray
        2D array with region labels
    x_bins, y_bins : np.ndarray
        Bin edges

    Returns
    -------
    dict
        Computed metrics
    """
    # Get inflammatory cells
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    # Digitize inflammatory cell positions
    x_idx = np.digitize(inflammatory["x"].values, x_bins) - 1
    y_idx = np.digitize(inflammatory["y"].values, y_bins) - 1

    # Clip to valid range
    x_idx = np.clip(x_idx, 0, region_grid.shape[0] - 1)
    y_idx = np.clip(y_idx, 0, region_grid.shape[1] - 1)

    # Get region for each inflammatory cell
    cell_regions = region_grid[x_idx, y_idx]

    # Calculate bin area
    bin_width = x_bins[1] - x_bins[0]
    bin_height = y_bins[1] - y_bins[0]
    bin_area = bin_width * bin_height

    # Count cells and area per region
    metrics = {
        "total_cells": len(points_df),
        "inflammatory_count": len(inflammatory),
        "inflammatory_pct_overall": len(inflammatory) / len(points_df) * 100
        if len(points_df) > 0
        else 0,
    }

    # Per-region metrics
    region_names = {0: "other", 1: "core", 2: "margin", 3: "stroma"}
    for region_id, region_name in region_names.items():
        n_bins_region = np.sum(region_grid == region_id)
        n_cells_region = np.sum(cell_regions == region_id)
        area_region = n_bins_region * bin_area

        metrics[f"inflammatory_count_{region_name}"] = n_cells_region
        metrics[f"n_bins_{region_name}"] = n_bins_region
        metrics[f"area_{region_name}"] = area_region
        metrics[f"inflammatory_density_{region_name}"] = (
            n_cells_region / area_region if area_region > 0 else 0
        )

    # Compute ratios
    core_density = metrics["inflammatory_density_core"]
    stroma_density = metrics["inflammatory_density_stroma"]
    margin_density = metrics["inflammatory_density_margin"]

    metrics["infiltration_ratio"] = (
        core_density / stroma_density if stroma_density > 0 else float("inf")
    )
    metrics["margin_core_ratio"] = (
        margin_density / core_density if core_density > 0 else float("inf")
    )

    return metrics


def detect_inflammatory_clusters(
    points_df: pd.DataFrame,
    eps: int = 100,
    min_samples: int = 30,
    min_cluster_size: int = 30,
) -> tuple[int, list]:
    """
    Find clusters of inflammatory cells (potential TLS).

    Uses DBSCAN for efficient clustering (~350x faster than AgglomerativeClustering).

    Parameters
    ----------
    points_df : pd.DataFrame
        DataFrame with x, y, classification_name
    eps : int
        Maximum distance between samples in a cluster (DBSCAN epsilon)
    min_samples : int
        Minimum samples in neighborhood for core point
    min_cluster_size : int
        Minimum cells for valid cluster

    Returns
    -------
    tuple
        (n_valid_clusters, cluster_info_list)
    """
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    if len(inflammatory) < min_cluster_size:
        return 0, []

    coords = inflammatory[["x", "y"]].values

    # Use DBSCAN for efficient density-based clustering
    clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
    labels = clustering.fit_predict(coords)

    # Count valid clusters and get centroids
    cluster_info = []
    unique_labels = set(labels)
    unique_labels.discard(-1)  # Remove noise label

    for label in unique_labels:
        cluster_mask = labels == label
        cluster_size = np.sum(cluster_mask)

        if cluster_size >= min_cluster_size:
            cluster_coords = coords[cluster_mask]
            centroid = cluster_coords.mean(axis=0)
            cluster_info.append({
                "label": int(label),
                "size": int(cluster_size),
                "centroid_x": float(centroid[0]),
                "centroid_y": float(centroid[1]),
            })

    return len(cluster_info), cluster_info


def compute_cohort_thresholds(results: list[dict]) -> dict:
    """
    Compute percentile-based thresholds from cohort distribution.

    Uses tertile-based cutoffs for robust classification that adapts
    to cohort characteristics rather than arbitrary absolute values.

    Parameters
    ----------
    results : list[dict]
        List of metrics dicts from all samples

    Returns
    -------
    dict
        Threshold values for classification
    """
    import numpy as np

    inf_pcts = [r["inflammatory_pct_overall"] for r in results]
    infiltrations = [r["infiltration_ratio"] for r in results]

    return {
        "inflammatory_low": np.percentile(inf_pcts, CONFIG["inflammatory_low_percentile"]),
        "inflammatory_high": np.percentile(inf_pcts, CONFIG["inflammatory_high_percentile"]),
        "infiltration_low": np.percentile(infiltrations, CONFIG["infiltration_low_percentile"]),
    }


def classify_phenotype(metrics: dict, thresholds: dict) -> str:
    """
    Classify sample into immune phenotype based on spatial metrics.

    Uses a simplified 2-axis classification:
    - Axis 1: Inflammatory content (low/medium/high based on cohort tertiles)
    - Axis 2: Infiltration ratio (tumor core vs stroma density)

    Classification logic:
    - Immune-desert: Low inflammatory + poor infiltration
    - Altered-excluded: Low/medium inflammatory + poor infiltration (immune cells blocked from tumor)
    - Altered-immunosuppressed: Medium inflammatory + moderate infiltration (present but not active)
    - Hot: High inflammatory + good infiltration (active anti-tumor response)

    Parameters
    ----------
    metrics : dict
        Computed metrics from compute_regional_metrics
    thresholds : dict
        Cohort-derived percentile thresholds

    Returns
    -------
    str
        Phenotype classification
    """
    inf_pct = metrics["inflammatory_pct_overall"]
    infiltration = metrics["infiltration_ratio"]

    inf_low = thresholds["inflammatory_low"]
    inf_high = thresholds["inflammatory_high"]
    infil_low = thresholds["infiltration_low"]

    # Low inflammatory content (bottom tertile)
    if inf_pct < inf_low:
        if infiltration < infil_low:
            return "Immune-desert"
        else:
            return "Altered-excluded"

    # High inflammatory content (top tertile) with good infiltration = Hot
    elif inf_pct >= inf_high and infiltration >= infil_low:
        return "Hot"

    # Everything else: check infiltration pattern
    else:
        if infiltration < infil_low:
            return "Altered-excluded"
        else:
            return "Altered-immunosuppressed"


def process_sample(zarr_path: Path, sample_id: str) -> dict:
    """
    Process a single sample: load data and compute spatial metrics.

    Note: Phenotype classification is done separately after all samples
    are processed, using cohort-derived percentile thresholds.

    Parameters
    ----------
    zarr_path : Path
        Path to the .zarr directory
    sample_id : str
        Sample identifier

    Returns
    -------
    dict
        Metrics dict (phenotype added later by classify_phenotype)
    """
    # Load spatial data
    points_df = load_sample_spatial_data(zarr_path)

    # Compute region grid
    region_grid, x_bins, y_bins, histograms = compute_region_grid(
        points_df,
        n_bins=CONFIG["n_bins"],
        min_cells_per_bin=CONFIG["min_cells_per_bin"],
    )

    # Compute regional metrics
    metrics = compute_regional_metrics(points_df, region_grid, x_bins, y_bins)

    # Detect inflammatory clusters (supplementary metric, not used for classification)
    n_clusters, cluster_info = detect_inflammatory_clusters(
        points_df,
        eps=CONFIG["cluster_eps"],
        min_samples=CONFIG["cluster_min_samples"],
        min_cluster_size=CONFIG["cluster_min_size"],
    )
    metrics["n_inflammatory_clusters"] = n_clusters
    metrics["sample_id"] = sample_id

    return metrics


def find_available_zarrs() -> list[tuple[str, Path]]:
    """
    Find all available SpatialData zarrs.

    Returns
    -------
    list[tuple[str, Path]]
        List of (sample_id, zarr_path) tuples
    """
    sample_data_dir = SAMPLE_DATA_DIR

    available = []
    if sample_data_dir.exists():
        for sample_dir in sample_data_dir.iterdir():
            if sample_dir.is_dir() and sample_dir.name.startswith("dg_4DFC_"):
                zarr_path = sample_dir / f"{sample_dir.name}_consep_spatialdata.zarr"
                if zarr_path.exists():
                    available.append((sample_dir.name, zarr_path))

    return sorted(available)


def plot_phenotype_distribution(results_df: pd.DataFrame, output_path: Path):
    """Create bar chart of phenotype distribution."""
    phenotype_order = [
        "Immune-desert",
        "Altered-excluded",
        "Altered-immunosuppressed",
        "Hot",
    ]
    phenotype_colors = {
        "Immune-desert": "#2166ac",
        "Altered-excluded": "#92c5de",
        "Altered-immunosuppressed": "#f4a582",
        "Hot": "#b2182b",
    }

    counts = results_df["phenotype"].value_counts()
    counts = counts.reindex(phenotype_order, fill_value=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(
        counts.index,
        counts.values,
        color=[phenotype_colors[p] for p in counts.index],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add count labels
    for bar, count in zip(bars, counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            str(count),
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    ax.set_xlabel("Immune Phenotype", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.set_title("Tumor Immune Phenotype Distribution\nHTAN CRC Cohort", fontsize=14)
    ax.set_ylim(0, max(counts.values) * 1.2)

    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_phenotype_scatter(results_df: pd.DataFrame, output_path: Path, thresholds: dict):
    """Create scatter plot of infiltration ratio vs inflammatory % colored by phenotype."""
    phenotype_colors = {
        "Immune-desert": "#2166ac",
        "Altered-excluded": "#92c5de",
        "Altered-immunosuppressed": "#f4a582",
        "Hot": "#b2182b",
    }

    fig, ax = plt.subplots(figsize=(10, 7))

    for phenotype in phenotype_colors:
        mask = results_df["phenotype"] == phenotype
        if mask.sum() > 0:
            subset = results_df[mask]
            # Cap infiltration ratio for visualization
            infiltration = subset["infiltration_ratio"].clip(upper=5)
            ax.scatter(
                subset["inflammatory_pct_overall"],
                infiltration,
                c=phenotype_colors[phenotype],
                label=f"{phenotype} (n={mask.sum()})",
                s=80,
                alpha=0.7,
                edgecolor="black",
                linewidth=0.5,
            )

    # Show percentile-based thresholds
    ax.axhline(
        y=thresholds["infiltration_low"],
        color="gray",
        linestyle="--",
        alpha=0.5,
        label=f"Infiltration P{CONFIG['infiltration_low_percentile']} ({thresholds['infiltration_low']:.2f})",
    )
    ax.axvline(
        x=thresholds["inflammatory_low"],
        color="gray",
        linestyle=":",
        alpha=0.5,
        label=f"Inflammatory P{CONFIG['inflammatory_low_percentile']} ({thresholds['inflammatory_low']:.1f}%)",
    )
    ax.axvline(
        x=thresholds["inflammatory_high"],
        color="gray",
        linestyle="-.",
        alpha=0.5,
        label=f"Inflammatory P{CONFIG['inflammatory_high_percentile']} ({thresholds['inflammatory_high']:.1f}%)",
    )

    ax.set_xlabel("Overall Inflammatory Cell %", fontsize=12)
    ax.set_ylabel("Infiltration Ratio (core/stroma density)", fontsize=12)
    ax.set_title(
        "Immune Phenotype Classification\nPercentile-Based Thresholds (Cohort-Relative)",
        fontsize=14,
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0, None)
    ax.set_ylim(0, 5.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    """Main entry point for immune phenotype classification."""
    project_root = get_project_root()
    output_dir = project_root / "data"
    figures_dir = project_root / "figures"
    output_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Tumor Immune Phenotype Classification")
    print("=" * 70)

    # Find available zarrs
    available_zarrs = find_available_zarrs()
    print(f"\nFound {len(available_zarrs)} samples with SpatialData zarrs")

    if not available_zarrs:
        print("\nNo SpatialData zarrs found in sample_data/")
        print("Run scripts/fetch_all_spatialdata.sh to download data first.")
        sys.exit(1)

    # Load sample metadata for participant mapping
    print("\nLoading sample metadata...")
    consep_df = load_cell_densities("consep")
    sample_meta = get_sample_metadata(consep_df)
    sample_to_participant = dict(
        zip(sample_meta["sample_id"], sample_meta["participant_id"])
    )
    sample_to_biospecimen = dict(
        zip(sample_meta["sample_id"], sample_meta["biospecimen_id"])
    )

    # Process each sample (compute metrics only)
    print(f"\nProcessing {len(available_zarrs)} samples...")
    results = []

    for sample_id, zarr_path in tqdm(available_zarrs, desc="Computing metrics"):
        try:
            metrics = process_sample(zarr_path, sample_id)
            metrics["participant_id"] = sample_to_participant.get(sample_id, "")
            metrics["biospecimen_id"] = sample_to_biospecimen.get(sample_id, "")
            results.append(metrics)
        except Exception as e:
            print(f"\n  Error processing {sample_id}: {e}")
            continue

    if not results:
        print("\nNo samples were successfully processed.")
        sys.exit(1)

    # Compute cohort-based thresholds
    print("\nComputing cohort-based thresholds...")
    thresholds = compute_cohort_thresholds(results)
    print(f"  Inflammatory low (P{CONFIG['inflammatory_low_percentile']}): {thresholds['inflammatory_low']:.1f}%")
    print(f"  Inflammatory high (P{CONFIG['inflammatory_high_percentile']}): {thresholds['inflammatory_high']:.1f}%")
    print(f"  Infiltration low (P{CONFIG['infiltration_low_percentile']}): {thresholds['infiltration_low']:.2f}")

    # Classify each sample using cohort thresholds
    print("\nClassifying samples...")
    for metrics in results:
        metrics["phenotype"] = classify_phenotype(metrics, thresholds)

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Reorder columns for clarity
    key_cols = [
        "sample_id",
        "participant_id",
        "biospecimen_id",
        "phenotype",
        "inflammatory_pct_overall",
        "infiltration_ratio",
        "n_inflammatory_clusters",
    ]
    other_cols = [c for c in results_df.columns if c not in key_cols]
    results_df = results_df[key_cols + other_cols]

    # Sort by inflammatory percentage
    results_df = results_df.sort_values("inflammatory_pct_overall", ascending=False)

    # Save full metrics
    metrics_path = output_dir / "immune_phenotype_metrics.csv"
    results_df.to_csv(metrics_path, index=False)
    print(f"\nSaved full metrics to: {metrics_path}")

    # Save classification summary
    classification_cols = [
        "sample_id",
        "participant_id",
        "phenotype",
        "inflammatory_pct_overall",
        "infiltration_ratio",
        "n_inflammatory_clusters",
    ]
    classification_df = results_df[classification_cols].copy()
    classification_path = output_dir / "immune_phenotype_classification.csv"
    classification_df.to_csv(classification_path, index=False)
    print(f"Saved classification to: {classification_path}")

    # Generate visualizations
    print("\nGenerating visualizations...")

    bar_path = figures_dir / "immune_phenotype_distribution.png"
    plot_phenotype_distribution(results_df, bar_path)
    print(f"  Saved: {bar_path}")

    scatter_path = figures_dir / "immune_phenotype_scatter.png"
    plot_phenotype_scatter(results_df, scatter_path, thresholds)
    print(f"  Saved: {scatter_path}")

    # Print summary
    print("\n" + "=" * 70)
    print("Classification Summary")
    print("=" * 70)

    phenotype_counts = results_df["phenotype"].value_counts()
    for phenotype in ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]:
        count = phenotype_counts.get(phenotype, 0)
        pct = count / len(results_df) * 100
        print(f"  {phenotype:25s}: {count:3d} ({pct:5.1f}%)")

    print("\n" + "-" * 70)
    print("Sample Classifications (sorted by inflammatory %):")
    print("-" * 70)
    print(f"{'Participant':<12} {'Inflammatory %':>14} {'Infiltration':>12} {'Clusters':>10} {'Phenotype':<25}")
    print("-" * 70)

    for _, row in results_df.iterrows():
        inf_ratio = row["infiltration_ratio"]
        inf_ratio_str = f"{inf_ratio:.2f}" if inf_ratio < 100 else ">100"
        print(
            f"{row['participant_id']:<12} "
            f"{row['inflammatory_pct_overall']:>13.1f}% "
            f"{inf_ratio_str:>12} "
            f"{row['n_inflammatory_clusters']:>10} "
            f"{row['phenotype']:<25}"
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
