#!/usr/bin/env python3
"""
Compute truly spatial features from cell coordinates for survival analysis.

Features computed:
1. Immune-tumor nearest-neighbor distances — continuous measure of immune
   penetration using cKDTree on malignant vs inflammatory cell coordinates
2. DBSCAN cluster features — density-based clustering with fixed eps
   (32.5um neighborhood) on all and tumor-proximal inflammatory cells
3. HDBSCAN cluster features — hierarchical density-based clustering with
   no eps parameter, on all and tumor-proximal inflammatory cells
4. Tumor-immune Moran's I — spatial autocorrelation of immune cell density
   within the tumor region

Usage:
    python scripts/compute_spatial_features.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from fast_hdbscan import HDBSCAN
from sklearn.cluster import DBSCAN
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)

sys.path.insert(0, str(Path(__file__).parent))
from classify_immune_phenotypes import find_available_zarrs, load_sample_spatial_data
from load_cell_densities import get_project_root, load_cell_densities, get_sample_metadata
from utils.spatial_utils import UM_PER_PIXEL

# Clustering parameters
MIN_CLUSTER_SIZE = 30  # minimum immune aggregate size (in cells)
DBSCAN_EPS = 100  # pixels (~32.5um, ~2-3 cell diameters)
DBSCAN_MIN_SAMPLES = 30


def compute_nn_distances(points_df: pd.DataFrame) -> dict:
    """
    Compute nearest-neighbor distances from malignant cells to inflammatory cells.

    Parameters
    ----------
    points_df : pd.DataFrame
        DataFrame with x, y, classification_name columns

    Returns
    -------
    dict
        nn_dist_median, nn_dist_mean, nn_dist_10pct (all in microns)
    """
    malignant = points_df[points_df["classification_name"] == "Malignant epithelial"]
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    if len(malignant) == 0 or len(inflammatory) == 0:
        return {
            "nn_dist_median": np.nan,
            "nn_dist_mean": np.nan,
            "nn_dist_10pct": np.nan,
        }

    # Build KD-tree on inflammatory cell coordinates
    tree = cKDTree(inflammatory[["x", "y"]].values)

    # Query nearest neighbor for every malignant cell
    distances, _ = tree.query(malignant[["x", "y"]].values, k=1)

    # Convert pixels to microns
    distances_um = distances * UM_PER_PIXEL

    return {
        "nn_dist_median": float(np.median(distances_um)),
        "nn_dist_mean": float(np.mean(distances_um)),
        "nn_dist_10pct": float(np.percentile(distances_um, 10)),
    }


def _run_dbscan(coords: np.ndarray, prefix: str = "db_") -> dict:
    """Run DBSCAN on cell coordinates and return cluster summary features."""
    empty = {
        f"{prefix}n_clusters": 0,
        f"{prefix}cluster_fraction": 0.0,
        f"{prefix}mean_cluster_size": 0.0,
        f"{prefix}largest_cluster_size": 0,
    }

    if len(coords) < MIN_CLUSTER_SIZE:
        return empty

    labels = DBSCAN(
        eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES, n_jobs=-1
    ).fit_predict(coords)

    unique_labels = set(labels)
    unique_labels.discard(-1)

    cluster_sizes = [
        int(np.sum(labels == label))
        for label in unique_labels
        if np.sum(labels == label) >= MIN_CLUSTER_SIZE
    ]

    if len(cluster_sizes) == 0:
        return empty

    clustered_cells = sum(cluster_sizes)
    return {
        f"{prefix}n_clusters": len(cluster_sizes),
        f"{prefix}cluster_fraction": float(clustered_cells / len(coords)),
        f"{prefix}mean_cluster_size": float(np.mean(cluster_sizes)),
        f"{prefix}largest_cluster_size": max(cluster_sizes),
    }


def compute_dbscan_features(points_df: pd.DataFrame) -> dict:
    """DBSCAN cluster features for ALL inflammatory cells."""
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]
    return _run_dbscan(inflammatory[["x", "y"]].values)


def compute_tumor_dbscan_features(
    points_df: pd.DataFrame,
    distance_threshold_um: float = 100,
) -> dict:
    """DBSCAN cluster features for tumor-proximal inflammatory cells only."""
    malignant = points_df[points_df["classification_name"] == "Malignant epithelial"]
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    if len(malignant) < 10 or len(inflammatory) < 10:
        return _run_dbscan(np.empty((0, 2)), prefix="db_tumor_")

    mal_tree = cKDTree(malignant[["x", "y"]].values)
    inf_coords = inflammatory[["x", "y"]].values
    dists, _ = mal_tree.query(inf_coords, k=1)

    distance_threshold_px = distance_threshold_um / UM_PER_PIXEL
    proximal_coords = inf_coords[dists <= distance_threshold_px]

    return _run_dbscan(proximal_coords, prefix="db_tumor_")


def _run_hdbscan(coords: np.ndarray, prefix: str = "") -> dict:
    """
    Run HDBSCAN on cell coordinates and return cluster summary features.

    HDBSCAN builds a cluster hierarchy across all density scales and extracts
    the most stable clusters. Unlike DBSCAN, there is no eps parameter — only
    min_cluster_size, which maps directly to the biological question of minimum
    immune aggregate size.

    Parameters
    ----------
    coords : np.ndarray
        (N, 2) array of cell coordinates in pixels.
    prefix : str
        Prefix for feature names (e.g. "tumor_" for tumor-proximal features).

    Returns
    -------
    dict
        {prefix}n_clusters, {prefix}cluster_fraction,
        {prefix}mean_cluster_size, {prefix}largest_cluster_size
    """
    empty = {
        f"{prefix}n_clusters": 0,
        f"{prefix}cluster_fraction": 0.0,
        f"{prefix}mean_cluster_size": 0.0,
        f"{prefix}largest_cluster_size": 0,
    }

    if len(coords) < MIN_CLUSTER_SIZE:
        return empty

    labels = HDBSCAN(min_cluster_size=MIN_CLUSTER_SIZE).fit_predict(coords)

    unique_labels = set(labels)
    unique_labels.discard(-1)

    cluster_sizes = [int(np.sum(labels == label)) for label in unique_labels]

    if len(cluster_sizes) == 0:
        return empty

    clustered_cells = sum(cluster_sizes)
    return {
        f"{prefix}n_clusters": len(cluster_sizes),
        f"{prefix}cluster_fraction": float(clustered_cells / len(coords)),
        f"{prefix}mean_cluster_size": float(np.mean(cluster_sizes)),
        f"{prefix}largest_cluster_size": max(cluster_sizes),
    }


def compute_cluster_features(points_df: pd.DataFrame) -> dict:
    """HDBSCAN cluster features for ALL inflammatory cells."""
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]
    return _run_hdbscan(inflammatory[["x", "y"]].values)


def compute_tumor_cluster_features(
    points_df: pd.DataFrame,
    distance_threshold_um: float = 100,
) -> dict:
    """
    HDBSCAN cluster features for tumor-proximal inflammatory cells only.

    Excludes immune aggregates in healthy tissue (e.g. Peyer's patches) by
    restricting to inflammatory cells within distance_threshold_um of the
    nearest malignant cell.
    """
    malignant = points_df[points_df["classification_name"] == "Malignant epithelial"]
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    if len(malignant) < 10 or len(inflammatory) < 10:
        return _run_hdbscan(np.empty((0, 2)), prefix="tumor_")

    mal_tree = cKDTree(malignant[["x", "y"]].values)
    inf_coords = inflammatory[["x", "y"]].values
    dists, _ = mal_tree.query(inf_coords, k=1)

    distance_threshold_px = distance_threshold_um / UM_PER_PIXEL
    proximal_coords = inf_coords[dists <= distance_threshold_px]

    return _run_hdbscan(proximal_coords, prefix="tumor_")


def _morans_i_masked_grid(values: np.ndarray, grid_positions: list[tuple[int, int]]) -> float:
    """
    Compute Moran's I with queen contiguity on a subset of grid cells.

    Parameters
    ----------
    values : np.ndarray
        1D array of values at each valid grid cell (e.g. immune cell counts).
    grid_positions : list of (row, col) tuples
        Grid coordinates for each value entry.

    Returns
    -------
    float
        Moran's I statistic, or NaN if fewer than 10 valid bins.
    """
    n = len(values)
    if n < 10:
        return np.nan

    mean = values.mean()
    deviations = values - mean
    ss = np.sum(deviations ** 2)
    if ss == 0:
        return np.nan

    # Build lookup from grid position to index
    pos_to_idx = {pos: i for i, pos in enumerate(grid_positions)}

    # Queen contiguity: 8 neighbors
    offsets = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    weighted_sum = 0.0
    total_weights = 0
    for i, (r, c) in enumerate(grid_positions):
        for dr, dc in offsets:
            j = pos_to_idx.get((r + dr, c + dc))
            if j is not None:
                weighted_sum += deviations[i] * deviations[j]
                total_weights += 1

    if total_weights == 0:
        return np.nan

    return (n / total_weights) * (weighted_sum / ss)


def compute_tumor_immune_morans_i(
    points_df: pd.DataFrame,
    distance_threshold_um: float = 100,
    n_bins: int = 50,
) -> dict:
    """
    Compute Moran's I for tumor-proximal immune cell density.

    Identifies inflammatory cells within `distance_threshold_um` of the nearest
    malignant cell, grids the tumor region, counts proximal inflammatory cells
    per bin, and computes Moran's I with queen contiguity on the tumor-region
    bins.

    Parameters
    ----------
    points_df : pd.DataFrame
        DataFrame with x, y, classification_name columns.
    distance_threshold_um : float
        Maximum distance (microns) from nearest malignant cell to be considered
        tumor-proximal. Default 100um (~308px, ~7-10 cell diameters).
    n_bins : int
        Number of bins per dimension for the grid.

    Returns
    -------
    dict
        {"tumor_immune_morans_i": float}
    """
    malignant = points_df[points_df["classification_name"] == "Malignant epithelial"]
    inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

    if len(malignant) < 10 or len(inflammatory) < 10:
        return {"tumor_immune_morans_i": np.nan}

    mal_coords = malignant[["x", "y"]].values
    inf_coords = inflammatory[["x", "y"]].values

    # Find tumor-proximal inflammatory cells
    distance_threshold_px = distance_threshold_um / UM_PER_PIXEL
    mal_tree = cKDTree(mal_coords)
    dists, _ = mal_tree.query(inf_coords, k=1)
    proximal_mask = dists <= distance_threshold_px
    proximal_coords = inf_coords[proximal_mask]

    if len(proximal_coords) < 10:
        return {"tumor_immune_morans_i": np.nan}

    # Grid the tumor region (bounding box of malignant cells + buffer)
    buffer = distance_threshold_px
    x_min = mal_coords[:, 0].min() - buffer
    x_max = mal_coords[:, 0].max() + buffer
    y_min = mal_coords[:, 1].min() - buffer
    y_max = mal_coords[:, 1].max() + buffer

    x_edges = np.linspace(x_min, x_max, n_bins + 1)
    y_edges = np.linspace(y_min, y_max, n_bins + 1)

    # Bin malignant cells to identify valid (tumor-region) bins
    mal_x_bin = np.clip(np.digitize(mal_coords[:, 0], x_edges) - 1, 0, n_bins - 1)
    mal_y_bin = np.clip(np.digitize(mal_coords[:, 1], y_edges) - 1, 0, n_bins - 1)
    valid_bins = set(zip(mal_y_bin, mal_x_bin))  # (row, col) format

    if len(valid_bins) < 10:
        return {"tumor_immune_morans_i": np.nan}

    # Count proximal inflammatory cells per bin
    inf_x_bin = np.clip(np.digitize(proximal_coords[:, 0], x_edges) - 1, 0, n_bins - 1)
    inf_y_bin = np.clip(np.digitize(proximal_coords[:, 1], y_edges) - 1, 0, n_bins - 1)

    # Only count in valid (tumor-region) bins
    counts = {}
    for pos in valid_bins:
        counts[pos] = 0
    for r, c in zip(inf_y_bin, inf_x_bin):
        pos = (r, c)
        if pos in counts:
            counts[pos] += 1

    grid_positions = sorted(counts.keys())
    values = np.array([counts[pos] for pos in grid_positions], dtype=float)

    morans_i = _morans_i_masked_grid(values, grid_positions)
    return {"tumor_immune_morans_i": float(morans_i) if not np.isnan(morans_i) else np.nan}


def process_sample(sample_id: str, zarr_path: Path) -> dict:
    """
    Compute all spatial features for one sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier (dg_4DFC_UUID)
    zarr_path : Path
        Path to the .zarr directory

    Returns
    -------
    dict
        All spatial features plus sample_id
    """
    points_df = load_sample_spatial_data(zarr_path)

    features = {"sample_id": sample_id}
    features.update(compute_nn_distances(points_df))
    features.update(compute_dbscan_features(points_df))
    features.update(compute_tumor_dbscan_features(points_df))
    features.update(compute_cluster_features(points_df))
    features.update(compute_tumor_cluster_features(points_df))
    features.update(compute_tumor_immune_morans_i(points_df))

    return features


def main():
    project_root = get_project_root()
    data_dir = project_root / "data"
    data_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Computing Spatial Features (NN distances + DBSCAN + HDBSCAN + Moran's I)")
    print("=" * 70)

    # Find available zarrs
    zarrs = find_available_zarrs()
    print(f"\nFound {len(zarrs)} samples")

    # Process each sample
    results = []
    for sample_id, zarr_path in tqdm(zarrs, desc="Processing samples"):
        features = process_sample(sample_id, zarr_path)
        results.append(features)

    results_df = pd.DataFrame(results)

    # Get sample → participant mapping
    consep = load_cell_densities("consep")
    metadata = get_sample_metadata(consep)
    results_df = results_df.merge(metadata[["sample_id", "participant_id"]], on="sample_id", how="left")

    # Aggregate to participant level (mean for HTA7_919 which has 2 samples)
    feature_cols = [
        "nn_dist_median", "nn_dist_mean", "nn_dist_10pct",
        "db_n_clusters", "db_cluster_fraction", "db_mean_cluster_size", "db_largest_cluster_size",
        "db_tumor_n_clusters", "db_tumor_cluster_fraction", "db_tumor_mean_cluster_size", "db_tumor_largest_cluster_size",
        "n_clusters", "cluster_fraction", "mean_cluster_size", "largest_cluster_size",
        "tumor_n_clusters", "tumor_cluster_fraction", "tumor_mean_cluster_size", "tumor_largest_cluster_size",
        "tumor_immune_morans_i",
    ]
    participant_df = results_df.groupby("participant_id")[feature_cols].mean().reset_index()

    # Save
    output_path = data_dir / "spatial_features.csv"
    participant_df.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path} ({len(participant_df)} participants)")

    # Print summary
    print("\nFeature ranges:")
    for col in feature_cols:
        vals = participant_df[col].dropna()
        print(f"  {col:25s}: {vals.min():.2f} — {vals.max():.2f} (median {vals.median():.2f})")

    print(f"\nNaN check:")
    for col in feature_cols:
        n_nan = participant_df[col].isna().sum()
        if n_nan > 0:
            print(f"  {col}: {n_nan} NaN values")
    if participant_df[feature_cols].isna().sum().sum() == 0:
        print("  No NaN values in any feature column")


if __name__ == "__main__":
    main()
