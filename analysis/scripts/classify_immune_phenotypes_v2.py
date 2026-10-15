#!/usr/bin/env python3
"""
Tumor Immune Phenotype Classification - Version 2

Uses KNN-based spatial graph analysis:
1. KNN graph construction with sklearn/scipy
2. Neighborhood composition analysis
3. Cell-by-cell microenvironment classification

Key improvements over v1:
- Cell-graph based region definition (not grid binning)
- Organic Tumor-Stroma Interface (TSI) from neighborhood composition
- Absolute density thresholds (cells/mm²) from immunoscore literature

Reference thresholds:
- Galon et al. (2006, 2014): CD3/CD8 densities in CRC
- Pages et al. (2018): International validation of Immunoscore
"""

from pathlib import Path
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import get_project_root, load_cell_densities, get_sample_metadata


# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG = {
    # KNN spatial graph parameters
    "n_neighbors": 30,  # K for KNN graph

    # Pixel to micron conversion (typical 0.25 µm/pixel at 40x)
    "um_per_pixel": 0.25,

    # Absolute density thresholds (cells/mm²)
    # ConSeP "Inflammatory" = all immune cells (~3x CD3/CD8 counts)
    # Scaled from immunoscore: 1500 × 3 = 4500, 500 × 3 = 1500
    "immune_density_high": 6000,  # cells/mm² - "Hot" threshold
    "immune_density_low": 2500,   # cells/mm² - "Desert" threshold

    # Tumor-associated threshold (% of immune cells in tumor core + TSI)
    # >= 20% means good infiltration into tumor regions
    "tumor_associated_good": 20.0,

    # Neighborhood composition thresholds
    "tumor_core_malignant_frac": 0.60,  # >60% malignant neighbors
    "tsi_malignant_min": 0.15,          # TSI: 15-60% malignant
    "stroma_stroma_frac": 0.40,         # >40% stroma neighbors
    "stroma_malignant_max": 0.15,       # <15% malignant in stroma
}

REGION_CODES = {
    "other": 0,
    "tumor_core": 1,
    "tsi": 2,
    "stroma": 3,
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
    2: "TSI",
    3: "Stroma",
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_sample_data(sample_id: str) -> pd.DataFrame:
    """
    Load sample data from parquet.

    Returns DataFrame with x, y, classification_name columns.
    """
    parquet_path = (
        SAMPLE_DATA_DIR / sample_id
        / f"{sample_id}_consep_spatialdata.zarr"
        / "points" / "cell_centroids" / "points.parquet"
    )

    return pd.read_parquet(parquet_path, columns=["x", "y", "classification_name"])


# =============================================================================
# KNN-BASED MICROENVIRONMENT CLASSIFICATION
# =============================================================================

def classify_microenvironments_knn(df: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    """
    Classify cell microenvironments using KNN spatial graphs.

    Steps:
    1. Build KNN graph using cKDTree (faster than sklearn for 2D)
    2. Compute neighborhood cell type fractions for each cell
    3. Classify based on local composition

    Returns:
    - regions: array of region codes for each cell
    - df: DataFrame with added neighborhood fraction columns
    """
    coords = df[["x", "y"]].values
    # Convert to numpy array (required for 2D indexing with neighbor_indices)
    cell_types = np.array(df["classification_name"].tolist())
    n_cells = len(df)

    # Build KNN graph using cKDTree (much faster for 2D spatial data)
    k = CONFIG["n_neighbors"]
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1, workers=-1)

    # Exclude self (first neighbor is self)
    neighbor_indices = indices[:, 1:]

    # Compute neighborhood composition for each cell
    # Vectorized for speed
    neighbor_types = cell_types[neighbor_indices]

    malignant_frac = (neighbor_types == "Malignant epithelial").sum(axis=1) / k
    stroma_frac = (neighbor_types == "Stroma").sum(axis=1) / k
    inflammatory_frac = (neighbor_types == "Inflammatory").sum(axis=1) / k

    # Store in DataFrame
    df = df.copy()
    df["neighbor_malignant_frac"] = malignant_frac
    df["neighbor_stroma_frac"] = stroma_frac
    df["neighbor_inflammatory_frac"] = inflammatory_frac

    # Classify regions
    regions = np.zeros(n_cells, dtype=int)

    # Tumor core: high malignant neighborhood
    tumor_core_mask = malignant_frac >= CONFIG["tumor_core_malignant_frac"]
    regions[tumor_core_mask] = REGION_CODES["tumor_core"]

    # Stroma: high stroma, low malignant
    stroma_mask = (
        (stroma_frac >= CONFIG["stroma_stroma_frac"]) &
        (malignant_frac < CONFIG["stroma_malignant_max"]) &
        ~tumor_core_mask
    )
    regions[stroma_mask] = REGION_CODES["stroma"]

    # TSI: interface region with mixed composition
    tsi_mask = (
        (malignant_frac >= CONFIG["tsi_malignant_min"]) &
        (malignant_frac < CONFIG["tumor_core_malignant_frac"]) &
        ~stroma_mask
    )
    regions[tsi_mask] = REGION_CODES["tsi"]

    df["region"] = regions
    df["region_label"] = [REGION_LABELS[r] for r in regions]

    return regions, df


def compute_region_areas(df: pd.DataFrame, regions: np.ndarray) -> dict:
    """
    Estimate region areas using Voronoi-like cell territories.

    Each cell's territory approximated as area around it based on
    neighbor distances. Returns areas in mm².
    """
    coords = df[["x", "y"]].values

    # Build KD-tree for nearest neighbor distances
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=2)
    nn_distances = distances[:, 1]

    # Territory radius = half NN distance (capped)
    territory_radius = np.clip(nn_distances / 2, 1, 100)
    territory_area_px = np.pi * territory_radius ** 2

    # Convert to mm²
    um_per_px = CONFIG["um_per_pixel"]
    territory_area_mm2 = territory_area_px * (um_per_px / 1000) ** 2

    # Sum by region
    areas = {}
    for region_name, region_code in REGION_CODES.items():
        mask = regions == region_code
        areas[region_name] = territory_area_mm2[mask].sum()

    return areas


# =============================================================================
# METRICS COMPUTATION
# =============================================================================

def compute_spatial_metrics(df: pd.DataFrame, regions: np.ndarray) -> dict:
    """Compute immune density metrics with absolute thresholds (cells/mm²)."""

    # Convert to numpy array to avoid Arrow issues
    cell_types = np.array(df["classification_name"].tolist())
    is_inflammatory = cell_types == "Inflammatory"

    areas = compute_region_areas(df, regions)

    metrics = {
        "total_cells": len(df),
        "inflammatory_count": is_inflammatory.sum(),
    }

    # Per-region metrics
    for region_name, region_code in REGION_CODES.items():
        region_mask = regions == region_code
        n_cells = region_mask.sum()
        n_inflammatory = (region_mask & is_inflammatory).sum()
        area_mm2 = areas[region_name]

        metrics[f"n_cells_{region_name}"] = n_cells
        metrics[f"n_inflammatory_{region_name}"] = n_inflammatory
        metrics[f"area_mm2_{region_name}"] = area_mm2
        metrics[f"density_{region_name}"] = n_inflammatory / area_mm2 if area_mm2 > 0 else 0

    # Overall metrics
    total_area = sum(areas.values())
    metrics["total_area_mm2"] = total_area
    metrics["density_overall"] = metrics["inflammatory_count"] / total_area if total_area > 0 else 0
    metrics["inflammatory_pct"] = metrics["inflammatory_count"] / metrics["total_cells"] * 100

    # Infiltration ratio (density-based, kept for reference)
    core_density = metrics["density_tumor_core"]
    stroma_density = metrics["density_stroma"]
    metrics["infiltration_ratio"] = core_density / stroma_density if stroma_density > 0 else float("inf")

    # Tumor-associated % (count-based - more intuitive)
    # % of immune cells in tumor core + TSI regions
    total_inflammatory = metrics["inflammatory_count"]
    if total_inflammatory > 0:
        tumor_inflammatory = metrics["n_inflammatory_tumor_core"] + metrics["n_inflammatory_tsi"]
        metrics["tumor_associated_pct"] = tumor_inflammatory / total_inflammatory * 100
    else:
        metrics["tumor_associated_pct"] = 0

    return metrics


def compute_colocalization_score(df: pd.DataFrame) -> dict:
    """
    Compute co-localization scores using neighborhood fractions.

    Returns metrics for inflammatory-tumor proximity.
    """
    # Get inflammatory cells (convert to numpy for proper boolean operations)
    cell_types = np.array(df["classification_name"].tolist())
    inflammatory_mask = cell_types == "Inflammatory"

    if inflammatory_mask.sum() == 0:
        return {
            "inflammatory_in_tumor_neighborhood": 0,
            "inflammatory_in_stroma_neighborhood": 0,
        }

    # Average malignant fraction in neighborhoods of inflammatory cells
    # Use .values to get numpy arrays for proper indexing
    avg_malignant_near_inflammatory = df["neighbor_malignant_frac"].values[inflammatory_mask].mean()
    avg_stroma_near_inflammatory = df["neighbor_stroma_frac"].values[inflammatory_mask].mean()

    return {
        "inflammatory_in_tumor_neighborhood": avg_malignant_near_inflammatory,
        "inflammatory_in_stroma_neighborhood": avg_stroma_near_inflammatory,
    }


# =============================================================================
# PHENOTYPE CLASSIFICATION
# =============================================================================

def classify_phenotype(metrics: dict) -> str:
    """
    Classify using density and tumor-associated % thresholds.

    Two axes:
    1. Density (cells/mm²): High (>6000), Medium (2500-6000), Low (<2500)
    2. Tumor infiltration: % of immune cells in tumor core + TSI regions

    Classification:
    - Hot: High density + good infiltration (>=20% tumor-associated)
    - Altered-immunosuppressed: Medium density + good infiltration
    - Altered-excluded: Any density + poor infiltration (<20% tumor-associated)
    - Immune-desert: Low density (regardless of infiltration)
    """
    density = metrics["density_overall"]
    tumor_assoc = metrics["tumor_associated_pct"]

    high_thresh = CONFIG["immune_density_high"]
    low_thresh = CONFIG["immune_density_low"]
    infil_thresh = CONFIG["tumor_associated_good"]

    if density < low_thresh:
        return "Immune-desert"
    elif density >= high_thresh:
        return "Hot" if tumor_assoc >= infil_thresh else "Altered-excluded"
    else:
        return "Altered-immunosuppressed" if tumor_assoc >= infil_thresh else "Altered-excluded"


# =============================================================================
# SAMPLE PROCESSING
# =============================================================================

def process_sample(sample_id: str) -> tuple[dict, pd.DataFrame]:
    """Process single sample with KNN analysis."""

    df = load_sample_data(sample_id)
    regions, df = classify_microenvironments_knn(df)
    metrics = compute_spatial_metrics(df, regions)

    # Add co-localization scores
    coloc = compute_colocalization_score(df)
    metrics.update(coloc)

    metrics["phenotype"] = classify_phenotype(metrics)
    metrics["sample_id"] = sample_id

    return metrics, df


def find_available_zarrs() -> list[tuple[str, Path]]:
    """Find available SpatialData zarrs."""
    sample_data_dir = SAMPLE_DATA_DIR

    available = []
    if sample_data_dir.exists():
        for sample_dir in sample_data_dir.iterdir():
            if sample_dir.is_dir() and sample_dir.name.startswith("dg_4DFC_"):
                zarr_path = sample_dir / f"{sample_dir.name}_consep_spatialdata.zarr"
                if zarr_path.exists():
                    available.append((sample_dir.name, zarr_path))

    return sorted(available)


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_phenotype_distribution(results_df: pd.DataFrame, output_path: Path):
    """Bar chart of phenotype distribution."""
    phenotype_order = ["Immune-desert", "Altered-excluded", "Altered-immunosuppressed", "Hot"]
    colors = {"Immune-desert": "#2166ac", "Altered-excluded": "#92c5de",
              "Altered-immunosuppressed": "#f4a582", "Hot": "#b2182b"}

    counts = results_df["phenotype"].value_counts().reindex(phenotype_order, fill_value=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(counts.index, counts.values,
                  color=[colors[p] for p in counts.index],
                  edgecolor="black", linewidth=0.5)

    for bar, count in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                str(count), ha="center", va="bottom", fontsize=12, fontweight="bold")

    ax.set_xlabel("Immune Phenotype", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.set_title("Tumor Immune Phenotype Distribution (v2)\nSquidpy KNN graphs + Absolute thresholds", fontsize=13)
    ax.set_ylim(0, max(counts.values) * 1.2 if counts.max() > 0 else 10)
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_density_scatter(results_df: pd.DataFrame, output_path: Path):
    """Scatter plot: density vs tumor-associated %."""
    colors = {"Immune-desert": "#2166ac", "Altered-excluded": "#92c5de",
              "Altered-immunosuppressed": "#f4a582", "Hot": "#b2182b"}

    fig, ax = plt.subplots(figsize=(10, 7))

    for phenotype, color in colors.items():
        mask = results_df["phenotype"] == phenotype
        if mask.sum() > 0:
            subset = results_df[mask]
            ax.scatter(subset["density_overall"],
                       subset["tumor_associated_pct"],
                       c=color, label=f"{phenotype} (n={mask.sum()})",
                       s=80, alpha=0.7, edgecolor="black", linewidth=0.5)

    # Threshold lines
    ax.axhline(y=CONFIG["tumor_associated_good"], color="gray", linestyle="--", alpha=0.7,
               label=f"Infiltration threshold ({CONFIG['tumor_associated_good']}%)")
    ax.axvline(x=CONFIG["immune_density_low"], color="gray", linestyle=":", alpha=0.7)
    ax.axvline(x=CONFIG["immune_density_high"], color="gray", linestyle="-.", alpha=0.7)

    ax.set_xlabel("Inflammatory Density (cells/mm²)", fontsize=12)
    ax.set_ylabel("% Immune Cells in Tumor Regions (Core + TSI)", fontsize=12)
    ax.set_title("Immune Phenotype Classification (v2)\nDensity vs Tumor Infiltration", fontsize=14)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0, None)
    ax.set_ylim(0, 50)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    project_root = get_project_root()
    output_dir = project_root / "data"
    figures_dir = project_root / "figures"

    print("=" * 70)
    print("Tumor Immune Phenotype Classification - Version 2")
    print("KNN spatial graphs + Absolute density thresholds")
    print("=" * 70)

    print(f"\nConfiguration:")
    print(f"  K neighbors: {CONFIG['n_neighbors']}")
    print(f"  High density: {CONFIG['immune_density_high']} cells/mm²")
    print(f"  Low density: {CONFIG['immune_density_low']} cells/mm²")
    print(f"  Tumor-associated threshold: {CONFIG['tumor_associated_good']}%")

    available_zarrs = find_available_zarrs()
    print(f"\nFound {len(available_zarrs)} samples")

    if not available_zarrs:
        print("No samples found!")
        sys.exit(1)

    # Load metadata
    consep_df = load_cell_densities("consep")
    sample_meta = get_sample_metadata(consep_df)
    sample_to_participant = dict(zip(sample_meta["sample_id"], sample_meta["participant_id"]))

    # Process samples
    print("\nProcessing samples...")
    results = []

    for sample_id, _ in tqdm(available_zarrs, desc="Classifying"):
        try:
            metrics, df = process_sample(sample_id)
            metrics["participant_id"] = sample_to_participant.get(sample_id, "")
            results.append(metrics)
        except Exception as e:
            print(f"\n  Error: {sample_id}: {e}")
            continue

    results_df = pd.DataFrame(results)

    # Reorder and sort
    key_cols = ["sample_id", "participant_id", "phenotype", "density_overall",
                "tumor_associated_pct", "infiltration_ratio", "density_tumor_core", "density_tsi", "density_stroma"]
    other_cols = [c for c in results_df.columns if c not in key_cols]
    results_df = results_df[[c for c in key_cols if c in results_df.columns] + other_cols]
    results_df = results_df.sort_values("density_overall", ascending=False)

    # Save
    results_df.to_csv(output_dir / "immune_phenotype_metrics_v2.csv", index=False)
    results_df[["sample_id", "participant_id", "phenotype", "density_overall", "tumor_associated_pct"]].to_csv(
        output_dir / "immune_phenotype_classification_v2.csv", index=False)

    # Visualizations
    print("\nGenerating visualizations...")
    plot_phenotype_distribution(results_df, figures_dir / "immune_phenotype_distribution_v2.png")
    plot_density_scatter(results_df, figures_dir / "immune_phenotype_scatter_v2.png")

    # Summary
    print("\n" + "=" * 70)
    print("Classification Summary (v2)")
    print("=" * 70)

    for phenotype in ["Hot", "Altered-immunosuppressed", "Altered-excluded", "Immune-desert"]:
        count = (results_df["phenotype"] == phenotype).sum()
        print(f"  {phenotype:25s}: {count:3d} ({count/len(results_df)*100:5.1f}%)")

    print("\n" + "-" * 75)
    print(f"{'Participant':<12} {'Density':>12} {'Tumor-Assoc%':>12} {'Phenotype':<25}")
    print("-" * 75)

    for _, row in results_df.iterrows():
        print(f"{row['participant_id']:<12} {row['density_overall']:>10.0f}/mm² {row['tumor_associated_pct']:>11.1f}% {row['phenotype']:<25}")

    print("\nDone!")


if __name__ == "__main__":
    main()
