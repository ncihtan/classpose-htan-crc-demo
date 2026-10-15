#!/usr/bin/env python3
"""
DBSCAN parameter sensitivity analysis for immune clustering features.

Tests whether the prognostic signal from immune cluster_fraction is robust
across a range of DBSCAN parameters (eps, min_samples), or whether it depends
on a specific parameter choice. Runs both all-cell and tumor-proximal DBSCAN
variants.

Usage:
    python scripts/dbscan_sensitivity.py
"""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*delta method.*")

sys.path.insert(0, str(Path(__file__).parent))
from classify_immune_phenotypes import find_available_zarrs, load_sample_spatial_data
from load_cell_densities import get_project_root, load_cell_densities, get_sample_metadata
from utils.spatial_utils import UM_PER_PIXEL

# Parameter grid
EPS_VALUES = [75, 100, 125, 150]
MIN_SAMPLES_VALUES = [20, 30, 50]

# Tumor-proximity threshold
DISTANCE_THRESHOLD_UM = 100
MIN_CLUSTER_SIZE = 30


def load_all_samples():
    """Load spatial data and pre-compute tumor-proximal coords for all samples."""
    zarrs = find_available_zarrs()
    print(f"Found {len(zarrs)} samples")

    samples = []
    for sample_id, zarr_path in tqdm(zarrs, desc="Loading samples"):
        points_df = load_sample_spatial_data(zarr_path)

        malignant = points_df[points_df["classification_name"] == "Malignant epithelial"]
        inflammatory = points_df[points_df["classification_name"] == "Inflammatory"]

        all_inf_coords = inflammatory[["x", "y"]].values

        # Compute tumor-proximal inflammatory coords
        if len(malignant) >= 10 and len(inflammatory) >= 10:
            mal_coords = malignant[["x", "y"]].values
            mal_tree = cKDTree(mal_coords)
            dists, _ = mal_tree.query(all_inf_coords, k=1)
            threshold_px = DISTANCE_THRESHOLD_UM / UM_PER_PIXEL
            proximal_coords = all_inf_coords[dists <= threshold_px]
        else:
            proximal_coords = np.empty((0, 2))

        samples.append({
            "sample_id": sample_id,
            "all_inf_coords": all_inf_coords,
            "proximal_coords": proximal_coords,
        })

    return samples


def run_dbscan_grid(coords: np.ndarray, eps: float, min_samples: int) -> float:
    """Run DBSCAN and return cluster_fraction."""
    if len(coords) < MIN_CLUSTER_SIZE:
        return 0.0

    labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(coords)

    unique_labels = set(labels)
    unique_labels.discard(-1)

    clustered = sum(
        np.sum(labels == lab)
        for lab in unique_labels
        if np.sum(labels == lab) >= MIN_CLUSTER_SIZE
    )

    return float(clustered / len(coords))


def run_cox_pfs(feature_series: pd.Series, survival_df: pd.DataFrame) -> dict:
    """Run univariate Cox for PFS on a feature series aligned by participant_id."""
    df = survival_df[["participant_id", "pfs_time", "pfs_event"]].copy()
    df = df.merge(feature_series.rename("feature").reset_index(), on="participant_id")
    df = df[["feature", "pfs_time", "pfs_event"]].dropna()

    if len(df) < 10:
        return {"hr": np.nan, "p": np.nan, "c_index": np.nan, "n": len(df)}

    std = df["feature"].std()
    if std == 0:
        return {"hr": np.nan, "p": np.nan, "c_index": np.nan, "n": len(df)}

    df["feature"] = (df["feature"] - df["feature"].mean()) / std

    try:
        cph = CoxPHFitter()
        cph.fit(df, duration_col="pfs_time", event_col="pfs_event")
        s = cph.summary
        return {
            "hr": s.loc["feature", "exp(coef)"],
            "ci_lower": s.loc["feature", "exp(coef) lower 95%"],
            "ci_upper": s.loc["feature", "exp(coef) upper 95%"],
            "p": s.loc["feature", "p"],
            "c_index": cph.concordance_index_,
            "n": len(df),
        }
    except Exception:
        return {"hr": np.nan, "p": np.nan, "c_index": np.nan, "n": len(df)}


def main():
    project_root = get_project_root()
    figures_dir = project_root / "figures"
    data_dir = project_root / "data"
    figures_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("DBSCAN Parameter Sensitivity Analysis")
    print("=" * 70)

    # Load survival dataset
    survival_path = data_dir / "survival_analysis_dataset.csv"
    if not survival_path.exists():
        print("ERROR: Run survival_analysis.py first to build survival dataset")
        return
    survival_df = pd.read_csv(survival_path)

    # Get sample → participant mapping
    consep = load_cell_densities("consep")
    metadata = get_sample_metadata(consep)

    # Load all sample coordinates
    print("\n1. Loading sample data and computing tumor-proximal coords...")
    samples = load_all_samples()

    # Run parameter grid
    print(f"\n2. Running DBSCAN grid: {len(EPS_VALUES)} eps × {len(MIN_SAMPLES_VALUES)} min_samples = {len(EPS_VALUES) * len(MIN_SAMPLES_VALUES)} combinations")

    results = []
    for eps in EPS_VALUES:
        for min_samples in MIN_SAMPLES_VALUES:
            # Compute cluster_fraction for each sample at this parameter combo
            rows = []
            for s in samples:
                all_frac = run_dbscan_grid(s["all_inf_coords"], eps, min_samples)
                tumor_frac = run_dbscan_grid(s["proximal_coords"], eps, min_samples)
                rows.append({
                    "sample_id": s["sample_id"],
                    "cluster_fraction": all_frac,
                    "tumor_cluster_fraction": tumor_frac,
                })

            sample_df = pd.DataFrame(rows)
            sample_df = sample_df.merge(metadata[["sample_id", "participant_id"]], on="sample_id")

            # Aggregate to participant level
            participant_df = sample_df.groupby("participant_id")[
                ["cluster_fraction", "tumor_cluster_fraction"]
            ].mean()

            # Cox PFS for both variants
            all_cox = run_cox_pfs(participant_df["cluster_fraction"], survival_df)
            tumor_cox = run_cox_pfs(participant_df["tumor_cluster_fraction"], survival_df)

            results.append({
                "eps": eps,
                "min_samples": min_samples,
                "all_hr": all_cox.get("hr"),
                "all_p": all_cox.get("p"),
                "all_c": all_cox.get("c_index"),
                "tumor_hr": tumor_cox.get("hr"),
                "tumor_p": tumor_cox.get("p"),
                "tumor_c": tumor_cox.get("c_index"),
            })

            sig_all = "*" if all_cox.get("p", 1) < 0.05 else " "
            sig_tumor = "*" if tumor_cox.get("p", 1) < 0.05 else " "
            print(
                f"  eps={eps:3d} min_samples={min_samples:2d}  |  "
                f"{sig_all} all: HR={all_cox.get('hr', np.nan):.3f} p={all_cox.get('p', np.nan):.4f} C={all_cox.get('c_index', np.nan):.3f}  |  "
                f"{sig_tumor} tumor: HR={tumor_cox.get('hr', np.nan):.3f} p={tumor_cox.get('p', np.nan):.4f} C={tumor_cox.get('c_index', np.nan):.3f}"
            )

    results_df = pd.DataFrame(results)
    results_path = data_dir / "dbscan_sensitivity.csv"
    results_df.to_csv(results_path, index=False)
    print(f"\n  Saved: {results_path}")

    # Plot heatmaps
    print("\n3. Plotting heatmaps...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    for row, (prefix, label) in enumerate([("all", "All cells"), ("tumor", "Tumor-proximal")]):
        for col, (metric, metric_label, cmap, vrange) in enumerate([
            ("p", "p-value (PFS)", "RdYlGn_r", (0, 0.15)),
            ("hr", "HR per SD (PFS)", "RdBu", (0.2, 0.8)),
            ("c", "C-index (PFS)", "RdYlGn", (0.55, 0.80)),
        ]):
            ax = axes[row, col]
            col_name = f"{prefix}_{metric}"

            pivot = results_df.pivot(index="min_samples", columns="eps", values=col_name)
            pivot = pivot.sort_index(ascending=False)  # Larger min_samples at top

            im = ax.imshow(
                pivot.values,
                cmap=cmap,
                vmin=vrange[0],
                vmax=vrange[1],
                aspect="auto",
            )

            # Annotate cells
            for i in range(pivot.shape[0]):
                for j in range(pivot.shape[1]):
                    val = pivot.values[i, j]
                    if metric == "p":
                        text = f"{val:.4f}"
                        bold = val < 0.05
                    elif metric == "hr":
                        text = f"{val:.3f}"
                        bold = False
                    else:
                        text = f"{val:.3f}"
                        bold = False

                    fontweight = "bold" if bold else "normal"
                    ax.text(j, i, text, ha="center", va="center",
                            fontsize=9, fontweight=fontweight)

            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels(pivot.columns)
            ax.set_yticks(range(len(pivot.index)))
            ax.set_yticklabels(pivot.index)
            ax.set_xlabel("eps (pixels)")
            ax.set_ylabel("min_samples")
            ax.set_title(f"{label} — {metric_label}")
            plt.colorbar(im, ax=ax, shrink=0.8)

    plt.suptitle("DBSCAN Parameter Sensitivity — cluster_fraction vs PFS", fontsize=14)
    plt.tight_layout()
    output_path = figures_dir / "dbscan_sensitivity.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    n_sig_all = (results_df["all_p"] < 0.05).sum()
    n_sig_tumor = (results_df["tumor_p"] < 0.05).sum()
    n_total = len(results_df)
    print(f"All-cell DBSCAN:      {n_sig_all}/{n_total} parameter combos significant (p<0.05)")
    print(f"Tumor-proximal DBSCAN: {n_sig_tumor}/{n_total} parameter combos significant (p<0.05)")
    print(f"All-cell p range:      {results_df['all_p'].min():.4f} — {results_df['all_p'].max():.4f}")
    print(f"Tumor-proximal p range: {results_df['tumor_p'].min():.4f} — {results_df['tumor_p'].max():.4f}")


if __name__ == "__main__":
    main()
