#!/usr/bin/env python3
"""
Match ClassPose cell centroids to Orion mIF (CellRingMask) data and generate
validation figures.

The Orion mIF data for one sample (P37/S29) lives in SAMPLE_DATA_DIR as
(override with the MIF_CSV env var, see config.py):
  P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv

ClassPose (ConSeP) centroids for the matching sample live in:
  sample_data/dg_4DFC_CD2C5DCB-9095-4E90-94F4-1FD37578830D/
    dg_4DFC_CD2C5DCB-9095-4E90-94F4-1FD37578830D_consep_spatialdata.zarr/
      points/cell_centroids/points.parquet

Both coordinate systems are in pixels at the same resolution, so matching is
done by nearest-neighbour centroid distance with a 10 px threshold (~3.25 µm).

Outputs (saved to figures/):
  match_distance_diagnostic.png  — histogram + CDF of NN distances
  match_area_correlation.png     — ClassPose area vs CellRingMask area scatter
  marker_boxplots_by_class.png   — boxplots of 16 mIF markers per cell class
  marker_violins_by_class.png    — violin plots of 16 mIF markers per cell class
  cd45_scatter.png               — spatial scatter coloured by CD45 intensity

Usage:
  uv run python scripts/match_mif.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from config import MIF_CSV, SAMPLE_DATA_DIR

# ── Paths ────────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# MIF_CSV comes from config.py (default: SAMPLE_DATA_DIR / <cellRingMask csv>)

SAMPLE_ID = "dg_4DFC_CD2C5DCB-9095-4E90-94F4-1FD37578830D"
CLASSPOSE_PARQUET = (
    SAMPLE_DATA_DIR
    / SAMPLE_ID
    / f"{SAMPLE_ID}_consep_spatialdata.zarr"
    / "points"
    / "cell_centroids"
    / "points.parquet"
)

FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────

MAX_DIST_PX = 10.0  # ~3.25 µm at 0.325 µm/px
MAX_VIOLIN_PER_CLASS = 50_000  # subsample ceiling for violin plotting

MARKERS = [
    "CD45", "CD3e", "CD4", "CD8a", "CD20", "CD68", "CD163", "PD-L1",
    "PD-1", "FOXP3", "CD45RO", "Pan-CK", "E-cadherin", "SMA", "Ki67", "CD31",
]

CLASS_ORDER = [
    "Inflammatory", "Stroma", "Healthy epithelial",
    "Malignant epithelial", "Muscle", "Other",
]

CLASS_COLORS = {
    "Inflammatory": "#e74c3c",
    "Stroma": "#3498db",
    "Healthy epithelial": "#2ecc71",
    "Malignant epithelial": "#9b59b6",
    "Muscle": "#e67e22",
    "Other": "#95a5a6",
}


# ── Loading ───────────────────────────────────────────────────────────────────

def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    print("Loading CellRingMask CSV...")
    mif = pd.read_csv(MIF_CSV)
    print(f"  {len(mif):,} cells, {len(mif.columns)} columns")

    print("Loading ClassPose parquet...")
    pts = pd.read_parquet(CLASSPOSE_PARQUET)
    for col in ("x", "y", "area"):
        pts[col] = pts[col].astype(float)
    pts = pts.reset_index(drop=True)
    print(f"  {len(pts):,} cells")
    print(f"  Classes: {pts['classification_name'].value_counts().to_dict()}")

    return mif, pts


# ── Matching ──────────────────────────────────────────────────────────────────

def match_centroids(mif: pd.DataFrame, pts: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return (distances, indices) arrays from KDTree NN query."""
    print("\nBuilding KD-tree on mIF centroids...")
    tree = cKDTree(mif[["X_centroid", "Y_centroid"]].values)

    print("Querying nearest neighbours for ClassPose cells...")
    distances, indices = tree.query(pts[["x", "y"]].values, k=1)

    for thresh in (2, 5, 10, 20):
        frac = (distances < thresh).mean() * 100
        print(f"  Within {thresh} px: {frac:.1f}%")

    return distances, indices


# ── Figures ───────────────────────────────────────────────────────────────────

def plot_distance_diagnostic(distances: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].hist(distances, bins=200, range=(0, 50), edgecolor="none", color="steelblue")
    axes[0].axvline(MAX_DIST_PX, color="red", ls="--", label=f"{MAX_DIST_PX:.0f} px threshold")
    axes[0].set_xlabel("Nearest-neighbor distance (px)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Match distance distribution")
    axes[0].legend()

    sorted_d = np.sort(distances)
    cdf = np.arange(1, len(sorted_d) + 1) / len(sorted_d)
    axes[1].plot(sorted_d, cdf, color="steelblue", lw=1.5)
    axes[1].axvline(MAX_DIST_PX, color="red", ls="--", label=f"{MAX_DIST_PX:.0f} px threshold")
    axes[1].axhline((distances < MAX_DIST_PX).mean(), color="red", ls=":", alpha=0.5)
    axes[1].set_xlim(0, 50)
    axes[1].set_xlabel("Distance (px)")
    axes[1].set_ylabel("Cumulative fraction")
    axes[1].set_title("CDF of match distances")
    axes[1].legend()

    plt.tight_layout()
    out = FIGURES_DIR / "match_distance_diagnostic.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")


def plot_area_correlation(
    pts: pd.DataFrame, mif: pd.DataFrame, distances: np.ndarray, indices: np.ndarray
) -> None:
    good = distances < MAX_DIST_PX
    cp_area = pts["area"].values
    mif_area = mif["Area"].values[indices]

    n = good.sum()
    r = np.corrcoef(cp_area[good], mif_area[good])[0, 1]
    print(f"\nArea correlation (n={n:,}, r={r:.4f})")

    rng = np.random.default_rng(42)
    idx = rng.choice(np.where(good)[0], size=min(50_000, n), replace=False)

    fig, ax = plt.subplots(figsize=(7, 10))
    ax.scatter(cp_area[idx], mif_area[idx], s=0.5, alpha=0.3, c="steelblue", edgecolors="none")
    lim = max(np.percentile(cp_area[good], 99), np.percentile(mif_area[good], 99))
    ax.plot([0, lim], [0, lim], "r--", lw=1)
    ax.set_xlim(0, np.percentile(cp_area[good], 99))
    ax.set_ylim(0, np.percentile(mif_area[good], 99))
    ax.set_xlabel("ClassPose area")
    ax.set_ylabel("CellRingMask area")
    ax.set_title(f"Area correlation (r={r:.3f}, n={n:,})")

    plt.tight_layout()
    out = FIGURES_DIR / "match_area_correlation.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")


def _build_marker_df(
    pts: pd.DataFrame, mif: pd.DataFrame, distances: np.ndarray, indices: np.ndarray
) -> pd.DataFrame:
    good = distances < MAX_DIST_PX
    matched_pts = pts[good].reset_index(drop=True)
    matched_mif_idx = indices[good]
    marker_data = mif.loc[matched_mif_idx, MARKERS].reset_index(drop=True)
    marker_data["cell_class"] = matched_pts["classification_name"].values
    return marker_data


def plot_boxplots(marker_data: pd.DataFrame, n_matched: int) -> None:
    present = [c for c in CLASS_ORDER if c in marker_data["cell_class"].unique()]

    fig, axes = plt.subplots(4, 4, figsize=(20, 16))
    for i, marker in enumerate(MARKERS):
        ax = axes.flat[i]
        data = [np.log1p(marker_data.loc[marker_data["cell_class"] == c, marker].values) for c in present]
        bp = ax.boxplot(data, labels=present, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=1.5))
        for patch, c in zip(bp["boxes"], present):
            patch.set_facecolor(CLASS_COLORS[c])
            patch.set_alpha(0.7)
        ax.set_title(marker, fontsize=11, fontweight="bold")
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.set_ylabel("log(intensity + 1)", fontsize=8)

    plt.suptitle(
        f"Marker intensities by ClassPose cell class (n={n_matched:,} matched cells)",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIGURES_DIR / "marker_boxplots_by_class.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")


def plot_violins(marker_data: pd.DataFrame, n_matched: int) -> None:
    present = [c for c in CLASS_ORDER if c in marker_data["cell_class"].unique()]

    # Subsample per class
    rng = np.random.default_rng(42)
    frames = []
    for c in present:
        sub = marker_data[marker_data["cell_class"] == c]
        if len(sub) > MAX_VIOLIN_PER_CLASS:
            sub = sub.sample(n=MAX_VIOLIN_PER_CLASS, random_state=42)
        frames.append(sub)
    sampled = pd.concat(frames, ignore_index=True)
    print(f"\nViolin subsample: {len(sampled):,} cells")

    fig, axes = plt.subplots(4, 4, figsize=(22, 18))
    for i, marker in enumerate(MARKERS):
        ax = axes.flat[i]
        data = [np.log1p(sampled.loc[sampled["cell_class"] == c, marker].values) for c in present]
        vp = ax.violinplot(data, positions=range(len(present)), showmedians=True, showextrema=False)
        for j, body in enumerate(vp["bodies"]):
            body.set_facecolor(CLASS_COLORS[present[j]])
            body.set_alpha(0.7)
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.5)
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels(present, rotation=45, ha="right", fontsize=8)
        ax.set_title(marker, fontsize=11, fontweight="bold")
        ax.set_ylabel("log(intensity + 1)", fontsize=8)

    plt.suptitle(
        f"Marker intensities by ClassPose cell class (n={n_matched:,} matched cells)",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIGURES_DIR / "marker_violins_by_class.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")


def plot_cd45_scatter(mif: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))
    sc = ax.scatter(
        mif["X_centroid"], mif["Y_centroid"],
        c=mif["CD45"], cmap="magma", s=0.5, edgecolors="none",
        norm=mcolors.LogNorm(),
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("X_centroid")
    ax.set_ylabel("Y_centroid")
    ax.set_title("Cells colored by CD45 intensity (log scale)")
    plt.colorbar(sc, ax=ax, label="CD45")
    plt.tight_layout()
    out = FIGURES_DIR / "cd45_scatter.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    mif, pts = load_data()
    distances, indices = match_centroids(mif, pts)

    n_matched = (distances < MAX_DIST_PX).sum()
    print(f"\nMatched {n_matched:,} / {len(pts):,} ClassPose cells ({n_matched/len(pts)*100:.1f}%)")

    print("\nGenerating figures...")
    plot_distance_diagnostic(distances)
    plot_area_correlation(pts, mif, distances, indices)

    marker_data = _build_marker_df(pts, mif, distances, indices)
    plot_boxplots(marker_data, n_matched)
    plot_violins(marker_data, n_matched)
    plot_cd45_scatter(mif)

    print("\nDone.")


if __name__ == "__main__":
    main()
