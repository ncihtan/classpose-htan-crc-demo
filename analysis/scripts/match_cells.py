#!/usr/bin/env python3
"""
Match cells between CoNIC and ConSep segmentation models.

Both models were run independently on the same H&E whole-slide images,
producing ~1.2M cell detections each. This script matches cells detected
by one model to cells detected by the other using:
  1. Centroid distance matching (fast, default)
  2. Contour IoU matching (slower, optional)

Usage:
    python scripts/match_cells.py --sample SAMPLE_ID [--contours] [--max-distance 10.0] [--min-iou 0.3]
    python scripts/match_cells.py --all [--contours]
"""

import argparse
import subprocess
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import STRtree
from tqdm import tqdm

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))
from config import AWS_PROFILE, S3_BASE, SAMPLE_DATA_DIR
from load_cell_densities import get_project_root
from utils.spatial_utils import UM_PER_PIXEL

# Constants
MAX_CENTROID_DISTANCE_PX = 10.0  # ~3.25 um, ~1 cell radius
MIN_IOU = 0.3

# S3_BASE and AWS_PROFILE come from config.py (env-var overridable)


# ── S3 Data Sync ───────────────────────────────────────────────────────────


def _s3_zarr_prefix(sample_id: str, model: str) -> str:
    """Return the S3 prefix for a sample's spatialdata zarr."""
    return (
        f"{S3_BASE}/{sample_id}/{model}/"
        f"{sample_id}_{model}_spatialdata.zarr"
    )


def _local_zarr_dir(sample_id: str, model: str) -> Path:
    """Return the local path for a sample's spatialdata zarr."""
    return (
        SAMPLE_DATA_DIR
        / sample_id
        / f"{sample_id}_{model}_spatialdata.zarr"
    )


def sync_parquet_from_s3(
    sample_id: str, model: str, include_contours: bool = False
) -> None:
    """
    Download centroid (and optionally contour) parquets from S3 if not local.

    Only fetches the specific parquet files needed, not the full zarr.
    Centroids are Dask-partitioned (points.parquet/part.0.parquet).
    Contours are a single file (shapes.parquet).

    Parameters
    ----------
    sample_id : str
        Sample identifier
    model : str
        Model name ('conic' or 'consep')
    include_contours : bool
        Whether to also fetch the contour shapes parquet
    """
    s3_base = _s3_zarr_prefix(sample_id, model)
    local_base = _local_zarr_dir(sample_id, model)

    # (s3_rel_path, local_rel_path, label, is_recursive)
    files_to_sync: list[tuple[str, str, str, bool]] = [
        # Centroids: Dask-partitioned directory — sync the whole directory
        (
            "points/cell_centroids/points.parquet/",
            "points/cell_centroids/points.parquet",
            "centroids",
            True,
        ),
    ]
    if include_contours:
        files_to_sync.append(
            (
                "shapes/cell_contours/shapes.parquet",
                "shapes/cell_contours/shapes.parquet",
                "contours",
                False,
            ),
        )

    for s3_rel, local_rel, label, recursive in files_to_sync:
        local_path = local_base / local_rel
        if local_path.exists():
            continue

        s3_uri = f"{s3_base}/{s3_rel}"
        local_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"    Downloading {model} {label} from S3...")

        cmd = ["aws", "s3", "cp"]
        if AWS_PROFILE:
            cmd.extend(["--profile", AWS_PROFILE])
        if recursive:
            cmd.append("--recursive")
        cmd.extend([s3_uri, str(local_path) + ("/" if recursive else "")])

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise FileNotFoundError(
                f"Failed to download {s3_uri}: {result.stderr.strip()}"
            )


# ── Data Loading ────────────────────────────────────────────────────────────


def get_sample_ids() -> list[str]:
    """Read sample IDs from data/sample_manifest.txt."""
    manifest = get_project_root() / "data" / "sample_manifest.txt"
    return manifest.read_text().strip().splitlines()


def load_centroids(sample_id: str, model: str) -> pd.DataFrame:
    """
    Load cell centroids from a SpatialData zarr via direct parquet read.

    Parameters
    ----------
    sample_id : str
        Sample identifier (e.g., 'dg_4DFC_...')
    model : str
        Model name ('conic' or 'consep')

    Returns
    -------
    pd.DataFrame
        DataFrame with x, y, classification_name, area columns.
        Index is reset for row-alignment with contours.
    """
    parquet_path = (
        SAMPLE_DATA_DIR
        / sample_id
        / f"{sample_id}_{model}_spatialdata.zarr"
        / "points"
        / "cell_centroids"
        / "points.parquet"
    )
    df = pd.read_parquet(
        parquet_path, columns=["x", "y", "classification_name", "area"]
    )
    return df.reset_index(drop=True)


def load_contours(sample_id: str, model: str) -> gpd.GeoDataFrame:
    """
    Load cell contours from a SpatialData zarr via direct parquet read.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    model : str
        Model name ('conic' or 'consep')

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with geometry (WKB polygons) and classification_name.
        Index is reset for row-alignment with centroids.
    """
    parquet_path = (
        SAMPLE_DATA_DIR
        / sample_id
        / f"{sample_id}_{model}_spatialdata.zarr"
        / "shapes"
        / "cell_contours"
        / "shapes.parquet"
    )
    gdf = gpd.read_parquet(parquet_path)
    return gdf.reset_index(drop=True)


# ── Centroid Distance Matching ──────────────────────────────────────────────


def match_by_centroid_distance(
    xy_a: np.ndarray, xy_b: np.ndarray, max_distance: float = MAX_CENTROID_DISTANCE_PX
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Match cells by nearest-centroid distance.

    Parameters
    ----------
    xy_a : np.ndarray
        (N, 2) centroid coordinates for model A
    xy_b : np.ndarray
        (M, 2) centroid coordinates for model B
    max_distance : float
        Maximum distance in pixels to consider a match

    Returns
    -------
    tuple of (a_indices, b_indices, distances)
        Matched pairs after filtering by distance and resolving many-to-one.
    """
    tree_b = cKDTree(xy_b)
    distances, b_indices = tree_b.query(xy_a, k=1)

    # Filter pairs exceeding max distance
    mask = distances <= max_distance
    a_idx = np.where(mask)[0]
    b_idx = b_indices[mask]
    dist = distances[mask]

    # Resolve many-to-one: when multiple A-cells match the same B-cell,
    # keep only the closest pair
    pairs = pd.DataFrame({"a_idx": a_idx, "b_idx": b_idx, "distance": dist})
    best = pairs.loc[pairs.groupby("b_idx")["distance"].idxmin()]

    return (
        best["a_idx"].values.astype(np.int64),
        best["b_idx"].values.astype(np.int64),
        best["distance"].values,
    )


# ── Contour Overlap Matching ───────────────────────────────────────────────


def match_by_contour_iou(
    contours_a: gpd.GeoDataFrame,
    contours_b: gpd.GeoDataFrame,
    min_iou: float = MIN_IOU,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Match cells by contour Intersection-over-Union.

    Parameters
    ----------
    contours_a : gpd.GeoDataFrame
        Contours for model A (reset index)
    contours_b : gpd.GeoDataFrame
        Contours for model B (reset index)
    min_iou : float
        Minimum IoU to consider a match

    Returns
    -------
    tuple of (a_indices, b_indices, ious)
        Matched pairs after filtering by IoU and resolving many-to-one.
    """
    geoms_b = contours_b.geometry.values
    tree_b = STRtree(geoms_b)

    a_indices = []
    b_indices = []
    ious = []

    for a_i, geom_a in enumerate(tqdm(contours_a.geometry, desc="IoU matching", leave=False)):
        if geom_a is None or geom_a.is_empty:
            continue

        # Query STRtree for bounding-box overlaps
        candidates = tree_b.query(geom_a)
        if len(candidates) == 0:
            continue

        best_iou = 0.0
        best_b_i = -1

        for b_i in candidates:
            geom_b = geoms_b[b_i]
            if geom_b is None or geom_b.is_empty:
                continue
            try:
                intersection = geom_a.intersection(geom_b).area
                union = geom_a.area + geom_b.area - intersection
                if union > 0:
                    iou = intersection / union
                    if iou > best_iou:
                        best_iou = iou
                        best_b_i = b_i
            except Exception:
                continue

        if best_iou >= min_iou and best_b_i >= 0:
            a_indices.append(a_i)
            b_indices.append(best_b_i)
            ious.append(best_iou)

    if not a_indices:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([])

    # Resolve many-to-one: keep highest IoU per B-cell
    pairs = pd.DataFrame(
        {"a_idx": a_indices, "b_idx": b_indices, "iou": ious}
    )
    best = pairs.loc[pairs.groupby("b_idx")["iou"].idxmax()]

    return (
        best["a_idx"].values.astype(np.int64),
        best["b_idx"].values.astype(np.int64),
        best["iou"].values,
    )


# ── Combined Matching ──────────────────────────────────────────────────────


def process_sample(
    sample_id: str,
    include_contours: bool = False,
    max_distance: float = MAX_CENTROID_DISTANCE_PX,
    min_iou: float = MIN_IOU,
) -> dict:
    """
    Match cells between CoNIC and ConSep for a single sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    include_contours : bool
        Whether to also run IoU-based contour matching
    max_distance : float
        Maximum centroid distance in pixels
    min_iou : float
        Minimum contour IoU threshold

    Returns
    -------
    dict
        Summary statistics for this sample
    """
    project_root = get_project_root()
    output_dir = project_root / "data" / "cell_matching"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sync parquets from S3 if not present locally
    print(f"  Syncing data from S3...")
    sync_parquet_from_s3(sample_id, "conic", include_contours=include_contours)
    sync_parquet_from_s3(sample_id, "consep", include_contours=include_contours)

    # Load centroids
    print(f"  Loading centroids...")
    conic_centroids = load_centroids(sample_id, "conic")
    consep_centroids = load_centroids(sample_id, "consep")
    print(f"    CoNIC: {len(conic_centroids):,} cells")
    print(f"    ConSep: {len(consep_centroids):,} cells")

    # Centroid distance matching
    print(f"  Running centroid distance matching (max {max_distance} px)...")
    xy_conic = conic_centroids[["x", "y"]].values
    xy_consep = consep_centroids[["x", "y"]].values
    cent_a, cent_b, cent_dist = match_by_centroid_distance(
        xy_conic, xy_consep, max_distance
    )
    print(f"    Matched: {len(cent_a):,} pairs")

    centroid_df = pd.DataFrame({
        "conic_idx": cent_a,
        "consep_idx": cent_b,
        "centroid_distance": cent_dist,
    })

    # Optional contour IoU matching
    iou_df = None
    if include_contours:
        print(f"  Loading contours...")
        conic_contours = load_contours(sample_id, "conic")
        consep_contours = load_contours(sample_id, "consep")
        print(f"  Running contour IoU matching (min IoU {min_iou})...")
        iou_a, iou_b, iou_vals = match_by_contour_iou(
            conic_contours, consep_contours, min_iou
        )
        print(f"    Matched: {len(iou_a):,} pairs")
        iou_df = pd.DataFrame({
            "conic_idx": iou_a,
            "consep_idx": iou_b,
            "iou": iou_vals,
        })

    # Merge match sets
    if iou_df is not None and len(iou_df) > 0:
        matched = centroid_df.merge(
            iou_df, on=["conic_idx", "consep_idx"], how="outer"
        )
        # Tag match method
        has_centroid = matched["centroid_distance"].notna()
        has_iou = matched["iou"].notna()
        matched["match_method"] = np.where(
            has_centroid & has_iou,
            "both",
            np.where(has_centroid, "centroid_only", "iou_only"),
        )
    else:
        matched = centroid_df.copy()
        matched["iou"] = np.nan
        matched["match_method"] = "centroid_only"

    # Annotate with classification names, areas, coordinates
    matched["conic_class"] = conic_centroids.loc[
        matched["conic_idx"].values, "classification_name"
    ].values
    matched["consep_class"] = consep_centroids.loc[
        matched["consep_idx"].values, "classification_name"
    ].values
    matched["conic_x"] = conic_centroids.loc[matched["conic_idx"].values, "x"].values
    matched["conic_y"] = conic_centroids.loc[matched["conic_idx"].values, "y"].values
    matched["conic_area"] = conic_centroids.loc[
        matched["conic_idx"].values, "area"
    ].values
    matched["consep_area"] = consep_centroids.loc[
        matched["consep_idx"].values, "area"
    ].values

    # Reorder columns to match output schema
    col_order = [
        "conic_idx",
        "consep_idx",
        "conic_class",
        "consep_class",
        "centroid_distance",
        "iou",
        "match_method",
        "conic_x",
        "conic_y",
        "conic_area",
        "consep_area",
    ]
    matched = matched[col_order]

    # Save
    output_path = output_dir / f"{sample_id}_matched_cells.parquet"
    matched.to_parquet(output_path, index=False)
    print(f"  Saved: {output_path}")

    # Compute and return summary
    summary = compute_summary(
        sample_id, len(conic_centroids), len(consep_centroids), matched
    )
    return summary


# ── Summary Statistics ─────────────────────────────────────────────────────


def compute_summary(
    sample_id: str, conic_n: int, consep_n: int, matched_df: pd.DataFrame
) -> dict:
    """
    Compute matching summary statistics.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    conic_n : int
        Total CoNIC cells
    consep_n : int
        Total ConSep cells
    matched_df : pd.DataFrame
        Matched cell pairs DataFrame

    Returns
    -------
    dict
        Summary statistics
    """
    n_matched = len(matched_df)
    method_counts = matched_df["match_method"].value_counts()

    summary = {
        "sample_id": sample_id,
        "conic_total": conic_n,
        "consep_total": consep_n,
        "n_matched": n_matched,
        "n_centroid_only": int(method_counts.get("centroid_only", 0)),
        "n_iou_only": int(method_counts.get("iou_only", 0)),
        "n_both": int(method_counts.get("both", 0)),
        "conic_match_rate": n_matched / conic_n if conic_n > 0 else 0.0,
        "consep_match_rate": n_matched / consep_n if consep_n > 0 else 0.0,
    }

    # Distance statistics (for centroid matches)
    dist_vals = matched_df["centroid_distance"].dropna()
    if len(dist_vals) > 0:
        summary["mean_distance_px"] = dist_vals.mean()
        summary["median_distance_px"] = dist_vals.median()
        summary["mean_distance_um"] = dist_vals.mean() * UM_PER_PIXEL
        summary["median_distance_um"] = dist_vals.median() * UM_PER_PIXEL

    # IoU statistics (for contour matches)
    iou_vals = matched_df["iou"].dropna()
    if len(iou_vals) > 0:
        summary["mean_iou"] = iou_vals.mean()
        summary["median_iou"] = iou_vals.median()

    return summary


def save_summary(summaries: list[dict]) -> None:
    """Append summaries to the matching summary CSV."""
    output_path = get_project_root() / "data" / "cell_matching" / "matching_summary.csv"
    df = pd.DataFrame(summaries)
    df.to_csv(output_path, index=False)
    print(f"\nSummary saved to: {output_path}")


# ── CLI ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Match cells between CoNIC and ConSep segmentations."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample", type=str, help="Process a single sample ID")
    group.add_argument("--all", action="store_true", help="Process all samples")
    parser.add_argument(
        "--contours", action="store_true", help="Include contour IoU matching"
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=MAX_CENTROID_DISTANCE_PX,
        help=f"Max centroid distance in pixels (default: {MAX_CENTROID_DISTANCE_PX})",
    )
    parser.add_argument(
        "--min-iou",
        type=float,
        default=MIN_IOU,
        help=f"Min contour IoU threshold (default: {MIN_IOU})",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Cell-Cell Matching: CoNIC vs ConSep")
    print("=" * 70)
    print(f"  Max centroid distance: {args.max_distance} px ({args.max_distance * UM_PER_PIXEL:.2f} um)")
    print(f"  Contour matching: {'enabled' if args.contours else 'disabled'}")
    if args.contours:
        print(f"  Min IoU: {args.min_iou}")
    print()

    if args.sample:
        sample_ids = [args.sample]
    else:
        sample_ids = get_sample_ids()

    summaries = []
    for sample_id in tqdm(sample_ids, desc="Samples", disable=len(sample_ids) == 1):
        print(f"Processing: {sample_id}")
        try:
            summary = process_sample(
                sample_id,
                include_contours=args.contours,
                max_distance=args.max_distance,
                min_iou=args.min_iou,
            )
            summaries.append(summary)

            # Print per-sample summary
            print(f"    Matched: {summary['n_matched']:,} pairs")
            print(f"    CoNIC match rate: {summary['conic_match_rate']:.1%}")
            print(f"    ConSep match rate: {summary['consep_match_rate']:.1%}")
            if "median_distance_px" in summary:
                print(f"    Median distance: {summary['median_distance_px']:.2f} px ({summary['median_distance_um']:.2f} um)")
            if "median_iou" in summary:
                print(f"    Median IoU: {summary['median_iou']:.3f}")
            print()
        except FileNotFoundError as e:
            print(f"  Skipped (data not found): {e}")
            print()
            continue
        except Exception as e:
            print(f"  Error: {e}")
            print()
            continue

    if summaries:
        save_summary(summaries)

    print("Done!")


if __name__ == "__main__":
    main()
