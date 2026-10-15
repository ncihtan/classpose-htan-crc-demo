#!/usr/bin/env python3
"""
Select representative samples for detailed visualization in Figures 2-3.

Identifies immune-hot and immune-cold samples based on inflammatory cell
percentage from ConSeP model outputs.
"""

import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from config import SAMPLE_DATA_DIR
from load_cell_densities import (
    load_cell_densities,
    get_sample_metadata,
    pivot_for_stacked_bar,
    CONSEP_IMMUNE_CLASSES,
    get_project_root,
)


def calculate_immune_percentage(consep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate inflammatory cell percentage for each sample.

    Returns DataFrame with sample_id, participant_id, inflammatory_pct, total_cells.
    """
    meta = get_sample_metadata(consep_df)

    # Pivot to get counts per cell type
    pivoted = pivot_for_stacked_bar(consep_df, by="sample_id")

    # Calculate totals and percentages
    results = []
    for sample_id in pivoted.index:
        total = pivoted.loc[sample_id].sum()
        inflammatory = pivoted.loc[sample_id].get("Inflammatory", 0)
        inflammatory_pct = (inflammatory / total * 100) if total > 0 else 0

        # Get participant info
        sample_meta = meta[meta["sample_id"] == sample_id].iloc[0]

        results.append({
            "sample_id": sample_id,
            "participant_id": sample_meta["participant_id"],
            "biospecimen_id": sample_meta["biospecimen_id"],
            "inflammatory_count": int(inflammatory),
            "total_cells": int(total),
            "inflammatory_pct": inflammatory_pct,
        })

    return pd.DataFrame(results).sort_values("inflammatory_pct", ascending=False)


def load_drs_mapping() -> dict:
    """Load mapping from sample_id to DRS metadata."""
    project_root = get_project_root()
    drs_df = pd.read_csv(project_root / "metadata" / "files_drs.csv")

    mapping = {}
    for _, row in drs_df.iterrows():
        drs_uri = row["drs_uri"]
        uuid = drs_uri.split("/")[-1]
        sample_id = f"dg_4DFC_{uuid}"
        mapping[sample_id] = {
            "drs_uri": drs_uri,
            "filename": row["name"],
            "data_file_id": row["data_file_id"],
            "biospecimen": row["biospecimen"],
        }
    return mapping


def select_samples():
    """Select representative immune-hot and immune-cold samples."""
    print("=" * 70)
    print("Sample Selection for Figures 2-3")
    print("=" * 70)

    # Load ConSeP data
    print("\nLoading ConSeP cell densities...")
    consep_df = load_cell_densities("consep")

    # Calculate immune percentages
    print("Calculating inflammatory percentages...")
    immune_df = calculate_immune_percentage(consep_df)

    # Load DRS mapping for file info
    drs_mapping = load_drs_mapping()

    # Add DRS info to dataframe
    immune_df["filename"] = immune_df["sample_id"].map(
        lambda x: drs_mapping.get(x, {}).get("filename", "")
    )
    immune_df["data_file_id"] = immune_df["sample_id"].map(
        lambda x: drs_mapping.get(x, {}).get("data_file_id", "")
    )

    # Print full ranking
    print("\n" + "=" * 70)
    print("All Samples Ranked by Inflammatory Cell Percentage (ConSeP)")
    print("=" * 70)
    print(f"{'Rank':<5} {'Participant':<12} {'Inflammatory %':<15} {'Cells':<12} {'File ID'}")
    print("-" * 70)

    for rank, (_, row) in enumerate(immune_df.iterrows(), 1):
        print(
            f"{rank:<5} {row['participant_id']:<12} "
            f"{row['inflammatory_pct']:>10.1f}%     "
            f"{row['total_cells']:>10,} "
            f" {row['data_file_id']}"
        )

    # Select top (immune-hot) and bottom (immune-cold)
    immune_hot = immune_df.iloc[0]
    immune_cold = immune_df.iloc[-1]

    print("\n" + "=" * 70)
    print("Selected Representative Samples")
    print("=" * 70)

    print("\nIMMUNE-HOT (highest inflammatory %):")
    print(f"  Sample ID: {immune_hot['sample_id']}")
    print(f"  Participant: {immune_hot['participant_id']}")
    print(f"  Inflammatory: {immune_hot['inflammatory_pct']:.1f}% ({immune_hot['inflammatory_count']:,} cells)")
    print(f"  Total cells: {immune_hot['total_cells']:,}")
    print(f"  File: {immune_hot['filename']}")

    print("\nIMMUNE-COLD (lowest inflammatory %):")
    print(f"  Sample ID: {immune_cold['sample_id']}")
    print(f"  Participant: {immune_cold['participant_id']}")
    print(f"  Inflammatory: {immune_cold['inflammatory_pct']:.1f}% ({immune_cold['inflammatory_count']:,} cells)")
    print(f"  Total cells: {immune_cold['total_cells']:,}")
    print(f"  File: {immune_cold['filename']}")

    # Check for available data in SAMPLE_DATA_DIR (see config.py)

    print("\n" + "=" * 70)
    print("Data Availability Check")
    print("=" * 70)

    # Check for HTA7_954 (sample used in prototypes)
    prototype_sample = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
    prototype_row = immune_df[immune_df["sample_id"] == prototype_sample]

    if not prototype_row.empty:
        proto = prototype_row.iloc[0]
        print(f"\nPrototype sample (HTA7_954) status:")
        print(f"  Inflammatory: {proto['inflammatory_pct']:.1f}%")
        print(f"  Rank: {list(immune_df['sample_id']).index(prototype_sample) + 1} of {len(immune_df)}")

        # Check if data files exist
        he_path = SAMPLE_DATA_DIR / prototype_sample / f"{prototype_sample}.tif"
        zarr_path = SAMPLE_DATA_DIR / prototype_sample / f"{prototype_sample}_consep_spatialdata.zarr"

        print(f"  H&E tif exists: {he_path.exists()}")
        print(f"  ConSeP zarr exists: {zarr_path.exists()}")

        if he_path.exists() and zarr_path.exists():
            print("\n  >> RECOMMENDATION: Use prototype sample (HTA7_954) for Figures 2-3")
            print("     Full SpatialData and H&E image already available locally.")

    # Save summary to data directory
    output_dir = get_project_root() / "data"
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / "sample_immune_ranking.csv"
    immune_df.to_csv(output_path, index=False)
    print(f"\nSaved ranking to: {output_path}")

    # Save selected samples
    selected = pd.DataFrame([
        {"role": "immune_hot", **immune_hot.to_dict()},
        {"role": "immune_cold", **immune_cold.to_dict()},
        {"role": "prototype", **proto.to_dict()} if not prototype_row.empty else {},
    ])
    selected_path = output_dir / "selected_samples.csv"
    selected.to_csv(selected_path, index=False)
    print(f"Saved selected samples to: {selected_path}")

    return immune_df, immune_hot, immune_cold


if __name__ == "__main__":
    immune_df, hot, cold = select_samples()
