#!/usr/bin/env python3
"""
Shared utility module for loading cell density data from classpose outputs.

Data is indexed by sample_id (unique file identifier) rather than participant_id,
since some participants have multiple H&E images.
"""

import pandas as pd
from pathlib import Path
import re


def get_project_root() -> Path:
    """Get the project root directory."""
    return Path(__file__).parent.parent


def load_cell_densities(model: str = "conic") -> pd.DataFrame:
    """
    Load all cell density CSVs for a given model (conic or consep).

    Parameters
    ----------
    model : str
        Either 'conic' or 'consep'

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: sample_id, participant_id, biospecimen_id,
        cell_class, count, density, model
    """
    data_dir = get_project_root() / "classpose_outputs" / "data"

    all_data = []

    for csv_file in data_dir.rglob(f"*/{model}/*_cell_densities.csv"):
        # Extract sample_id from filename (dg_4DFC_UUID format)
        match = re.search(r"(dg_4DFC_[A-F0-9-]+)", csv_file.name)
        sample_id = match.group(1) if match else csv_file.stem

        # Extract IDs from path: data/{participant}/{biospecimen}/{model}/
        parts = csv_file.relative_to(data_dir).parts
        participant_id = parts[0]
        biospecimen_id = parts[1]

        df = pd.read_csv(csv_file)
        df["sample_id"] = sample_id
        df["participant_id"] = participant_id
        df["biospecimen_id"] = biospecimen_id
        df["model"] = model

        all_data.append(df)

    if not all_data:
        raise FileNotFoundError(
            f"No cell density files found for model '{model}' in {data_dir}"
        )

    combined = pd.concat(all_data, ignore_index=True)

    # Sort by sample_id for consistent ordering
    combined = combined.sort_values(["sample_id", "cell_class"]).reset_index(
        drop=True
    )

    return combined


def load_both_models() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load cell densities for both conic and consep models.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (conic_df, consep_df)
    """
    conic_df = load_cell_densities("conic")
    consep_df = load_cell_densities("consep")
    return conic_df, consep_df


def pivot_for_stacked_bar(
    df: pd.DataFrame,
    sort_by_immune: bool = False,
    immune_classes: list[str] | None = None,
    by: str = "sample_id",
) -> pd.DataFrame:
    """
    Pivot data for stacked bar chart: samples as rows, cell classes as columns.

    Parameters
    ----------
    df : pd.DataFrame
        Output from load_cell_densities()
    sort_by_immune : bool
        If True, sort by percentage of immune cells (descending)
    immune_classes : list[str] | None
        List of cell classes considered "immune". Required if sort_by_immune=True.
    by : str
        Column to use as index ('sample_id' or 'participant_id'). Default: 'sample_id'

    Returns
    -------
    pd.DataFrame
        Pivoted DataFrame with sample_id as index, cell_class as columns, count as values
    """
    pivoted = df.pivot_table(
        index=by, columns="cell_class", values="count", aggfunc="sum"
    ).fillna(0)

    if sort_by_immune and immune_classes:
        # Calculate immune percentage
        total = pivoted.sum(axis=1)
        immune_cols = [c for c in immune_classes if c in pivoted.columns]
        immune_count = pivoted[immune_cols].sum(axis=1)
        pivoted["_immune_pct"] = (immune_count / total) * 100
        pivoted = pivoted.sort_values("_immune_pct", ascending=False)
        pivoted = pivoted.drop(columns=["_immune_pct"])
    else:
        # Sort by total cell count
        pivoted["_total"] = pivoted.sum(axis=1)
        pivoted = pivoted.sort_values("_total", ascending=False)
        pivoted = pivoted.drop(columns=["_total"])

    return pivoted


# Define immune cell classes for each model
CONIC_IMMUNE_CLASSES = ["Lymphocyte", "Plasma cell", "Neutrophil", "Eosinophil"]
CONSEP_IMMUNE_CLASSES = ["Inflammatory"]


def get_total_cells_per_sample(df: pd.DataFrame, by: str = "sample_id") -> pd.Series:
    """
    Get total cell count per sample.

    Parameters
    ----------
    df : pd.DataFrame
        Output from load_cell_densities()
    by : str
        Column to group by ('sample_id' or 'participant_id'). Default: 'sample_id'

    Returns
    -------
    pd.Series
        Series with sample_id as index, total count as values
    """
    return df.groupby(by)["count"].sum().sort_values(ascending=False)


def get_sample_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get metadata for each sample (mapping sample_id to participant/biospecimen).

    Parameters
    ----------
    df : pd.DataFrame
        Output from load_cell_densities()

    Returns
    -------
    pd.DataFrame
        DataFrame with sample_id, participant_id, biospecimen_id (one row per sample)
    """
    return df[["sample_id", "participant_id", "biospecimen_id"]].drop_duplicates()


def create_sample_labels(df: pd.DataFrame, style: str = "short") -> dict:
    """
    Create display labels for samples.

    Parameters
    ----------
    df : pd.DataFrame
        Output from load_cell_densities()
    style : str
        'short': participant number only (e.g., '919')
        'medium': participant_id (e.g., 'HTA7_919')
        'full': biospecimen_id (e.g., 'HTA7_919_6')

    Returns
    -------
    dict
        Mapping from sample_id to display label
    """
    meta = get_sample_metadata(df)
    labels = {}

    for _, row in meta.iterrows():
        if style == "short":
            # Extract just the number (e.g., 'HTA7_919' -> '919')
            labels[row["sample_id"]] = row["participant_id"].replace("HTA7_", "")
        elif style == "medium":
            labels[row["sample_id"]] = row["participant_id"]
        else:  # full
            labels[row["sample_id"]] = row["biospecimen_id"]

    return labels


if __name__ == "__main__":
    # Test the loading functions
    print("Loading conic data...")
    conic = load_cell_densities("conic")
    print(f"  Loaded {len(conic)} rows")
    print(f"  Unique samples (files): {conic['sample_id'].nunique()}")
    print(f"  Unique participants: {conic['participant_id'].nunique()}")
    print(f"  Cell classes: {conic['cell_class'].unique().tolist()}")

    print("\nLoading consep data...")
    consep = load_cell_densities("consep")
    print(f"  Loaded {len(consep)} rows")
    print(f"  Unique samples (files): {consep['sample_id'].nunique()}")
    print(f"  Unique participants: {consep['participant_id'].nunique()}")
    print(f"  Cell classes: {consep['cell_class'].unique().tolist()}")

    print("\nPivoted conic data shape (by sample):", pivot_for_stacked_bar(conic).shape)
    print("Pivoted consep data shape (by sample):", pivot_for_stacked_bar(consep).shape)

    print("\nSample metadata:")
    print(get_sample_metadata(consep).to_string())
