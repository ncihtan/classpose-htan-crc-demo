#!/usr/bin/env python3
"""
Create flat (denormalized) files from HTAN metadata tables.

Relationships:
- participant (1) → (many) biospecimen (via Participant ID)
- biospecimen (1) → (many) file (via Biospecimen)
- file (1) → (1) files_drs (via Data File ID)
"""

import pandas as pd
from pathlib import Path

# Set up paths
metadata_dir = Path(__file__).parent
output_dir = metadata_dir / "flat"
output_dir.mkdir(exist_ok=True)

# Load all tables
print("Loading metadata files...")
participant = pd.read_csv(metadata_dir / "participant.tsv", sep="\t")
biospecimen = pd.read_csv(metadata_dir / "biospecimen.tsv", sep="\t")
file_df = pd.read_csv(metadata_dir / "file.tsv", sep="\t")
files_drs = pd.read_csv(metadata_dir / "files_drs.csv")

print(f"Loaded: participant ({len(participant)} rows), biospecimen ({len(biospecimen)} rows), "
      f"file ({len(file_df)} rows), files_drs ({len(files_drs)} rows)")

# ============================================================================
# 1. PARTICIPANT FLAT FILE
# Already flat, just save with prefix on columns for clarity when joined
# ============================================================================
print("\n1. Creating participant_flat.tsv...")
participant_flat = participant.copy()
participant_flat.columns = [f"participant_{col}" if col != "HTAN Participant ID" else col
                            for col in participant_flat.columns]
participant_flat.to_csv(output_dir / "participant_flat.tsv", sep="\t", index=False)
print(f"   Saved {len(participant_flat)} rows")

# ============================================================================
# 2. BIOSPECIMEN FLAT FILE
# Join biospecimen with participant information
# ============================================================================
print("\n2. Creating biospecimen_flat.tsv...")

# Prepare participant for join (add prefix to avoid column conflicts)
participant_for_join = participant.copy()
participant_for_join.columns = [f"participant_{col}" if col != "HTAN Participant ID" else col
                                 for col in participant_for_join.columns]

# Join biospecimen with participant
biospecimen_flat = biospecimen.merge(
    participant_for_join,
    left_on="Participant ID",
    right_on="HTAN Participant ID",
    how="left"
)

# Drop duplicate ID column
if "HTAN Participant ID" in biospecimen_flat.columns:
    biospecimen_flat = biospecimen_flat.drop(columns=["HTAN Participant ID"])

biospecimen_flat.to_csv(output_dir / "biospecimen_flat.tsv", sep="\t", index=False)
print(f"   Saved {len(biospecimen_flat)} rows")

# ============================================================================
# 3. FILE FLAT FILE
# Join file with biospecimen and participant information
# ============================================================================
print("\n3. Creating file_flat.tsv...")

# Prepare biospecimen for join (add prefix)
biospecimen_for_join = biospecimen.copy()
biospecimen_for_join.columns = [f"biospecimen_{col}" if col != "HTAN Biospecimen ID" else col
                                 for col in biospecimen_for_join.columns]

# Join file with biospecimen
file_flat = file_df.merge(
    biospecimen_for_join,
    left_on="Biospecimen",
    right_on="HTAN Biospecimen ID",
    how="left"
)

# Extract Participant ID from biospecimen for participant join
# The Participant ID is in the biospecimen_Participant ID column after the merge
participant_id_col = "biospecimen_Participant ID"

# Join with participant
file_flat = file_flat.merge(
    participant_for_join,
    left_on=participant_id_col,
    right_on="HTAN Participant ID",
    how="left"
)

# Clean up duplicate columns
cols_to_drop = ["HTAN Biospecimen ID", "HTAN Participant ID"]
for col in cols_to_drop:
    if col in file_flat.columns:
        file_flat = file_flat.drop(columns=[col])

file_flat.to_csv(output_dir / "file_flat.tsv", sep="\t", index=False)
print(f"   Saved {len(file_flat)} rows")

# ============================================================================
# 4. FILES_DRS FLAT FILE
# Join files_drs with file, biospecimen, and participant information
# ============================================================================
print("\n4. Creating files_drs_flat.tsv...")

# First join files_drs with file on data_file_id = Data File ID
file_for_join = file_df.copy()
file_for_join.columns = [f"file_{col}" if col not in ["Data File ID", "Biospecimen"] else col
                          for col in file_for_join.columns]

files_drs_flat = files_drs.merge(
    file_for_join,
    left_on="data_file_id",
    right_on="Data File ID",
    how="left"
)

# Join with biospecimen (using the biospecimen column from files_drs)
files_drs_flat = files_drs_flat.merge(
    biospecimen_for_join,
    left_on="biospecimen",
    right_on="HTAN Biospecimen ID",
    how="left"
)

# Join with participant
files_drs_flat = files_drs_flat.merge(
    participant_for_join,
    left_on="biospecimen_Participant ID",
    right_on="HTAN Participant ID",
    how="left"
)

# Clean up duplicate columns
cols_to_drop = ["Data File ID", "Biospecimen", "HTAN Biospecimen ID", "HTAN Participant ID"]
for col in cols_to_drop:
    if col in files_drs_flat.columns:
        files_drs_flat = files_drs_flat.drop(columns=[col])

files_drs_flat.to_csv(output_dir / "files_drs_flat.tsv", sep="\t", index=False)
print(f"   Saved {len(files_drs_flat)} rows")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*60)
print("SUMMARY - Flat files created in metadata/flat/")
print("="*60)
for f in sorted(output_dir.glob("*.tsv")):
    df = pd.read_csv(f, sep="\t")
    print(f"  {f.name}: {len(df)} rows x {len(df.columns)} columns")

print("\nColumn counts by source table:")
print(f"  participant: {len(participant.columns)} original columns")
print(f"  biospecimen: {len(biospecimen.columns)} original columns")
print(f"  file: {len(file_df.columns)} original columns")
print(f"  files_drs: {len(files_drs.columns)} original columns")
