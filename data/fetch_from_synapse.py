"""Download the analysis inputs from Synapse into analysis/sample_data/.

Uses data/synapse_manifest.tsv. Zarr zips are unpacked to
sample_data/<sample_id>/<sample_id>_<model>_spatialdata.zarr, the layout the
analysis scripts expect. Other files go to the paths noted below.

    uv run --no-project --with synapseclient python data/fetch_from_synapse.py            # zarrs + figure 2/3 inputs
    uv run --no-project --with synapseclient python data/fetch_from_synapse.py --all      # everything in the manifest
"""
import argparse
import csv
import os
import re
import zipfile
from pathlib import Path

import synapseclient

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "data" / "synapse_manifest.tsv"
DEFAULT_FOLDERS = {"spatialdata_zarr", "reference_slide", "mif_orion"}
SAMPLE_ID = re.compile(r"dg_4DFC_[0-9A-F-]{36}")


def target_dir(row, sample_data):
    if row["folder"] in ("spatialdata_zarr", "reference_slide"):
        return sample_data / SAMPLE_ID.match(row["name"]).group(0)
    if row["folder"] == "mif_orion" and row["name"].endswith(".parquet"):
        return REPO / "analysis" / "data" / "cell_matching"
    if row["folder"] == "mif_orion":
        return sample_data
    return REPO / "data" / "synapse" / row["folder"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="also fetch single-slide run, PanNuke checkpoint, full-res figures")
    ap.add_argument("--sample-data", type=Path,
                    default=Path(os.environ.get("SAMPLE_DATA_DIR", REPO / "analysis" / "sample_data")))
    args = ap.parse_args()

    syn = synapseclient.login(silent=True)
    with open(MANIFEST) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        if not args.all and row["folder"] not in DEFAULT_FOLDERS:
            continue
        dest = target_dir(row, args.sample_data)
        dest.mkdir(parents=True, exist_ok=True)
        ent = syn.get(row["synapse_id"], downloadLocation=str(dest), ifcollision="overwrite.local")
        if row["name"].endswith(".zarr.zip"):
            with zipfile.ZipFile(ent.path) as zf:
                zf.extractall(dest)
            Path(ent.path).unlink()
        print(row["synapse_id"], "->", dest)


if __name__ == "__main__":
    main()
