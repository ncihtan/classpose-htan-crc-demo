"""Upload the salvaged ClassPose outputs to Synapse (one-off, run 2026-10-15).

Reads from the old working folder (SOURCE), never writes to it. Zarr stores are
zipped one at a time into a staging folder, uploaded, then deleted. Skips files
already in Synapse with the same name and MD5, so it is safe to rerun.

Writes data/synapse_manifest.tsv.

    uv run --no-project --with synapseclient python data/upload_to_synapse.py \
        --source ~/Documents/learn-classpose --staging /tmp/classpose_staging
"""
import argparse
import csv
import hashlib
import os
import shutil
import sys
import zipfile
from pathlib import Path

import synapseclient
from synapseclient import File, Folder

PROJECT_ID = "syn77984172"
REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "data" / "synapse_manifest.tsv"
S3_MANIFEST = REPO / "pipeline" / "s3_manifest.tsv"
NF_CLASSPOSE = "https://github.com/adamjtaylor/nf-classpose/tree/dcb045250a6df043ba0b4eb5912d8b016855aeeb"
COLUMNS = ["synapse_id", "folder", "name", "size_bytes", "md5", "source_path",
           "htan_data_file_id", "htan_biospecimen_id", "htan_participant_id", "model", "note"]


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def zip_dir(src, dest):
    """Store-only zip (zarr chunks are already compressed). Paths inside start with the zarr name."""
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for p in sorted(src.rglob("*")):
            if p.is_file() and p.name != ".DS_Store":
                zf.write(p, Path(src.name) / p.relative_to(src))
    return dest


def htan_ids():
    ids = {}
    with open(S3_MANIFEST) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ids[row["sample_id"]] = {k: row[k] for k in
                                     ("htan_data_file_id", "htan_biospecimen_id", "htan_participant_id")}
    return ids


def plan(source):
    """Yield (folder, local_path, source_rel, zip_it, annotations, note)."""
    dm = source / "htan_hms_crc_he_41_dual_mode"
    ref = source / "htan_hms_crc_he_41"
    ids = htan_ids()

    for zarr in sorted(dm.glob("sample_data/*/*_spatialdata.zarr")):
        sample = zarr.parent.name
        model = "consep" if "_consep_" in zarr.name else "conic"
        ann = {**ids.get(sample, {}), "model": model, "file_type": "spatialdata_zarr_zip",
               "pipeline": NF_CLASSPOSE}
        yield "spatialdata_zarr", zarr, True, ann, "nf-classpose SpatialData output, zipped"

    sid = "dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A"
    yield ("reference_slide", ref / f"{sid}.tif", False,
           {**ids.get(sid, {}), "file_type": "image_tif"},
           "VIPS-converted H&E for HTA7_954, read by figure 2 and 3 scripts. Also in S3.")

    yield ("single_slide_conic_run", ref / sid / f"{sid}_spatialdata.zarr", True,
           {**ids.get(sid, {}), "model": "conic", "file_type": "spatialdata_zarr_zip"},
           "Earlier CoNIC-only run on HTA7_954 (2026-01-25)")
    yield ("single_slide_conic_run", ref / sid / f"{sid}_cell_centroids.geojson", False,
           {**ids.get(sid, {}), "model": "conic", "file_type": "cell_centroids"},
           "Earlier CoNIC-only run on HTA7_954 (2026-01-25)")
    for name in ("cell_morphometrics.parquet", "hierarchical_cell_classification.parquet"):
        yield ("single_slide_conic_run", ref / name, False, {**ids.get(sid, {})},
               "Single-slide morphometrics / CoNIC+CoNSeP hierarchical classification")

    yield ("mif_orion", dm / "P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv", False,
           {"htan_biospecimen_id": "HTA7_926_1013", "source": "HTAN HMS Orion Level 4",
            "file_type": "mcmicro_cell_table"},
           "Orion mIF single-cell table used by match_mif.py. Downloaded from HTAN 2026-02-10; exact file ID not recorded")
    yield ("mif_orion", dm / "data/cell_matching/dg_4DFC_06CED89D-6754-40E2-A418-4F2C288119A2_matched_cells.parquet",
           False, {**ids.get("dg_4DFC_06CED89D-6754-40E2-A418-4F2C288119A2", {}), "file_type": "matched_cells"},
           "ClassPose to mIF cell matches (match_cells.py)")

    yield ("pannuke", source / "pannuke/checkpoints/demo/pannuke_demo/pannuke_demo_best", False,
           {"file_type": "model_checkpoint", "dataset": "PanNuke"},
           "Demo-quality ClassPose fine-tune on 2% of PanNuke, 10 epochs; per-class PQ near 0")

    for name in ("spatial_934.png", "figure2_hot_cold_comparison.png", "malignant_mosaic_all_samples.png"):
        yield ("figures_full_res", dm / "figures" / name, False, {"file_type": "figure"},
               "Full-resolution figure; a downsized preview is in git")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--staging", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    source = args.source.expanduser().resolve()
    args.staging.mkdir(parents=True, exist_ok=True)

    syn = synapseclient.login(silent=True)
    folders = {}
    done = {}
    if MANIFEST.exists():
        with open(MANIFEST) as fh:
            done = {(r["folder"], r["name"]): r for r in csv.DictReader(fh, delimiter="\t")}

    rows = list(done.values())
    for folder, path, zip_it, ann, note in plan(source):
        name = path.name + (".zip" if zip_it else "")
        if (folder, name) in done:
            print("skip (in manifest)", folder, name, flush=True)
            continue
        if not path.exists():
            print("MISSING", path, file=sys.stderr, flush=True)
            continue
        if args.dry_run:
            print("would upload", folder, name, flush=True)
            continue
        if folder not in folders:
            folders[folder] = syn.store(Folder(folder, parent=PROJECT_ID)).id
        local = zip_dir(path, args.staging / name) if zip_it else path
        md5 = md5sum(local)
        ent = syn.store(File(str(local), name=name, parent=folders[folder],
                             annotations={k: v for k, v in ann.items() if v}))
        remote_md5 = ent._file_handle["contentMd5"]
        if remote_md5 != md5:
            print("MD5 MISMATCH", name, md5, remote_md5, file=sys.stderr, flush=True)
        row = {"synapse_id": ent.id, "folder": folder, "name": name, "size_bytes": os.path.getsize(local),
               "md5": md5, "source_path": str(path.relative_to(source)),
               "htan_data_file_id": ann.get("htan_data_file_id", ""),
               "htan_biospecimen_id": ann.get("htan_biospecimen_id", ""),
               "htan_participant_id": ann.get("htan_participant_id", ""),
               "model": ann.get("model", ""), "note": note}
        rows.append(row)
        if zip_it:
            local.unlink()
        with open(MANIFEST, "w", newline="") as fh:
            w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print("uploaded", ent.id, folder, name, flush=True)

    shutil.rmtree(args.staging, ignore_errors=True)


if __name__ == "__main__":
    main()
