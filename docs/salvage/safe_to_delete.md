# Safe to delete

This page lists each large item in the old working folder `~/Documents/learn-classpose` (31 GB), with the reason it can be deleted. Nothing has been deleted. An item is safe to delete if it is junk, if it is input data with a known source, or if a verified copy exists in git, Synapse or S3. "Verified" for Synapse means that the MD5 Synapse recorded matches the local MD5 (see `data/synapse_manifest.tsv`).

| Item | Size | Class | Copy | Safe? |
|---|---|---|---|---|
| `htan_hms_crc_he_41_dual_mode/sample_data/*/*_spatialdata.zarr` (42) | 16 GB | COPIED | Synapse `spatialdata_zarr` (zipped, MD5 checked) and S3 | Yes |
| `htan_hms_crc_he_41_dual_mode/sample_data/` tifs and geojsons for `98B41CE1` and `1A504F88` | 7.6 GB | COPIED | S3, listed in `pipeline/s3_manifest.tsv` with matching sizes | Yes |
| `htan_hms_crc_he_41_dual_mode/P37_S29_..._cellRingMask.csv` | 776 MB | INPUT | Synapse `mif_orion` (MD5 `2fcab9b0f173c9bb99cd4c8c31d87785`), and HTAN HMS Orion Level 4 | Yes |
| `htan_hms_crc_he_41_dual_mode/data/cell_matching/*.parquet` | 24 MB | UNIQUE | Synapse `mif_orion` | Yes |
| `htan_hms_crc_he_41_dual_mode/figures/` | 135 MB | UNIQUE | git `analysis/figures/`, plus 3 full-size PNGs in Synapse `figures_full_res` | Yes |
| `htan_hms_crc_he_41_dual_mode/` code, docs, metadata, data CSVs, Nextflow reports | 40 MB | UNIQUE | git | Yes |
| `htan_hms_crc_he_41/dg_4DFC_6C434FEA-....tif` | 1.9 GB | COPIED | Synapse `reference_slide` and S3 | Yes |
| `htan_hms_crc_he_41/dual_model_output/` | 2.4 GB | COPIED | S3 (the same sample's dual-model outputs). The zarr is also in Synapse. | Yes |
| `htan_hms_crc_he_41/dg_4DFC_6C434FEA-.../` (CoNIC-only run) | 802 MB | UNIQUE | Synapse `single_slide_conic_run` (zarr zipped, centroids geojson) | Yes |
| `htan_hms_crc_he_41/*.parquet` | 96 MB | UNIQUE | Synapse `single_slide_conic_run` | Yes |
| `htan_hms_crc_he_41/` scripts, TLS methods, small figures | 20 MB | UNIQUE | git `archive/single-slide-conic/`, `docs/tls_analysis_methods.md` | Yes. The large figures are listed in `what-we-did-not-retain.md`. |
| `pannuke/data/` | 3.4 GB | INPUT | Hugging Face `RationAI/PanNuke`, `gs://htan-dcc-pannuke/data/` | Yes |
| `pannuke/checkpoints/` | 1.1 GB | UNIQUE | Synapse `pannuke` (MD5 `c04093e8d20beb732eb786c33cca8566`) | Yes |
| `pannuke/` code and notes | under 1 MB | UNIQUE | git `pannuke/` | Yes |
| `htan1_examples/` | 1.3 GB | INPUT | HTAN Synapse IDs in `data/htan1_examples_manifest.csv` | Yes |
| `classpose/` | 60 MB | INPUT | GitHub `sohmandal/classpose` at `617508c`, weights on Hugging Face. Local scripts are in `patches/`. | Yes |
| `nf-classpose/` | 25 MB | COPIED | GitHub `adamjtaylor/nf-classpose`; example run record in `pipeline/example_run/` | Yes. Revoke the Seqera token in `example_run/logs/` first. |
| `output/` | 9 MB | JUNK | Test outputs on a public slide | Yes |
| `first_plan.txt`, `htan_study_plan.md` | 108 KB | UNIQUE | `docs/sessions.md` (summary), git `docs/htan_study_plan.md` | Yes |
| `.claude/`, `.DS_Store`, `__pycache__/` | small | JUNK | | Yes |

Total: about 31 GB, all of it safe to delete once the upload below is complete.

## Upload status

Filled in after the Synapse upload finished.
