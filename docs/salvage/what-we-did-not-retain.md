# What we did not retain

The original working folder was `~/Documents/learn-classpose` on Adam's laptop (31 GB, wiped 2026-10-16). This page lists what did not come into this repo or Synapse, and why. If you need one of these items, the "Where it lives" column says how to get it back. The full file list is in `file_inventory.tsv`, and sizes and classes are in `size_inventory.tsv`.

## Not copied anywhere

| Item | Size | Why we dropped it | Where it lives |
|---|---|---|---|
| `nf-classpose/example_run/logs/` | under 1 MB | It holds live Seqera access and refresh tokens | Nowhere. Revoke the token in Seqera. |
| `nf-classpose/.nextflow.log*`, `.nextflow/`, `results/` | 16 MB | Nextflow logs and local test output | Nowhere. Rerun if needed. |
| `nf-classpose/example_run/dg_4DFC_624693b0-.../` geojsons and zarr | about 5 MB | Test-run outputs for one slide | Rerun the example run. The density CSV and reports are in `pipeline/example_run/`. |
| `first_plan.txt` | 100 KB | A raw Claude Code transcript, which may hold secrets | Summarised in `docs/sessions.md` |
| Past Claude Code transcripts | none left | Claude Code had already deleted them | `docs/sessions.md` is rebuilt from the prompt history |
| `htan_hms_crc_he_41/dg_4DFC_*/` (41 folders of CoNIC-only density CSVs) | about 200 KB | Superseded by the dual-model run | Rerun nf-classpose with `--models conic` |
| `htan_hms_crc_he_41/` figures over 3 MB (`tls_revised_examples.png`, `hierarchical_classification_figure.{png,pdf}`, `cell_types_crops.png`, `cell_types_spatial_plot.png`, `hierarchical_classification_spatial.png`, `tls_graph_clusters.png`) and two stray copies (`" .png"`, `cell_type_stacked_bars copy.png`) | about 55 MB | Large, and the scripts that make them are in `archive/single-slide-conic/` | Regenerate from the archive scripts and the Synapse `single_slide_conic_run` folder |
| `htan_hms_crc_he_41_dual_mode/figures/figure2d_cell_gallery.png.zip` | 1 MB | A zipped duplicate of a figure in git | `analysis/figures/figure2d_cell_gallery.png` |
| `htan_hms_crc_he_41_dual_mode/main.py` | 106 B | The `uv init` stub | Nothing to keep |
| `output/` | 9 MB | ClassPose outputs on the public CMU-1-Small-Region test slide | Rerun ClassPose on the OpenSlide test slide |
| `.claude/settings.local.json` (5 copies) | small | Local Claude Code permissions | Nowhere |
| `.DS_Store`, `__pycache__/`, `*.egg-info`, `pannuke/dist/` | small | Junk and build artefacts | Nowhere |
| The original nf-classpose params file for the main run | small | It was a temporary Seqera upload and was never saved | `pipeline/params.reconstructed.yaml` is our best reconstruction |

## Kept outside git

| Item | Size | Where it lives |
|---|---|---|
| nf-classpose outputs: VIPS tifs, cell contour and centroid geojsons, tissue contours, `versions.yml`, reports | 234 GB | `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/`, listed in `pipeline/s3_manifest.tsv` |
| SpatialData zarrs (41 CoNSeP and 1 CoNIC) | about 16 GB | Synapse `syn77984172/spatialdata_zarr` and S3 |
| HTA7_954 slide used by figures 2 and 3 | 1.9 GB | Synapse `reference_slide` and S3 |
| Single-slide CoNIC run, morphometrics and hierarchical classification parquets | about 1 GB | Synapse `single_slide_conic_run` |
| Orion mIF cell table for HTA7_926_1013 and the matched cells | 800 MB | Synapse `mif_orion`. The table is HTAN HMS Orion Level 4 data. Its HTAN file ID was not recorded. |
| PanNuke demo checkpoint | 1.1 GB | Synapse `pannuke` |
| Three full-resolution figures | 70 MB | Synapse `figures_full_res`. Previews are in `analysis/figures/`. |
| Input data: 41 HTAN H&E OME-TIFFs | not stored locally | CRDC by DRS URI (`analysis/metadata/files_drs.csv`) |
| Input data: `htan1_examples/` (6 HTAN1 H&E OME-TIFFs) | 1.3 GB | HTAN on Synapse. The Synapse IDs are in `data/htan1_examples_manifest.csv`, and the MD5s matched. |
| Input data: `pannuke/data/` | 3.4 GB | Hugging Face `RationAI/PanNuke` and `gs://htan-dcc-pannuke/data/` |
| Input data: `classpose/models/` | 50 MB | ClassPose weights on Hugging Face (`classpose/classpose`) and Zenodo |
| Code: `classpose/` clone | | `github.com/sohmandal/classpose` at `617508c`. The three local scripts are in `patches/classpose-local/`. |
| Code: `nf-classpose/` | | `github.com/adamjtaylor/nf-classpose`. The run commit `dcb0452` is merged to `main` as `ce699b6`. |

## Known gaps

- Six entries in `docs/figure_map.tsv` have no matching script (`unknown`).
- Three plot scripts read the v1 immune phenotype CSV (`immune_phenotype_classification.csv`), while the docs use v2. Both CSVs are in `analysis/data/`.
- The S3 output listing was last checked in January 2026 (`pipeline/s3_manifest.tsv`). We could not check it again during the salvage, because the AWS login had expired.
