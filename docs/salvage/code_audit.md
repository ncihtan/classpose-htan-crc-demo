# Code audit: learn-classpose (read-only salvage)

- SOURCE: `~/Documents/learn-classpose`.
- Main analysis: `htan_hms_crc_he_41_dual_mode/`, abbreviated `DM/` below.
- Paths in this document are relative to `DM/`, which becomes `analysis/` in the demo repo, unless shown as absolute.
- Audit date: 2026-10-15.
- Nothing in SOURCE was modified.

Abbreviations used in the tables:

| Abbrev | Path |
|---|---|
| `CPO[m]` | `classpose_outputs/data/<participant>/<biospecimen>/<m>/*_cell_densities.csv`. `m` is conic, consep or both. Each file is about 290 B; 82 files in all. Small enough for git. |
| `ZPTS` | `sample_data/<id>/<id>_consep_spatialdata.zarr/points/cell_centroids/points.parquet`. These are read directly with pandas. |
| `ZARR` | The full `sample_data/<id>/<id>_consep_spatialdata.zarr`, read with `spatialdata.read_zarr`. |
| `TIF` | `sample_data/<id>/<id>.tif`, the vips-converted H&E (about 2 GB each). |
| `REF` | `~/Documents/learn-classpose/htan_hms_crc_he_41`, the older single-slide run. |

---

## 1. Script inventory

### 1a. Python analysis scripts (`scripts/`)

| Script | Purpose | Reads | Writes | Absolute paths / S3 / profile |
|---|---|---|---|---|
| `load_cell_densities.py` | Shared loader: densities CSVs to long DataFrame, plus sample/participant helpers. Defines `get_project_root()` as `Path(__file__).parent.parent` (line 16). | `CPO` via rglob (line 34-38) | - | none |
| `utils/__init__.py` | Re-exports spatial_utils | - | - | none |
| `utils/spatial_utils.py` | Shared constants: `UM_PER_PIXEL=0.325`, ConSeP colours, scale bar, contours, legend. Its own `get_project_root()` is `parents x3` (line 47). | (`get_sample_path_info` builds REF paths; it is never called) | - | **L68 `REF`** |
| `select_representative_samples.py` | Ranks the 41 slides by inflammatory %, picks hot/cold/prototype | `CPO[consep]`, `metadata/files_drs.csv`; checks REF tif/zarr exist (print only) | `data/sample_immune_ranking.csv`, `data/selected_samples.csv` | **L139 `REF`** (also used at L156-157) |
| `classify_immune_phenotypes.py` | v1 grid-based Immunoscore-like phenotype (desert / excluded / immunosuppressed / hot) | `ZPTS` (all 41), `CPO[consep]` | `data/immune_phenotype_metrics.csv`, `data/immune_phenotype_classification.csv`, `figures/immune_phenotype_{distribution,scatter}.png` | none |
| `classify_immune_phenotypes_v2.py` | v2 KNN-graph microenvironment phenotype with absolute density thresholds | `ZPTS` (all 41), `CPO[consep]` | `data/immune_phenotype_{metrics,classification}_v2.csv`, `figures/immune_phenotype_{distribution,scatter}_v2.png` | none |
| `compute_spatial_features.py` | NN immune-tumour distance, DBSCAN/HDBSCAN cluster features, Moran's I, all per participant | `ZPTS` (all 41), `CPO[consep]` | `data/spatial_features.csv` | none |
| `survival_analysis.py` | Builds the survival dataset, then runs univariate and stage-adjusted Cox, KM (median split) and a forest plot | `data/immune_phenotype_metrics_v2.csv`, `data/immune_phenotype_classification_v2.csv`, `data/spatial_features.csv`, `CPO[both]`, `metadata/participant.tsv` | `data/survival_analysis_dataset.csv`, `data/survival_cox_results.csv`, `data/survival_cox_multivariate_results.csv`, `figures/km_<feature>_{os,css,pfs}.png`, `figures/km_phenotype_*.png`, `figures/survival_forest_plot.png` | none |
| `dbscan_sensitivity.py` | DBSCAN eps x min_samples grid vs Cox HR/C-index | `data/survival_analysis_dataset.csv`, `CPO[consep]`, `ZPTS` (all 41) | `data/dbscan_sensitivity.csv`, `figures/dbscan_sensitivity.png` | none |
| `mutation_interaction.py` | Tests whether mutation status (KRAS/BRAF/TP53/APC) changes the prognostic value of immune clustering (`db_tumor_cluster_fraction`) | `data/survival_analysis_dataset.csv`, `metadata/bquxjob_6358239c_19c48c9d3ac.csv` (MolecularTest, L48), `metadata/bquxjob_4136c6be_19c48c7c3e6.csv` (Therapy, L72) | `data/mutation_{univariate,stratified,interaction,burden}_cox.csv`, `figures/mutation_forest_plot.png`, `figures/km_interaction_<gene>_<ep>.png` | none |
| `match_cells.py` | CoNIC vs ConSeP cell matching (centroid kNN, optional contour IoU). The only script with a CLI (`--sample`/`--all`). | `data/sample_manifest.txt`, `sample_data/<id>/<id>_{conic,consep}_spatialdata.zarr/points/...` and `shapes/...`. Missing parquets are auto-downloaded from S3. | `data/cell_matching/<id>_matched_cells.parquet`, `data/cell_matching/matching_summary.csv` | **L37-38** bucket `htan-project-tower-bucket` and prefix `nf-classpose/outputs/htan_hms_crc_he_41_dual_model`; **L39** `AWS_PROFILE="tower"`; **L114** `aws s3 cp --profile` |
| `match_mif.py` | Matches ClassPose centroids to Orion mIF CellRingMask (one slide, P37/S29); marker validation plots | `P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv` at project root (813 MB, L42); `ZPTS` for `dg_4DFC_CD2C5DCB-...` | `figures/match_distance_diagnostic.png`, `match_area_correlation.png`, `marker_boxplots_by_class.png`, `marker_violins_by_class.png`, `cd45_scatter.png` | none, but the 813 MB file is expected at the repo root |
| `plot_classpose_runtimes.py` | Runtimes of the CLASSPOSE_PREDICT_WSI process: raw, per 1k cells, per GB read | `nextflow_reports/execution_trace_*.txt`, `CPO[both]` | `figures/classpose_runtimes.png` | none |
| `plot_clustering_clr.py` | CLR-transformed hierarchical clustering, plus raw vs CLR correlation | `CPO[both]` | `figures/clustering_clr_{conic,consep}.png`, `figures/correlation_comparison_raw_vs_clr.png` | none |
| `plot_clustering_conic.py` / `plot_clustering_consep.py` | Proportional-composition clustermap | `CPO[conic]` / `CPO[consep]` | `figures/clustering_{conic,consep}.png` | none |
| `plot_scatter_conic_vs_consep.py` | Total cells per slide, CoNIC vs ConSeP | `CPO[both]` | `figures/scatter_conic_vs_consep.png` | none |
| `plot_stacked_bar_{conic,consep}.py` | Stacked bar of counts by type | `CPO[m]` | `figures/stacked_bar_{m}.png` | none |
| `plot_stacked_100_bar_{conic,consep}.py` | Stacked bar of proportions (100%) by type | `CPO[m]` | `figures/stacked_100_bar_{m}.png` | none |
| `plot_figure1.py` | Figure 1 panels: cohort heatmap, validation, heterogeneity, serial sections | `CPO[both]`, `metadata/participant.tsv`, `metadata/files_drs.csv` | `figures/figure1{a_cohort,b_validation,c_heterogeneity,d_serial_sections}.{png,pdf}` | none |
| `plot_figure2_hot_cold_comparison.py` | Hot (HTA7_917) vs cold (HTA7_932) overview and zooms with contours | `TIF` + `ZARR` (shapes) for `1A504F88`, `98B41CE1` | `figures/figure2_hot_cold_comparison.png` | none |
| `plot_figure2a_slide_overview.py` | HTA7_954 slide thumbnail with ROIs | REF `<6C434FEA>.tif` | `figures/figure2a_slide_overview.{png,pdf}` | **L26 `REF`** |
| `plot_figure2b_region_immune.py` | Immune-rich ROI with contours | REF tif + REF `dual_model_output/consep/<id>_consep_spatialdata.zarr` | `figures/figure2b_region_immune.{png,pdf}` | **L28 `REF`** |
| `plot_figure2c_region_epithelial.py` | Glandular ROI with contours | REF tif + zarr | `figures/figure2c_region_epithelial.{png,pdf}` | **L29 `REF`** |
| `plot_figure2d_cell_gallery.py` | Per-class single-cell crops | REF tif + zarr | `figures/figure2d_cell_gallery.{png,pdf}` | **L30 `REF`** |
| `plot_figure2e_regional_comparison.py` | Composition of ROI B vs ROI C | REF zarr | `figures/figure2e_regional_comparison.{png,pdf}` | **L28 `REF`** |
| `plot_figure3a_spatial_overview.py` | Whole-slide centroid map by class | REF tif + zarr | `figures/figure3a_spatial_overview.{png,pdf}` | **L33 `REF`** |
| `plot_figure3b_density_heatmaps.py` | KDE per cell type | REF tif + zarr | `figures/figure3b_density_heatmaps.{png,pdf}` | **L30 `REF`** |
| `plot_immune_compartment_bars.py` | Immune % in core/TSI/stroma/other per slide | `data/immune_phenotype_metrics_v2.csv`, `ZPTS` (all) | `figures/immune_compartment_stacked_bars.png`, `data/immune_compartment_distribution.csv` | none |
| `plot_knn_regions.py` | KNN region map, hot vs cold | `ZPTS` for `B0A48536`, `8C77314D` | `figures/knn_region_comparison.png` | none |
| `plot_malignant_mosaic.py` | Grid of all 41 slides coloured by class, sorted by phenotype | `ZARR` (all 41), `data/immune_phenotype_classification.csv` (v1) | `figures/malignant_mosaic_all_samples.{png,pdf}` | none |
| `plot_phenotype_paired_panels.py` | 4x2 grid: cell types vs KNN regions per phenotype | `data/immune_phenotype_metrics_v2.csv`, `ZPTS` for 4 hardcoded IDs | `figures/phenotype_paired_panels.png` | none |
| `plot_phenotype_spatial_maps.py` | Filled-contour density maps per phenotype | `data/immune_phenotype_classification.csv` (v1), `ZPTS` for 4 IDs | `figures/immune_phenotype_spatial_comparison.png`, `figures/spatial_map_<phenotype>.png` | none |
| `plot_region_overlay.py` | Region heatmap over greyscale H&E | `data/immune_phenotype_classification.csv` (v1), `ZPTS` + **`TIF`** for `B0A48536`, `96A5875F`, `1A504F88`, `8C77314D` | `figures/region_assignment_comparison.png`, `figures/region_overlay_<phenotype>.png` | none. **Only `1A504F88.tif` is present locally**, so the other 3 tifs must be fetched (`vips_converted/`). |
| `survival_analysis.md` | Write-up (doc, not code) | - | - | none |

### 1b. Shell scripts (`scripts/`)

See section 3 for the S3 details.

| Script | Purpose | Reads | Writes | Absolute paths / S3 / profile |
|---|---|---|---|---|
| `generate_s3_manifest.sh` | `aws s3 ls` of the run outdir joined to HTAN IDs | `metadata/files_drs.csv`, S3 listing | `classpose_outputs/s3_manifest.tsv` | **L9** profile `tower`; **L10** bucket path; **L180** hardcoded `s3://htan-project-tower-bucket/` |
| `fetch_from_s3.sh` | Filtered download driven by the manifest | `classpose_outputs/s3_manifest.tsv` | `classpose_outputs/data/<participant>/<biospecimen>/<model>/` | **L16** profile `tower`; **L281** `--profile` |
| `fetch_all_spatialdata.sh` | Downloads the consep zarr for every ID in the manifest | `data/sample_manifest.txt` (cwd-relative) | `sample_data/<id>/<id>_consep_spatialdata.zarr/` | **L10** bucket path; no profile flag |
| `fetch_hot_cold_samples.sh` | Downloads tif + consep geojsons for 2 samples | S3 | `sample_data/<id>/` | **L6 `~/.../sample_data`**; L23/28/33/46/51/56 S3 URIs |

### 1c. Root and other code

| File | Purpose | Reads | Writes | Absolute paths / S3 / profile |
|---|---|---|---|---|
| `build_slides.py` | Builds Google Slides from `classpose_meeting_deck.md` with the `gws` CLI; uploads 12 figures to Drive | `classpose_meeting_deck.md`, `figures/*.png` | Google Drive folder + Slides (remote) | No paths. **L132-136 grants every uploaded image `{"role":"reader","type":"anyone"}`** (public link). Needs `gws` auth. Not needed in the demo. |
| `main.py` | uv placeholder ("Hello from ...") | - | - | none; drop it |
| `metadata/create_flat_files.py` | Denormalises the HTAN participant/biospecimen/file/drs tables | `metadata/{participant,biospecimen,file}.tsv`, `metadata/files_drs.csv` | `metadata/flat/*_flat.tsv` | none (`Path(__file__).parent`) |

### 1d. Full grep result

Command: `grep -rn "/Users/\|s3://\|syn[0-9]\{6,\}\|profile\|expanduser\|home()\|[0-9]\{12\}\|~/"` over scripts/, build_slides.py and main.py.

- **`/Users/ataylor`:** 10 hits in 10 files:
  - plot_figure2a:26, 2b:28, 2c:29, 2d:30, 2e:28, 3a:33, 3b:30
  - select_representative_samples:139
  - utils/spatial_utils:68
  - fetch_hot_cold_samples.sh:6
- **`s3://`:** match_cells:48 (built from L37-38), fetch_all_spatialdata:10, fetch_hot_cold_samples ×6, generate_s3_manifest:10 and 180.
- **Docs and data that will be copied:**
  - `methods.md:36` names the bucket.
  - `classpose_outputs/s3_manifest.tsv` has the bucket URI on every one of its 460 rows.
- **AWS profile:** `tower` only (match_cells:39, fetch_from_s3:16, generate_s3_manifest:9).
- **No hits:**
  - Synapse IDs, AWS account IDs, `expanduser`/`home()` and `~/` in any script.
  - `/Users/` in the md docs.
- **Leftovers that should not be copied:** `scripts/__pycache__/`, `.DS_Store` files, `DM/.claude/settings.local.json` (local permission allowlist).

### 1e. Execution order (as implied by the data dependencies)

1. `generate_s3_manifest.sh` and then `fetch_from_s3.sh --file-type cell_densities --all` produce `classpose_outputs/`. This output is already small, so it can just be committed.
2. `fetch_all_spatialdata.sh` fetches the 41 consep zarrs (about 19 GB with the tifs). `fetch_hot_cold_samples.sh` fetches 2 tifs.
3. `select_representative_samples.py`
4. `classify_immune_phenotypes.py` (v1). Its output is still read by the mosaic, spatial_maps and region_overlay scripts.
5. `classify_immune_phenotypes_v2.py`
6. `compute_spatial_features.py`
7. `survival_analysis.py`
8. `dbscan_sensitivity.py` and `mutation_interaction.py`
9. All `plot_*.py`, `match_cells.py` and `match_mif.py`, in any order.

---

## 2. Path resolution and the minimal config proposal

### 2a. Current behaviour

| Mechanism | Where | Resolves to |
|---|---|---|
| `load_cell_densities.get_project_root()` = `Path(__file__).parent.parent` | Most analysis and plot scripts (23 import `load_cell_densities`) | DM root |
| `utils.spatial_utils.get_project_root()` = `Path(__file__).parent.parent.parent` | plot_figure2*, plot_figure3*, hot_cold | DM root |
| `Path(__file__).resolve().parents[1]` | match_mif.py:40, plot_classpose_runtimes.py:26 | DM root |
| `Path(__file__).parent` | build_slides.py:13, metadata/create_flat_files.py:15 | DM root / metadata |
| `SCRIPT_DIR=$(dirname BASH_SOURCE)` then `PROJECT_DIR=..` | fetch_from_s3.sh:17-18, generate_s3_manifest.sh:11-12 | DM root |
| **cwd-relative** (`sample_data`, `data/sample_manifest.txt`) | fetch_all_spatialdata.sh:11-12 | Must be run from the DM root |
| **Absolute** | the 10 `/Users/ataylor` hits above | breaks anywhere else |

- Imports work because each script does `sys.path.insert(0, Path(__file__).parent)`, which puts `scripts/` on the path.
- There is no config file and no environment variables. Only match_cells.py has argparse.
- **Because `analysis/` mirrors the DM root, every file-relative path keeps working after the copy without edits.**
- Only these break: the absolute `REF` paths, the root-level mIF CSV, and the hardcoded bucket and profile.

### 2b. Proposal: one file, `analysis/scripts/config.py`

The file goes in `scripts/` rather than `analysis/` because every script already puts `scripts/` on `sys.path`. Both `from config import ...` and `utils/spatial_utils.py` then work with no extra path hacks.

```python
"""Single source of truth for analysis paths. Every value is env-var overridable."""
import os
from pathlib import Path

def _p(var: str, default: Path) -> Path:
    return Path(os.environ.get(var, default)).expanduser().resolve()

ROOT            = _p("CLASSPOSE_ROOT", Path(__file__).resolve().parents[1])   # analysis/
SAMPLE_DATA_DIR = _p("SAMPLE_DATA_DIR", ROOT / "sample_data")   # zarrs + tifs + mIF csv (Synapse/S3, gitignored)
METADATA_DIR    = _p("METADATA_DIR", ROOT / "metadata")
DATA_DIR        = _p("DATA_DIR", ROOT / "data")
FIGURES_DIR     = _p("FIGURES_DIR", ROOT / "figures")
MIF_CSV         = _p("MIF_CSV", SAMPLE_DATA_DIR / "P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv")
S3_BASE         = os.environ.get("CLASSPOSE_S3_BASE",
                                 "s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model")
AWS_PROFILE     = os.environ.get("AWS_PROFILE")   # None -> default credential chain
```

The `S3_BASE` default is the current private bucket. Swap it for the public Synapse/S3 location once that exists.

### 2c. Required edits (Tier 1: without these the code breaks on any other machine)

In each Python file, add `from config import ...` just after the existing `sys.path.insert` line.

| File | Line(s) | Change |
|---|---|---|
| `scripts/plot_figure2a_slide_overview.py` | 26-27 | Delete `REF_PROJECT`. Set `HE_TIF_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}.tif"` |
| `scripts/plot_figure2b_region_immune.py` | 28-30 | Same for `HE_TIF_PATH`. Set `CONSEP_ZARR_PATH = SAMPLE_DATA_DIR / SAMPLE_ID / f"{SAMPLE_ID}_consep_spatialdata.zarr"` |
| `scripts/plot_figure2c_region_epithelial.py` | 29-31 | as 2b |
| `scripts/plot_figure2d_cell_gallery.py` | 30-32 | as 2b |
| `scripts/plot_figure2e_regional_comparison.py` | 28-29 | `CONSEP_ZARR_PATH` only |
| `scripts/plot_figure3a_spatial_overview.py` | 33-40 | as 2b (the multi-line `CONSEP_ZARR_PATH`) |
| `scripts/plot_figure3b_density_heatmaps.py` | 30-32 | as 2b |
| `scripts/select_representative_samples.py` | 139, 156-157 | Change `ref_project` to `SAMPLE_DATA_DIR / prototype_sample`, using the same filenames as above |
| `scripts/utils/spatial_utils.py` | 68-73 | `get_sample_path_info()` is dead code. Point it at `SAMPLE_DATA_DIR` the same way, or delete it together with `utils/__init__.py:9,18` |
| `scripts/match_mif.py` | 42 | `MIF_CSV = config.MIF_CSV`. The 813 MB file cannot live at the repo root in git. |
| `scripts/match_cells.py` | 37-39, 48, 114 | Use `config.S3_BASE` and `config.AWS_PROFILE`. Only append `--profile` when it is set. |
| `scripts/fetch_hot_cold_samples.sh` | 6; 23, 28, 33, 46, 51, 56 | `OUTPUT_DIR="${SAMPLE_DATA_DIR:-$PROJECT_DIR/sample_data}"`, with SCRIPT_DIR/PROJECT_DIR as in fetch_from_s3.sh:17-18. Use `S3_BASE="${CLASSPOSE_S3_BASE:-...}"` and build the URIs from it. |
| `scripts/fetch_all_spatialdata.sh` | 10-12 | Same `S3_BASE` env default. Make `OUTPUT_DIR` and `MANIFEST` SCRIPT_DIR-relative (they are cwd-relative today). |
| `scripts/fetch_from_s3.sh` | 16, 281 | `AWS_PROFILE="${AWS_PROFILE:-}"`. Pass `--profile` only if it is non-empty. |
| `scripts/generate_s3_manifest.sh` | 9-10, 84, 180 | Same profile handling. **L180 hardcodes `s3://htan-project-tower-bucket/` and ignores `--bucket` (bug).** Derive it from `$S3_BUCKET`. |

**Data move needed for Tier 1:** add `dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A.tif` (HTA7_954, 2.0 GB, currently only in `REF/`) to `sample_data/dg_4DFC_6C434FEA-.../` on Synapse/S3. The consep zarr is already in `DM/sample_data`; I confirmed it is the same size and that `points.parquet` has the same md5.

### 2d. Optional edits

**Tier 2: only needed if `SAMPLE_DATA_DIR` must be overridable by env var.**
- If data is placed or symlinked at `analysis/sample_data`, these lines need no change.
- Otherwise, replace `get_project_root() / "sample_data"` with `SAMPLE_DATA_DIR` at:
  - `classify_immune_phenotypes.py:448`
  - `classify_immune_phenotypes_v2.py:99, 348`
  - `match_cells.py:57, 154, 186`
  - `match_mif.py:47`
  - `plot_figure2_hot_cold_comparison.py:46`
  - `plot_immune_compartment_bars.py:52`
  - `plot_knn_regions.py:55`
  - `plot_malignant_mosaic.py:28`
  - `plot_phenotype_paired_panels.py:69`
  - `plot_phenotype_spatial_maps.py:66`
  - `plot_region_overlay.py:62, 76`
- That is 13 files and 17 lines. `figures/` (33 call sites), `data/` and `metadata/` can stay root-relative; they are in git.

**Tier 3 (cosmetic):** rename the two BigQuery job CSVs to `metadata/htan_molecular_test.csv` and `metadata/htan_therapy.csv`, and update `mutation_interaction.py:48, 72`.

---

## 3. S3 fetch scripts

All of them use one bucket and one profile:
- Bucket and prefix: `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/`
- AWS profile: `tower`

| Script | S3 objects | Profile | Path handling | Notes |
|---|---|---|---|---|
| `generate_s3_manifest.sh` | `aws s3 ls <base>/ --recursive` (L84); URI rebuilt at L180 | `${AWS_PROFILE:-tower}`, `--profile` flag | SCRIPT_DIR-relative | Skips `.zarr/` internals. L180 hardcodes the bucket. |
| `fetch_from_s3.sh` | Takes the `s3_uri` column from `classpose_outputs/s3_manifest.tsv`: `<base>/<id>/{conic,consep,vips_converted}/...` | `${AWS_PROFILE:-tower}`, `--profile` flag (L281) | SCRIPT_DIR-relative | Filters by sample, participant, biospecimen, model and file-type |
| `fetch_all_spatialdata.sh` | `<base>/<id>/consep/<id>_consep_spatialdata.zarr/` (recursive) | **none passed**: uses `$AWS_PROFILE` or the default chain | **cwd-relative** | Header says 42 samples; the manifest has 41 |
| `fetch_hot_cold_samples.sh` | `<base>/<id>/vips_converted/<id>.tif`, `<base>/<id>/consep/<id>_consep_cell_{contours,centroids}.geojson` for `1A504F88` (HTA7_917) and `98B41CE1` (HTA7_932) | **none passed** | **absolute** `~/...` | Fetches the geojsons, but the plotting script reads the zarr (from fetch_all_spatialdata) |
| `match_cells.py` (Python) | `<base>/<id>/<model>/<id>_<model>_spatialdata.zarr/{points/cell_centroids/points.parquet/, shapes/cell_contours/shapes.parquet}` | hardcoded `tower` | project-root-relative | Auto-fetches on demand |

The Nextflow reports also expose the scratch bucket `s3://htan-project-tower-scratch/work/` (section 6).

---

## 4. Figure map

See `docs/figure_map.tsv`:
- 73 rows covering all 124 files in `DM/figures/`.
- The km family is grouped into 4 rows: 16 features x 3 endpoints, phenotype, stale largest_cluster_size, interaction.
- Script paths are given as `analysis/scripts/...`. Input paths are relative to `analysis/`.

| Status | Figures |
|---|---|
| No producing script (`unknown`) | `conic_consep_confusion_matrix.png` (probably ad-hoc from `data/cell_matching/*_matched_cells.parquet`), `epithelial_inflammatory_investigation.png`, `spatial_934.png` (40 MB), `figure1.{png,pdf}` (earlier combined version of plot_figure1.py), `km_largest_cluster_size_{os,pfs,css}.png` (stale: feature dropped from `FEATURES`) |
| Not figures | `figure1_caption.md`, `immune_phenotype_figures.md` (docs), `figure2d_cell_gallery.png.zip` (manual zip) |
| Need data not currently local | `region_overlay_*` and `region_assignment_comparison.png` need tifs for `B0A48536`, `96A5875F`, `8C77314D`. figure2a-d and 3a-b need the HTA7_954 tif from `REF/`. |
| Note | `tumor_grade_high` is in `FEATURES` but has no KM plot. Its median is 0, so the low group is empty and the plot is skipped. |
| Used in the deck (`classpose_meeting_deck.md`) | classpose_runtimes, figure1a-d, figure2d_cell_gallery, immune_phenotype_scatter_v2, immune_phenotype_spatial_comparison, km_interaction_KRAS_os, malignant_mosaic_all_samples, marker_violins_by_class, survival_forest_plot |

---

## 5. data/ and metadata/

### 5a. data/ (24 MB; everything here is derived, nothing is raw)

| File | Size | Rows | What | Level |
|---|---|---|---|---|
| `sample_manifest.txt` | 1.8 KB | 41 | Slide IDs (dg_4DFC GUIDs) | IDs |
| `sample_immune_ranking.csv` | 6.9 KB | 41 | Inflammatory counts/% + filename + HTAN data_file_id | slide |
| `selected_samples.csv` | 645 B | 3 | Hot/cold/prototype picks | slide |
| `immune_phenotype_classification{,_v2}.csv` | 4.7 / 4.5 KB | 41 | Phenotype label | slide |
| `immune_phenotype_metrics{,_v2}.csv` | 15 / 17.5 KB | 41 | Region densities, counts, areas | slide |
| `immune_compartment_distribution.csv` | 7.6 KB | 41 | Immune % per compartment | slide |
| `spatial_features.csv` | 11.8 KB | 40 | NN/DBSCAN/HDBSCAN/Moran's I | participant (derived) |
| `survival_analysis_dataset.csv` | 23 KB | 40 | Image features joined to clinical fields (Vital Status, Days to Last Followup, Days to Last Known Disease Status, Progression or Recurrence, AJCC Stage, Cause of Death, Age at Diagnosis, Gender, Tumor Grade) + os/css/pfs | **participant clinical (open-access tier)** |
| `survival_cox_results.csv`, `survival_cox_multivariate_results.csv` | 8.3 / 3.0 KB | 51 / 15 | Cox results | aggregate |
| `mutation_{univariate,stratified,interaction,burden}_cox.csv` | 1-2.8 KB | 6-23 | Cox by gene | aggregate |
| `dbscan_sensitivity.csv` | 1.6 KB | 12 | Grid results | aggregate |
| `cell_matching/matching_summary.csv` | 440 B | 1 | Matching stats, 1 slide | slide |
| `cell_matching/dg_4DFC_06CED89D-..._matched_cells.parquet` | **25 MB** | ~750k | Per-cell match table | cell; **keep out of git** (Synapse/S3) |

### 5b. metadata/ (432 KB)

| File | Size | Rows | What | Notes |
|---|---|---|---|---|
| `participant.tsv` | 16.6 KB | 40 | HTAN portal Participant export (HTA7) | See flags in 5d |
| `biospecimen.tsv` | 15.9 KB | 40 | HTAN Biospecimen export | Collection/Processing Days from Index (about 11.8k-17.4k; the index looks like birth). FFPE resections. |
| `file.tsv` | 10.7 KB | 41 | HTAN Files export, H&E Level 2 OME-TIFF | `Data Access` = "CRDC-GC/SB-CGC (open access)", Release 4.0 |
| `files_drs.csv` | 7.3 KB | 41 | `drs://nci-crdc.datacommons.io/dg.4DFC/<GUID>` | open DRS |
| `bquxjob_6358239c_19c48c9d3ac.csv` | 194 KB | 1055 | HTAN BigQuery **MolecularTest**, 96 HTA7 participants, 102 genes | Gene-level positive/negative only ("Targeted Sequencing"). AA_Change, Chromosome, Locus and Variant_Origin are all empty, so there are no variants and no germline data. |
| `bquxjob_4136c6be_19c48c7c3e6.csv` | 15.4 KB | 95 | HTAN BigQuery **Therapy** | Regimen/agents, mostly empty |
| `create_flat_files.py` | 6.1 KB | - | Script | |
| `flat/*_flat.tsv` | 18-51 KB | 40-41 | Denormalised joins | Same content as the source tables |

### 5c. Other inputs

| File | Size | What | Notes |
|---|---|---|---|
| `P37_S29_..._cellRingMask.csv` (DM root) | 813 MB, 1.62M rows | Orion mIF per-cell quantification (19 marker means, centroid, morphology) | Derived from an open HTAN Orion dataset (Lin et al. 2023). Put it on Synapse/S3, not git. |
| `classpose_outputs/s3_manifest.tsv` | 249 KB, 460 rows | S3 object manifest | Bucket URI on every row (not a secret). Strip or replace before publishing. |
| `classpose_outputs/data/**/*_cell_densities.csv` | about 290 B x 82 | region, cell_class, count, density | Fine for git |

### 5d. Controlled-access verdict: none found

| Check | Result |
|---|---|
| Raw sequencing (BAM/FASTQ/VCF), variant-level MAF | None |
| Germline / genomic coordinates | None (MolecularTest is gene-level only) |
| dbGaP-only fields (Days To Birth, Year Of Birth, Days to Death, Country, Ancestry) | Present as columns but **empty** |
| Exact calendar dates, free text, MRNs, zip codes | None |
| **Year of Death (15 pts), Year Of Diagnosis (40 pts, 2014-2016)** in `participant.tsv` and `flat/` | **FLAG (low).** Year-only values, allowed under Safe Harbor and present in the HTAN portal export. No script uses them, so consider dropping the columns before committing. |
| **Age at Diagnosis = 90 (1 participant)**, `Age Is Obfuscated` blank; the value is also in `survival_analysis_dataset.csv` | **FLAG (low).** Consistent with top-coding of ages over 89. Acceptable as released by HTAN. |

---

## 6. nextflow_reports/ secret scan

### 6a. Files (9.6 MB)

- `execution_report_{2026-01-26_16-04-11,2026-01-27_10-03-08,2026-01-27_11-56-06}.html`, about 3.3 MB each
- `execution_trace_*.txt` (same 3 timestamps), 34 KB each
- `pipeline_dag_*.svg` (same 3 timestamps), about 12.7 KB each

### 6b. Hits

The same pattern appears in all 3 HTML reports. Line numbers are for the 01-26 report; the 01-27_11-56 report is offset by about 103 lines.

| Pattern | File:line | Context (redacted) | Verdict |
|---|---|---|---|
| `accessToken` (23 per HTML) | html:811 | `s.accessToken=u.accesstoken;var h=s.map=new n.Map(` | **FP**: bundled plotly/mapbox JS |
| `eyJ`, `AKIA`/`ASIA` (case-sensitive), `tw_`, `TOWER_ACCESS_TOKEN`, `Bearer`, `/Users/ataylor`, `arn:aws`, `dkr.ecr`, emails | - | none | clean |
| `password` | html:469, 471 | `h.open(b.type,b.url,b.async,b.username,b.password)` | **FP**: jQuery |
| `secret` | html:811 | `not a secret access token (sk.*). See https://www.mapbox...` | **FP**: mapbox error string |
| 12-digit numbers (177-184 per HTML, 1 per trace) | html:1153 trace JSON; trace.txt | `"vmem":"<12DIGITS>"`, `"peak_vmem":"<12DIGITS>"`, hash fragments | **FP**: byte counts and hashes. No AWS account ID. |
| **`-params-file` ephemeral URL** | html:277 (01-26, 01-27_10-03), html:174 (01-27_11-56) | `-params-file 'https://tower.sagebionetworks.org/api/ephemeral/<REDACTED>.yaml' -with-tower ... -profile tower_gpu` | **TP (low risk)**: one-time Seqera capability URL, probably expired. **Redact before publishing.** |
| `s3://` workdir | html:251, 288; html:1153 (x164) | `s3://htan-project-tower-scratch/work/fd/<hash>` | **TP, not a secret**: leaks the scratch bucket name |
| Batch queue | html:1153 | `TowerForge-<REDACTED>-work` | **TP, not a secret**: Tower Forge compute-env ID |
| Containers | html:254, 1153 | `wave.seqera.io/wt/<id>/adamjtaylor/nf-classpose:latest`, `ghcr.io/adamjtaylor/nf-classpose{,-gen3,-vips}` | FP (public) |
| gen3 | html:1153 (x41) | `gen3-client configure --profile="htan" --cred="credentials.json" ...` | **FP**: filename only, no credential content |

**Recommendation:**
- Commit only the 3 `execution_trace_*.txt` files. These are all that `plot_classpose_runtimes.py` reads, and they contain no bucket, token or queue strings.
- If the HTML reports are published, redact the `-params-file` URL. Optionally also scrub `htan-project-tower-scratch` and `TowerForge-*`.

---

## 7. nf-classpose and classpose clones

### 7a. nf-classpose (`~/Documents/learn-classpose/nf-classpose`)

| Item | Value |
|---|---|
| Remote | `https://github.com/adamjtaylor/nf-classpose.git` |
| Branch | `claude/issue-3-20260125-2131` |
| HEAD | `dcb045250a6df043ba0b4eb5912d8b016855aeeb` (2026-01-26 15:56 UTC) |
| Status | Clean working tree. Level with `origin/claude/issue-3-20260125-2131`. |
| vs main | **9 commits ahead, 0 behind, not merged.** main = origin/main = `eab0ffb` |
| Run provenance | Run 1 (`maniac_panini`, 2026-01-26 16:04) used `-r claude/issue-3-20260125-2131 -latest`, which resolved to dcb0452. Runs 2 and 3 pinned `-r dcb045250a6d` with `-resume`. All runs used `-profile tower_gpu`. |

`git log --oneline -5`:
```
dcb0452 Add vips_save parameter and update publishDir in VIPS_CONVERT process
d54f3d8 Update output file naming to include model name and add stub credentials
d5bdfe9 Fix model parsing and add stub test fixtures
8a549e4 Add stub blocks to all processes for dry-run testing
b084e46 Add optional id column support and refactor model handling
```

`git log main..HEAD` has 9 commits: the 5 above plus `b4cb198` (merge of PR #12), `ae41e91`, `e24defd`, and `13ff2aa` ("Add consep model and implement model matrix selection").

`git branch -a`:
- Local: `claude/issue-3-20260125-2131` (current), `claude/issue-5-20260121-2104`, `feature/drs-uri-support`, `fix/vips-resolution`, `initial-dev` (3 behind upstream), `main`.
- Remote only: `add-claude-github-actions-1768991628324`, `claude/issue-3-20260125-2138`, `claude/issue-3-20260125-2155`, `claude/issue-8-20260124-1745`.

**Params and samplesheets for the htan_hms_crc_he_41_dual_model run**

| Path | Tracked | Contents |
|---|---|---|
| `nf-classpose/test/htan_hms_crc_he_41.csv` | yes (also on main) | **The run's samplesheet.** `slide_path` column with 41 `drs://nci-crdc.datacommons.io/dg.4DFC/<GUID>` rows. No absolute paths, S3, Synapse IDs, account IDs or profiles. |
| Params file | n/a | **Not saved**: it was an ephemeral Seqera upload. Rebuild it from `nextflow.config` at dcb0452 plus `-profile tower_gpu`. outdir is presumably `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model`, inferred from the fetch scripts. |
| `example_run/fabulous_faggin-params.json` | no (gitignored) | Earlier DRS test params; useful as a template. Includes `gen3_credentials: s3://htan-project-tower-bucket/nf-classpose/credentials.json` and `gen3_profile: htan`. |
| `test/htan_samplesheet.csv`, `test/tower_samplesheet.csv` | yes | `s3://htan-dcc-hms/ORION202305_imaging_level_2/...` and `s3://htan-project-tower-bucket/histoqc/examples/*.svs` |
| `nextflow.config:113` | yes | `tower_test` profile outdir `s3://htan-project-tower-bucket/classpose/results`; `gen3_profile = 'htan'` |

- No nf-classpose file mentions `dual_model`.
- DM has no samplesheet or params files; `data/sample_manifest.txt` is the nearest equivalent.

**Secret scan of the 27 tracked files:**
- None of `/Users/`, `AKIA`, 12-digit account IDs or `accessToken`.
- `tower` appears only as profile names.
- `test/stub_credentials.json` holds a dummy key (FP).

**Untracked files:**
- **TRUE POSITIVE:** `nf-classpose/example_run/logs/nf-3stNqO29UUEC1M.log` lines 131 and 160 contain `towerAccessToken:eyJ<REDACTED>` and `towerRefreshToken:eyJ<REDACTED>`. These are real Seqera tokens from January 2026. The file is gitignored and untracked. **Never copy `example_run/`, and consider revoking the token.**
- Root `.nextflow.log*` files (ignored) contain 19 `/Users/` paths.

### 7b. classpose clone (`~/Documents/learn-classpose/classpose`)

| Item | Value |
|---|---|
| Remote | `https://github.com/sohmandal/classpose.git` (upstream) |
| Branch / HEAD | `main`, level with origin; `617508c07b5c2b40be56c262c2e2a453f9811b13` (merge of PR #6, bump-cp) |
| Modified tracked files | none |
| Untracked | `.DS_Store`, `cell_mosaic.py`, `models/`, `paper_experiments/.DS_Store`, `spatialdata_plots.py`, `visualize_cells.py` |

- The 3 untracked scripts are experiments on the public CMU-1-Small-Region slide. They hardcode `~/Documents/GitHub/nf-wsi-vignette/...` and `~/Documents/learn-classpose/output/`. Not needed for the demo.
- `models/` (50 MB, untracked) holds `grandqc_artefact.pt` (25.4 MB, md5 `ece5ac99...`) and `grandqc_tissue.pt` (26.6 MB, md5 `25ccb259...`). These look like renamed GrandQC weights; they were not checksum-verified against Zenodo.
- **All weights can be downloaded upstream, so do not ship `models/`:**
  - Classpose weights: HuggingFace `classpose/classpose` (`src/classpose/model_configs.py:22`). Files are conic/consep/glysac/monusac/nucls/puma `.pt`, auto-downloaded to `$CLASSPOSE_MODEL_DIR` (default `~/.classpose_models`). `nf-classpose/docker/Dockerfile:33-35` builds conic.pt and consep.pt into the image.
  - GrandQC tissue: `https://zenodo.org/records/14507273/files/Tissue_Detection_MPP10.pth`
  - GrandQC artefact: `https://zenodo.org/records/14041538/files/GrandQC_MPP1.pth`. This is the classpose default; the nf-classpose Dockerfile uses `GrandQC_MPP15.pth` from the same record.
  - timm `tf_efficientnet_b0_aa-827b6e33.pth` from GitHub releases (Dockerfile:49).
- Licences: Classpose and GrandQC are both non-commercial. Cite GrandQC in the demo README, as the classpose README requests.

---

## 8. Older run, pannuke, output/

### 8a. `htan_hms_crc_he_41/` (5.2 GB; single slide HTA7_954 = `dg_4DFC_6C434FEA-...`)

Every script there hardcodes `/Users/ataylor`; the count per script is in brackets.

| File | Size | Purpose | Keep? |
|---|---|---|---|
| `TLS_ANALYSIS_METHODS.md` | 10.3 KB | TLS methods: rationale, clustering, metrics, literature-based staging, composite score, limitations, refs | **Yes** (methods doc) |
| `tls_revised_analysis.py` | 22 KB | Revised literature-based TLS staging, matching the md [5] | Yes, if the TLS story is kept |
| `tls_radius_analysis.py` | 16 KB | 50/100/200 um radius comparison [3] | Optional |
| `find_tls.py`, `tls_graph.py`, `tls_heatmap.py` | 5-7 KB | Early TLS detection, graph, heatmap [2-3] | Superseded |
| `tls_maturity.py`, `tls_maturity_logscale.py` | 5-8 KB | Plasma/lymphocyte maturity, the interpretation the md corrects [2] | No |
| `hierarchical_classification_figure.py` | 23 KB | CoNIC+ConSeP hierarchical classification figure [1] | Maybe (explains the dual-model choice) |
| `morphometric_analysis.py` | 24 KB | Cell shape/size from contours [10] | Optional |
| `analyze_cell_densities.py`, `cell_density_visualizations.py`, `plot_cell_types.py`, `plot_cell_types_crops.py` | 4-10 KB | Density bars, clustermap/PCA, spatial plots and crops | Superseded by DM |

- There are about 40 PNG/PDF outputs, plus `cell_morphometrics.parquet` (82 MB) and `hierarchical_cell_classification.parquet` (14 MB).
- There is a stray file literally named `" .png"`.
- 41 `dg_4DFC_*/` directories hold CoNIC-only density CSVs.

**What DM reads from here:**

| File | Size | Read by | Also in DM `sample_data/`? |
|---|---|---|---|
| `dg_4DFC_6C434FEA-....tif` | 2.0 GB | 2a, 2b, 2c, 2d, 3a, 3b, select_representative_samples, spatial_utils | **No**: upload to Synapse/S3 |
| `dual_model_output/consep/dg_4DFC_6C434FEA-..._consep_spatialdata.zarr` | 357 MB | 2b-2e, 3a, 3b, select_representative_samples, spatial_utils | Yes (same size; `points.parquet` has the same md5) |

### 8b. `pannuke/` (4.5 GB: data 3.4 GB, checkpoints 1.1 GB; own git repo with no commits; unrelated to the CRC demo)

| File | Size | Purpose |
|---|---|---|
| `claude.md` | 2.3 KB | Notes on training Classpose on PanNuke (HF `RationAI/PanNuke`) |
| `convert_pannuke.py` | 4.0 KB | PanNuke to Classpose label format |
| `run_training.py` | 16 KB | Training script (from classpose `paper_experiments`) |
| `evaluate.py` | 9.9 KB | Per-class PQ on fold 3 vs the PUMA baseline |
| `vertex_job.py` | 2.8 KB | Vertex AI submit; hardcodes `PROJECT="htan-dcc"`, `gs://htan-dcc-pannuke` and the `us-central1-docker.pkg.dev/htan-dcc/pannuke/trainer` image |
| `trainer/task.py` | 2.0 KB | Vertex entrypoint with spot-preemption resume |
| `Dockerfile` | 1.6 KB | PyTorch CUDA + classpose |
| `explore.ipynb` | 6.0 KB | Exploration |
| `pyproject.toml`, `requirements.txt`, `setup.py`, `uv.lock` | small | Packaging |

### 8c. `output/` (9.2 MB)

- Contents: ClassPose outputs on the public CMU-1-Small-Region slide.
  - geojsons: centroids, contours, tissue, artefact
  - a densities CSV and a SpatialData zarr
  - 4 PNGs: cell_mosaic, cell_visualization, spatialdata_analysis, spatialdata_density_by_type
- No HTAN data.
- Not needed, though it could serve as a tiny smoke-test fixture.

### 8d. Other top-level items

| Item | Notes |
|---|---|
| `htan1_examples/` (1.3 GB) | 6 HTAN1 H&E OME-TIFFs + a Synapse manifest CSV (Synapse IDs, Synapse user IDs, participant/biospecimen IDs; open tier). Not needed. |
| `htan_study_plan.md` | Plan for a pan-cancer HTAN ClassPose phenotyping study |
| `first_plan.txt` (100 KB) | Raw Claude Code session transcript. **Do not publish.** |
