# Methods

## Dataset

Forty patients with colorectal adenocarcinoma were drawn from the HTAN HMS atlas (Human Tumor Atlas Network, Harvard Medical School). One patient (HTA7_919) contributed two near-serial H&E sections, yielding 41 whole-slide images (WSIs) in total. H&E images were acquired as part of the Orion multiplexed imaging platform and are available at HTAN Level 2 as registered OME-TIFF files. Clinical and demographic metadata, biospecimen records, and targeted sequencing results for 11 CRC driver genes were obtained from the HTAN data portal and are provided in `metadata/`. Survival endpoints (overall survival, cancer-specific survival, progression-free survival) were derived from harmonised HTAN participant-level data.

## Cell Segmentation and Classification

### nf-classpose pipeline

Cell segmentation and classification were performed using ClassPose (v1.0.0) via **nf-classpose** (github.com/ncihtan/nf-classpose), a Nextflow wrapper that handles WSI ingestion, tiling, inference, and output aggregation. The pipeline was executed on AWS Batch via Seqera Platform (run name: `maniac_panini_4`; resumed from session `5293201e-0771-47a2-bf8a-c3b06a1767fe`) on 27 January 2026.

**Input manifest:** `https://raw.githubusercontent.com/ncihtan/nf-classpose/refs/heads/main/test/htan_hms_crc_he_41.csv`

**Key pipeline parameters:**

| Parameter | Value |
|---|---|
| Models | `conic,consep` |
| Tile size | 1024 px |
| Tile overlap | 64 px |
| Batch size | 8 |
| Test-time augmentation | disabled |
| bfloat16 | disabled |
| Tissue detection model | `Tissue_Detection_MPP10.pth` |
| Artefact detection model | `GrandQC_MPP15.pth` |
| Artefact filtering | disabled |
| Output formats | CSV + SpatialData zarr |

**Compute resources (AWS Batch):**
- GPU processes: 8 CPUs, 128 GB RAM, 1 GPU accelerator, 8 h wall time
- Default processes: 4 CPUs × attempt, 16 GB × attempt, 4 h × attempt (max 2 retries)
- Container: `ghcr.io/ncihtan/nf-classpose:latest`
- Region: `us-east-1`

Outputs were written to `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/`. Pipeline execution report, trace, and DAG are stored at the same location under `pipeline_info/`.

### Cell classification models

Both the **ConSeP** and **CoNIC** pre-trained models were applied to each slide:

- **ConSeP** (primary model): trained on the Colorectal Adenocarcinoma dataset. Classifies cells into six categories: malignant epithelial, healthy epithelial, stroma, muscle, inflammatory, other. Total cells classified: 50,191,395.
- **CoNIC**: trained on the CoNSeP + Lizard datasets. Classifies cells into six categories: epithelial, connective, lymphocyte, plasma cell, neutrophil, eosinophil. Total cells classified: 39,283,683.

The difference in total cell count reflects model-level differences in detection sensitivity and what constitutes a classifiable cell, not a data quality issue. Cross-model concordance was assessed by mapping both schemas to three broad biological categories (immune, epithelial, stromal) and computing Pearson correlation across the 41 images (immune: r = 0.73; epithelial: r = 0.90; stromal: r = 0.82; all p < 10⁻⁷).

Per-sample cell density outputs (cells per tissue region, in counts and cells/mm²) were stored as CSV files under `classpose_outputs/data/HTA7_*/`.

---

## Downstream Analysis

All downstream analysis was performed in Python 3.13 using the packages listed in `pyproject.toml`. Key package versions (as pinned in `uv.lock`):

| Package | Version |
|---|---|
| pandas | 3.0.0 |
| numpy | 2.3.5 |
| scipy | 1.17.0 |
| scikit-learn | 1.8.0 |
| matplotlib | 3.10.8 |
| seaborn | 0.13.2 |
| spatialdata | 0.6.1 |
| squidpy | 1.7.0 |
| fast-hdbscan | 0.2.2 |
| lifelines | 0.30.0 |

Scripts are in `scripts/`; analysis was run from the project root using `uv run`.

### Cell density loading and composition metrics

`scripts/load_cell_densities.py` — shared utility module. Loads per-sample ConSeP and CoNIC cell density CSVs from `classpose_outputs/`, pivots to wide format (one row per sample), and computes:
- Percentage of each cell class relative to total cells per slide
- Centred log-ratio (CLR) transformed proportions (using geometric mean across classes, with pseudocount replacement for zeros) for clustering and ordination
- Neutrophil-to-lymphocyte ratio (NLR) from CoNIC counts

### Immune phenotype classification

`scripts/classify_immune_phenotypes_v2.py` — KNN-graph-based region definition with absolute density thresholds derived from the immunoscore literature.

**Region assignment (per cell):** A K-nearest-neighbour graph (k = 30) was constructed from cell centroid coordinates using Squidpy. Each cell was assigned a microenvironment region based on the composition of its 30 nearest neighbours:
- **Tumour core:** ≥ 60% malignant neighbours
- **Tumour-stroma interface (TSI):** 15–60% malignant neighbours
- **Stroma:** ≥ 40% stroma neighbours, < 10% malignant

**Immune density calculation:** For each cell, a territory radius was estimated as half the distance to its nearest neighbour (clipped to 1–100 px). Territory area was computed as πr² and converted to mm² using a pixel size of 0.25 µm/px. Immune cell density (cells/mm²) was summed across the slide.

**Phenotype thresholds:**

| Parameter | Value |
|---|---|
| Hot threshold (immune density) | 6,000 cells/mm² |
| Desert threshold (immune density) | 2,500 cells/mm² |
| Tumour-associated infiltration threshold | 20% of immune cells in tumour core + TSI |
| Pixel size | 0.25 µm/px |

**Classification rules:**

| Immune density | Tumour-associated % | Phenotype |
|---|---|---|
| ≥ 6,000 cells/mm² | ≥ 20% | Hot |
| ≥ 6,000 cells/mm² | < 20% | Altered-excluded |
| 2,500–6,000 cells/mm² | ≥ 20% | Altered-immunosuppressed |
| 2,500–6,000 cells/mm² | < 20% | Altered-excluded |
| < 2,500 cells/mm² | any | Immune-desert |

Outputs: `data/immune_phenotype_metrics_v2.csv`, `data/immune_phenotype_classification_v2.csv`.

### Spatial feature computation

`scripts/compute_spatial_features.py` — computes truly spatial (coordinate-based) features not derivable from density alone.

**Nearest-neighbour distance:** Median distance (px) from each inflammatory cell to its nearest inflammatory neighbour, computed using a KD-tree on cell centroids.

**DBSCAN clustering:**

| Parameter | Value |
|---|---|
| eps | 100 px (~32.5 µm, ~2–3 cell diameters) |
| min\_samples | 30 cells |
| n\_jobs | −1 (all cores) |

Cluster fraction = proportion of inflammatory cells assigned to a cluster (non-noise). Tumour cluster fraction = proportion of clustered inflammatory cells whose centroid falls within a tumour core bin.

**HDBSCAN clustering:**

| Parameter | Value |
|---|---|
| min\_cluster\_size | 30 cells |

**Moran's I:** Spatial autocorrelation of tumour-immune co-localisation, computed on a masked grid of immune cell density within tumour-proximal bins.

Outputs: `data/spatial_features.csv`.

### Survival analysis

`scripts/survival_analysis.py` — univariate Cox proportional-hazards regression.

**Endpoints:**

| Endpoint | Time variable | Event variable |
|---|---|---|
| Overall survival (OS) | `os_time` | `os_event` |
| Cancer-specific survival (CSS) | `css_time` | `css_event` |
| Progression-free survival (PFS) | `pfs_time` | `pfs_event` |

**Features tested (18 total):**
- Clinical: age at diagnosis, sex (male), stage (binary I–II vs III–IV), tumour grade (high)
- ConSeP spatial: inflammatory %, infiltration ratio, tumour-associated %, overall inflammatory density, tumour-stroma ratio
- CoNIC-derived: NLR, eosinophil density
- Spatial: median NN distance, DBSCAN cluster fraction, DBSCAN tumour cluster fraction, HDBSCAN cluster fraction, HDBSCAN tumour cluster fraction, tumour-immune Moran's I

Continuous features were standardised (z-score, per feature across the cohort) prior to fitting so that hazard ratios are expressed per standard deviation. Models were fit using the `lifelines` CoxPHFitter with default L2 penaliser. Patients with missing survival data were excluded per endpoint (n = 39 evaluable for OS and PFS; n varies for CSS due to cause-of-death availability). Kaplan-Meier curves were stratified at the cohort median. Outputs: `data/survival_analysis_dataset.csv`.

### Mutation interaction analysis

`scripts/mutation_interaction.py` — formal interaction Cox models testing whether mutation status modifies the prognostic effect of spatial immune features.

For each gene–feature–endpoint combination, an interaction term (mutation status × standardised feature) was added to the Cox model. The interaction p-value was extracted from the Wald test on the interaction coefficient. Stratified models were also fit separately within mutation-positive and mutation-negative subgroups to report stratum-specific hazard ratios. Sensitivity analyses excluded BRAF-mutant patients (likely MSI-H) and applied permutation testing (1,000 permutations) to assess specificity of the KRAS interaction. Outputs: `data/mutation_interaction_cox.csv`, `data/mutation_stratified_cox.csv`, `data/mutation_univariate_cox.csv`.
