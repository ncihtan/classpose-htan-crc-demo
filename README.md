# ClassPose on HTAN HMS colorectal cancer H&E

This repo is demo work from the HTAN Data Coordinating Center. We ran [ClassPose](https://github.com/sohmandal/classpose) cell segmentation and classification on 41 open-access H&E whole-slide images from 40 colorectal cancer patients in the HTAN HMS atlas. We then built spatial features and immune phenotypes from the cell calls, and tested them against survival and driver mutations. It is exploratory work, not a paper. The cohort is small (n = 39 to 40), and the p-values are not corrected for the number of features tested.

Start with [`analysis/classpose_meeting_deck.md`](analysis/classpose_meeting_deck.md) for the story, then [`docs/methods.md`](docs/methods.md).

## Headline results

| Result | Value | Details |
|---|---|---|
| Cells classified | CoNSeP 50.2M, CoNIC 39.3M | 41 slides, both models |
| Cross-model agreement (per-slide proportions) | immune r = 0.73, epithelial r = 0.90, stromal r = 0.82 | `figure1b_validation.png` |
| Near-serial sections (HTA7_919) | inflammatory 19.2% vs 18.1%, stroma 34.9% vs 33.8% | `figure1d_serial_sections.png` |
| Immune phenotypes | hot 22%, altered-excluded 44%, altered-immunosuppressed 22%, immune-desert 12% | `classify_immune_phenotypes_v2.py` |
| Inflammatory % vs PFS | HR 0.44 per SD, p = 0.003 (univariate Cox) | `docs/survival_analysis.md` |
| KRAS x tumour-proximal immune clustering (OS) | interaction p = 0.032; KRAS-mutant high-clustering HR 0.10 | `docs/mutation_analysis.md` |
| Check against Orion mIF (same section, one slide) | 1.2M cells matched; CD45 up in Inflammatory, Pan-CK up in epithelial classes | `match_mif.py`, `match_cells.py` |
| GPU time | median 1.5 h (CoNIC) and 1.7 h (CoNSeP) per slide, about 130 h in total | `classpose_runtimes.png` |

## Data and compute

| Item | Where | Size |
|---|---|---|
| Input H&E images (HTAN Level 2 OME-TIFF, open access) | CRDC, by DRS URI: [`analysis/metadata/files_drs.csv`](analysis/metadata/files_drs.csv) | 41 files |
| Clinical, biospecimen and mutation tables (HTAN open-access tier) | [`analysis/metadata/`](analysis/metadata/) | under 1 MB |
| nf-classpose outputs (tifs, geojsons, csvs) | `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/`, listed in [`pipeline/s3_manifest.tsv`](pipeline/s3_manifest.tsv) | 257 GB, including 82 zarrs |
| SpatialData zarrs and other files the analysis needs | Synapse [syn77984172](https://www.synapse.org/Synapse:syn77984172), listed in [`data/synapse_manifest.tsv`](data/synapse_manifest.tsv) | about 21 GB |
| Per-slide cell density CSVs and derived tables | [`analysis/classpose_outputs/`](analysis/classpose_outputs/), [`analysis/data/`](analysis/data/) | small |
| Compute | nf-classpose on Seqera Platform, AWS Batch GPU spot, January 2026 | see [`pipeline/`](pipeline/) |

Neither the S3 bucket nor the Synapse project is public yet. Ask the HTAN DCC for access.

## Layout

| Folder | What it holds |
|---|---|
| `pipeline/` | The nf-classpose commit, samplesheet, reconstructed params, S3 output manifest and Nextflow reports |
| `analysis/` | Analysis scripts, `pyproject.toml` and `uv.lock`, metadata, derived tables and figures. It mirrors the original working folder, so relative paths still work. |
| `analysis/scripts/config.py` | The one place for paths. Set `SAMPLE_DATA_DIR`, `MIF_CSV`, `CLASSPOSE_S3_BASE` or `AWS_PROFILE` to override. |
| `data/` | Synapse manifest, plus the upload and download scripts |
| `docs/` | Methods, survival, mutation and TLS notes, the figure map, the session history and the salvage records |
| `archive/single-slide-conic/` | The earlier CoNIC-only run on one slide (HTA7_954), including the TLS analysis. Kept for reference. |
| `pannuke/` | A small ClassPose training demo on the public PanNuke dataset. It did not produce a useful model. |
| `patches/classpose-local/` | Three helper scripts written against a local ClassPose clone |

## Rerun

1. Install [uv](https://docs.astral.sh/uv/), then run `cd analysis && uv sync`. Python 3.13 is required.
2. Get the zarrs and the other inputs from Synapse. You need Synapse access and a token in `~/.synapseConfig`.

   ```bash
   uv run --no-project --with synapseclient python data/fetch_from_synapse.py
   ```

   This puts the files in `analysis/sample_data/`, which git ignores. To keep them somewhere else, set `SAMPLE_DATA_DIR`.
3. Run a script from `analysis/`, for example `uv run python scripts/survival_analysis.py`. [`docs/figure_map.tsv`](docs/figure_map.tsv) lists the script and inputs for each figure.
4. To rerun ClassPose itself, see [`pipeline/README.md`](pipeline/README.md).

Six entries in the figure map have no matching script and are marked `unknown`. Three full-resolution figures are on Synapse only. Downsized `*_preview.jpg` copies are in `analysis/figures/`.

## Credits

- ClassPose: Mandal S, de Almeida JG, Papanikolaou N, Graham TA. Classpose: foundation model-driven whole slide image-scale cell phenotyping in H&E. bioRxiv 2025. [doi:10.64898/2025.12.18.695211](https://doi.org/10.64898/2025.12.18.695211). ClassPose is released under CC BY-NC 4.0, and its tissue and artefact models come from [GrandQC](https://github.com/cpath-ukk/grandqc), which also has a non-commercial licence.
- HTAN HMS colorectal cancer atlas: Lin JR, Wang S, Coy S, et al. Multiplexed 3D atlas of state transitions and immune interaction in colorectal cancer. Cell 2023;186(2):363-381. [doi:10.1016/j.cell.2022.12.028](https://doi.org/10.1016/j.cell.2022.12.028)
- Orion same-section H&E and immunofluorescence: Lin JR, Chen YA, Campton D, et al. High-plex immunofluorescence imaging and traditional histology of the same tissue section for discovering image-based biomarkers. Nat Cancer 2023;4(7):1036-1052. [doi:10.1038/s43018-023-00576-1](https://doi.org/10.1038/s43018-023-00576-1)
- The Human Tumor Atlas Network (HTAN), funded by the NCI Cancer Moonshot. Data from the [HTAN Data Portal](https://data.humantumoratlas.org).
- PanNuke (Gamper et al.), used through Hugging Face `RationAI/PanNuke`.

## AI assistance

Most of the code, docs and figures here were written with Claude Code (Anthropic), working with Adam Taylor, between January and March 2026. In October 2026, Claude Code also copied the work into this repo and wrote the README and the salvage docs. A person directed every analysis and reviewed the results, but the code has not been independently reviewed. Check any result before you rely on it.

## Licence

The code in this repo is under the MIT licence (see `LICENSE`). ClassPose, its models and the input data keep their own licences and terms.
