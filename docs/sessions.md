# Past Claude Code sessions for learn-classpose

This file summarises the Claude Code sessions that built nf-classpose and analysed ClassPose output for the 41 HTAN HMS CRC H&E slides, written on 2026-10-15 while the project was being salvaged. None of the past transcript `.jsonl` files exist any more; Claude Code had already deleted them. The only `.jsonl` left in the root folder is the current salvage session, which is excluded here. The summary therefore comes from five indirect sources: the `sessions-index.json` files (first prompt, auto summary, message count, git branch), the user prompts in `~/.claude/history.jsonl`, two `memory/` folders, the terminal export `first_plan.txt`, and the headings of `htan_study_plan.md`. A few points were checked read-only against the working tree, git and GitHub, and are marked as such. Bash commands, errors and tool output were not recorded in these sources, so anything here about errors and fixes is pieced together from prompts and commit messages.

## Sessions

Sizes are message counts from `sessions-index.json` where an index exists, and prompt counts from `history.jsonl` otherwise. No transcript files survive.

| Folder | Session | Dates | Size | First prompt (paraphrased) | Topic |
|---|---|---|---|---|---|
| root | 3ad6b241 | 2026-01-06 | 16 prompts; `first_plan.txt` export (1,493 lines) | Plan a ClassPose demo using a uv venv and an example SVS | Local ClassPose demo, GrandQC models, spatialdata bug, plots, `htan_study_plan.md` |
| root | afd0f742 | 2026-10-15 | current | Salvage plan; large data to a new Synapse project | Excluded (current session) |
| nf-classpose | e496073d | 2026-01-21 | 2 msgs | (none, exit) | Empty |
| nf-classpose | 7428b7e4 | 2026-01-21 | 9 msgs | Write a light Nextflow wrapper for ClassPose with a samplesheet | Wrapper design, Docker image with models baked in |
| nf-classpose | 0ed92453 | 2026-01-21 | 47 msgs, branch initial-dev | Implement the nf-classpose plan | First pipeline, PR, samplesheet params, profiles |
| nf-classpose | acfdd341 | 2026-01-21 | 79 msgs, initial-dev | Allow workflow_dispatch of the Docker build on a branch | Container builds, Seqera/AWS Batch runtime fixes |
| nf-classpose | 214e577d | 2026-01-21 | 12 msgs | Add a module to fetch DRS URIs with gen3-client | DRS download design |
| nf-classpose | c41a1f65 | 2026-01-21 to 01-22 | 25 msgs, claude/issue-5-20260121-2104 | Implement the DRS via Gen3 plan | DRS PR, gen3 container, VIPS branch |
| nf-classpose | d05a49f2 | 2026-01-23 | 31 msgs, issue-5 | Test the VIPS container converter on an HTAN OME-TIFF | OME-TIFF conversion testing |
| nf-classpose | 7bba86a8 | 2026-01-23 | 3 msgs, main | (opened tower samplesheet in IDE) | Tower samplesheet |
| nf-classpose | ad9edd9c | 2026-01-25 | 16 msgs, main | (opened example_run execution report in IDE) | Example run review |
| nf-classpose | 367d8f36 | 2026-01-26 | 35 msgs, claude/issue-3-20260125-2131 | (opened test/htan_hms_crc_he_41.csv in IDE) | Multi-model run, id column, stub tests |
| nf-classpose | b046542f | 2026-01-26 | 5 msgs, issue-3 | Add a param to save the VIPS-converted output | `vips_save` param |
| htan-hms-crc-he-41 | 3e1454ff | 2026-01-25 | 6 msgs | Stacked bar plots of cell type fractions and densities | Restarted straight away |
| htan-hms-crc-he-41 | 39af9d57 | 2026-01-25 | 24 msgs | Same stacked bar prompt | Composition, immune ranking, TLS idea |
| htan-hms-crc-he-41 | 139a94cb | 2026-01-25 to 01-26 | 31 msgs | Plot cell types from one spatialdata zarr in a venv | Spatial plots, TLS neighbourhood, graph, maturity |
| htan-hms-crc-he-41 | ad491579 | 2026-01-26 | 21 msgs | Implement the revised TLS analysis plan | TLS staging, radius, crops, morphometrics, cell totals for email |
| dual-mode | 9bf7f2a3 | 2026-01-27 | 13 prompts | Work out keys across metadata files and flatten them | Metadata join, S3 manifest, per-model bars, Figure 1 idea |
| dual-mode | 1fa4cbf6 | 2026-01-27 | 24 prompts | Implement the Figure 1 plan | Figure 1 panels, ConSeP as primary, 41 files |
| dual-mode | 64dc8ae9 | 2026-01-27 | 14 prompts | Implement the Figures 2 and 3 plan | Hot and cold examples, zooms, phenotype plan |
| dual-mode | ecaaa35b | 2026-01-27 to 01-28 | 22 prompts | Implement the tumour immune phenotype classification plan | Phenotype v1 then v2 (KNN, squidpy) |
| dual-mode | 210059d3 | 2026-02-09 to 02-10 | 18 prompts | Mosaic of all samples with the malignant component | Figure 3a overview, misclassification |
| dual-mode | 602ec9a2 | 2026-02-10 | 17 prompts | Move the cellRingMask CSV from Downloads into the project | Orion mIF matching to ClassPose |
| dual-mode | 2e72106f | 2026-02-10 | 2 prompts | Cell-to-cell matching of CoNIC and ConSeP segmentations | Plan only |
| dual-mode | 165e60ad | 2026-02-10 | 21 prompts | Spatialdata zarrs are on S3, use the tower profile | CoNIC vs ConSeP matching, survival idea |
| dual-mode | d9273f0e | 2026-02-10 | 6 prompts | Days to last follow-up looks risky | Survival endpoints and write-up |
| dual-mode | 579fd557 | 2026-02-10 | 7 prompts | Define immune cluster fraction in one sentence | DBSCAN vs Moran's I |
| dual-mode | 6ef4d677 | 2026-02-10 | 20 prompts | Something broke; what happened | Tumour-associated DBSCAN, HDBSCAN, covariates, mutation idea |
| dual-mode | f8deaa88 | 2026-02-10 to 02-11 | 15 prompts | Why is that missing from the mutation forest plot | Mutation interaction, literature check |
| dual-mode | e30d0e93 | 2026-03-12 | 16 prompts | Review work and write a "picking ClassPose back up" report | Meeting report, markdown deck, methods, Google Slides |
| dual-mode | 00299338 | 2026-03-12 | 30 prompts | Find in Claude logs where the mIF matching was done | Recovered `match_mif.py`, runtime figure, pan-cancer scope, PanNuke note |
| pannuke | 99c0ca04 | 2026-03-12 | 3 prompts | Read and refine the CLAUDE.md | Training plan, time estimates |
| pannuke | eb2f3e1e | 2026-03-12 to 03-14 | 59 prompts | What's next | 2% MPS demo, Vertex AI training, container builds |
| pannuke | 34155e93 | 2026-03-13 | 22 prompts | How much disk OrbStack images and builds use | Disk clean-up (Docker, uv, npm, Downloads) |

There are no transcripts or index files for the dual-mode and pannuke folders. The nf-classpose sessions after c41a1f65 (d05a49f2, 7bba86a8, ad9edd9c, 367d8f36, b046542f) have no prompts in `history.jsonl` and were probably run from the IDE extension. They are known only from the index.

## Date range

The past sessions run from 2026-01-06 to 2026-03-14. After that the next activity is the salvage session on 2026-10-15.

## Timeline

- 2026-01-06: ClassPose is run locally on MPS. The full-size SVS was too slow and was killed at about 11%, so the demo switched to `CMU-1-Small-Region.svs`. GrandQC models went into `classpose/models/`. spatialdata was patched by hand. `visualize_cells.py`, `cell_mosaic.py`, `spatialdata_plots.py` and `htan_study_plan.md` were written.
- 2026-01-21: nf-classpose was built (PR #1 initial-dev, PR #2 Claude GitHub actions). Containers were built on GitHub Actions and the pipeline was tested on Seqera Platform with AWS Batch GPU (spot). There were many runtime fixes. DRS support was planned and PR #4 opened. PR #6 (VIPS OME-TIFF conversion) was started on claude/issue-5-20260121-2104.
- 2026-01-23: VIPS converter tested on `htan1_examples/1027015.ome.tiff`. PR #6 merged.
- 2026-01-24: PR #7 (VIPS resolution from PhysicalSizeX) and PR #9 (converted file naming) merged.
- 2026-01-25: PR #4 (DRS) merged. Analysis of CoNIC-only outputs starts in `htan_hms_crc_he_41`. Branches claude/issue-3-20260125-2131, -2138 and -2155 created for issue #3 (bundle more models).
- 2026-01-26: TLS and morphometric analyses. A draft email to José and Soham (ClassPose developers) reports 41 HTAN Phase 1 CRC images run on Seqera with spot AWS Batch GPU. On issue-3: PR #12 (2155 into 2131) merged, id column, stub tests, model name in filenames, `vips_save`. First dual-model execution report at 16:04.
- 2026-01-27: PR #11 (issue-3-2131 into main) merged at 09:54 and issue #3 closed. Dual-model execution reports at 10:03 and 11:56. `htan_hms_crc_he_41_dual_mode` started: metadata flattening, S3 manifest, Figure 1, Figures 2 and 3, immune phenotypes v1 and v2.
- 2026-02-09 to 02-11: spatial overview mosaic, misclassification review, Orion mIF matching, CoNIC vs ConSeP matching, survival analysis, DBSCAN/HDBSCAN/Moran's I features, mutation interaction, literature comparison.
- 2026-03-12: preparation for the meeting with the ClassPose developers (report, markdown deck, `methods.md`, runtime plot, Google Slides through the Workspace CLI, image.sc forum reply, pan-cancer scoping). PanNuke 2% demo trained on MPS. Vertex AI L4 training job launched.
- 2026-03-13 to 03-14: Vertex jobs failed. Docker builds bloated OrbStack to about 164 GB, so container building was paused and the disk cleaned up.

## Findings

### Pipeline and compute (nf-classpose, AWS Batch, Seqera)

- Repo: `https://github.com/adamjtaylor/nf-classpose`. It takes a CSV samplesheet (simplified to the image path plus optional id) and has `tower`, `gpu` and `tower_test` profiles. Accelerator type is not pinned because it is not known on Seqera. Seqera needed a `main.nf` on main before it could launch from a branch, so an empty one was added.
- Containers: `ghcr.io/adamjtaylor/nf-classpose:latest` (also tagged `:buildcache`), `ghcr.io/adamjtaylor/nf-classpose-gen3:latest`, and an `nf-classpose-vips` image. The first request (`ghcr.io/sohmandal/...`) was changed to Adam's namespace.
- The GitHub runner (14 GB disk) ran out of space building with all models. The workaround was CoNIC only plus issue #3. ConSeP and a model matrix were added later through PR #11.
- Nextflow could not find `ps`, so procps was added to the Dockerfiles.
- Adam asked why EfficientNet weights (timm `tf_efficientnet_b0_aa`) were downloaded inside the process. The aim was to bake every model into the image. The final fix is not recorded.
- A test run with `output_type "csv spatialdata"` published only GeoJSON. The work dir was inspected. The resolution is not recorded.
- The retry strategy did not scale resources for GPU processes. Raised; the fix is not recorded.
- OOM despite large memory was traced to artefact detection and thumbnail generation. Artefact detection is now off by default behind a param (commit "Disable artefact detection by default").
- ClassPose does not read OME-TIFF, so a VIPS conversion module was added (PRs #6, #7, #9, plus `vips_save` to publish the converted file). Adam objected to deciding in the script whether to convert, and the logic moved into Nextflow branching (commit "Refactor vips conversion to use Nextflow branching").
- DRS URIs (`drs://nci-crdc.datacommons.io/dg.4DFC/<guid>`) are fetched with gen3-client in a separate container. The pipeline fails fast if credentials are missing; the credentials path is an optional param.
- `sample_id` derivation was changed from baseName to simpleName.
- The runtime analysis (2026-03-12) came from the Nextflow traces in `htan_hms_crc_he_41_dual_mode/nextflow_reports/`: box plots plus connected dots paired by sample, normalised by cell count. Wall clock includes scheduling time. It became a slide (`plot_classpose_runtimes.py`).
- Local hacks (2026-01-06): `openslide-bin` installed into the venv. A spatialdata points-write bug ("Scale is not JSON serializable") survived pyarrow 17 and 18 pins and an upgrade to spatialdata 0.6.1, and was fixed by editing `site-packages/spatialdata/_io/io_points.py` in the `classpose/.venv`. An earlier edit to ClassPose `outputs.py` was reverted. These patches exist only in that local venv.

### Data

- Cohort: HTAN HMS atlas, 40 CRC patients and 41 H&E whole-slide images. The H&E is from the Orion platform (Level 2 OME-TIFF under `ORION202305_imaging_level_2/`). HTA7_919 has two near-serial sections, used as a reproducibility panel (Figure 1d). At one point only 39 samples appeared; relisting the bucket by file rather than by participant fixed the count.
- Slides were addressed by CRDC DRS GUID. Outputs are named `dg_4DFC_<GUID>`. The test example used was `dg.4DFC/624693b0-7e68-11ee-a75b-033941d3e6da`.
- Metadata in `htan_hms_crc_he_41_dual_mode/metadata/`: `participant.tsv`, `biospecimen.tsv`, `file.tsv`, `files_drs.csv`, plus BigQuery exports `bquxjob_4136c6be_19c48c7c3e6.csv` (therapy) and `bquxjob_6358239c_19c48c9d3ac.csv` (targeted sequencing, 11 CRC genes). Plots are keyed on "Data File ID" (for example HTA9_122_123). Cohort plots use the short participant ID (HTA1_1).
- BigQuery, not the portal, was used to find H&E slides whose adjacent biospecimen has multiplexed imaging. The deck names TNP SARDANA (CRC, CyCIF) and HTAPP CODEX; Adam judged this "not exactly great" for a pan-cancer approach.
- S3 output paths in the prompts:
  - `s3://htan-project-tower-bucket/nf-classpose/outputs/testing/` (test run, 2026-01-21)
  - `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/` (the main CoNIC plus ConSeP run; layout `dg_4DFC_<GUID>/{conic,consep}/..._spatialdata.zarr`, read with AWS profile `tower`)
  - `s3://htan-hms-crc-he-41/classpose_outputs/` appears only inside the Figures 2 and 3 plan text and is probably wrong
- Seqera samplesheet: a dataset in the Sage Seqera Platform workspace (`tower.sagebionetworks.org`), mirrored by `nf-classpose/test/tower_samplesheet.csv`.
- No Synapse IDs appear anywhere in the session records.
- Other inputs: `PCGA-01-0001-029-19762-00897BX.svs` (too big for the local demo), `CMU-1-Small-Region.svs`, `htan1_examples/1027015.ome.tiff`.

### Analysis

- Composition: stacked and 100% stacked bars per sample for each model, with different palettes and the same ordering. There is a CoNIC vs ConSeP total-cell scatter. Adam asked whether lymphocyte vs plasma counts were just a tissue-size effect. ConSeP was chosen as the primary model because it is trained on CRC and separates malignant from healthy epithelium.
- TLS (in `htan_hms_crc_he_41`): rolling-window and graph-based detection of lymphocyte and plasma cell aggregates. The first plasma/lymphocyte maturity reading was inverted relative to the literature and was replaced by a revised staging (size, density, compactness, plasma localisation; stages 1 to 4). A radius scan (50, 100, 200 µm) settled on 50 µm, with 250 µm FOV examples. The method is written up in `TLS_ANALYSIS_METHODS.md`.
- Morphometrics from contours were explored for malignant vs normal epithelium and T vs B cells. Adam raised confounding by the direction the section runs through the mucosa.
- Immune phenotypes (desert, excluded, inflamed and so on, based on Immunoscore): v1 used grid binning with cohort-relative percentiles. v2 used a KNN cell graph with absolute density thresholds (squidpy, parallel workers). There are compartment bars, paired tumour/margin/stroma panels and spatial maps.
- Model agreement: CoNIC vs ConSeP cell matching (centroid plus contour) produced a Venn diagram and a confusion matrix with a "Not matched" class. Many inflammatory cells were called connective. This was judged "not super promising", and the conclusion was to run ConSeP only.
- Misclassification: a Peyer's patch in normal mucosa was largely called malignant (seen in sample 934).
- Survival (`scripts/survival_analysis.py` and `.md`): days to last follow-up was flagged as risky. The 4-level phenotypes were underpowered (p about 0.11). Continuous inflammatory percentage predicted PFS (HR 0.44 per SD, p 0.003, stage-adjusted). Eosinophil density also helped (PFS HR 0.63, p 0.043). No significant cancer-specific survival metric was found. Age and sex were added as covariates. Adam asked directly whether this was p-hacking.
- Spatial clustering: DBSCAN immune cluster fraction for all cells and for tumour-associated immune cells only (to exclude Peyer's patches), univariate Moran's I, a DBSCAN parameter scan, and HDBSCAN via `fast_hdbscan` with both kept (`compute_spatial_features.py`, `dbscan_sensitivity.py`).
- KRAS interaction (`mutation_interaction.py`, `MUTATION_ANALYSIS.md`): no single mutation predicts survival. Tumour-proximal immune clustering × KRAS interaction for OS has HR 0.085 and p 0.032 (KRAS-mutant n 14). There is no interaction for PFS. Adam asked why the forest plot looked null while the KM curves looked strong. TCGA CRC was discussed as a validation cohort.
- mIF matching: see the cellRingMask section below. Ki67, panCK and E-cadherin were judged the best malignant markers, and PD-1 was bimodal in healthy but not malignant cells. Adam doubted the cell-to-cell approach and suggested comparing neighbourhoods in mIF and H&E separately.
- Literature: the Orion paper (Lin et al. 2023, Nature Cancer, PMID 37349501), Cell 2022 (S0092-8674(22)01571-9), a Nature Cancer paper (s43018-023-00576-1) and PMID 41475845. The cohort is very likely the Orion paper's 40-patient BWH cohort 1.

### Figures and deck

- Figure 1: cohort metadata heatmap, CoNIC vs ConSeP scatters (shared limits per panel, square), stacked bars for both models, and the HTA7_919 serial-section panel. Each sub-figure was saved separately as PNG only, and a caption md was written.
- Figure 2: hot vs cold samples showing H&E alone, filled contours, and zooms into dense malignant, inflammatory and stromal regions. The stroma zooms were first chosen badly and were then re-picked by cell-type density. The zarr outputs were used, not GeoJSON.
- Figure 3: 3a spatial overview mosaic of all samples (spatialdata-plot, tiny points), 3b density heatmaps, and TLS panels.
- On 2026-03-12: a "picking ClassPose back up" report, `classpose_meeting_deck.md` (full-sentence titles, mIF section last), `methods.md` (nf-classpose config plus scripts, venv versions), `build_slides.py` (Google Slides through the Workspace CLI), a 40-minute meeting agenda, and an image.sc forum reply pointing users to nf-classpose.

### PanNuke training demo

- Plan: retrain ClassPose on PanNuke (HuggingFace `RationAI/PanNuke`) as a candidate pan-cancer tumour/immune/stroma model for HTAN.
- ClassPose 0.1.0 from GitHub main ships without the `transforms/` subpackage. Those files were copied into the venv by hand.
- 2% subsample (72 train, 21 val) on M2 Max MPS: stopped on purpose at epoch 10 with val loss 3.540 (best). Best checkpoint at `checkpoints/demo/pannuke_demo/pannuke_demo_best`, about 70 to 90 s per epoch. Adam asked to test on a 2% sample of fold 3.
- Full run on GCP project htan-dcc with Vertex AI custom jobs in us-central1: g2-standard-4 with an L4, prebuilt PyTorch 2.1 container, package `gs://htan-dcc-pannuke/dist/pannuke-0.1.0.tar.gz`, data at `gs://htan-dcc-pannuke/data/train`, checkpoints at `gs://htan-dcc-pannuke/checkpoints`, 500 epochs, batch 8. Spot, H100 and A100 options were discussed, along with whether 500 epochs were needed (the paper used 200 for CoNIC and NuCLS). The job failed with a non-zero exit, and the next attempt switched to the nf-classpose container. Further failures followed on 2026-03-13. These records show no full training run finishing.

## Dead ends

- Hierarchical clustering on raw proportions: compositional data that sums to 100%, and Adam was "not sure the clustering adds much". Dropped from Figure 1. `plot_clustering_clr.py` may be a follow-up attempt.
- Phenotype v1 grid binning with cohort-relative thresholds: replaced by the v2 KNN graph with absolute thresholds.
- Categorical phenotype vs survival: underpowered, so continuous features were used instead.
- The first TLS maturity reading: inverted relative to the literature, replaced by the revised staging.
- Morphometric separation of cell subtypes: possibly confounded by section orientation. Not pursued.
- CoNIC vs ConSeP cell matching as validation: poor agreement on inflammatory cells, so the decision was ConSeP only.
- Cell-to-cell mIF matching: robustness questioned in favour of neighbourhood-level comparison.
- Adjacent-section mIF plus H&E across HTAN as pan-cancer ground truth: too little suitable data.
- Baking all models into the Docker image on a 14 GB runner: out of disk. Later solved for ConSeP.
- Artefact detection on by default: OOM, so it is now off.
- Branch claude/issue-3-20260125-2138 (PR #10, model subfolders without model names in filenames): closed unmerged. The merged code went the other way and puts the model name in filenames.
- Pinning pyarrow for the spatialdata bug: did not help, so site-packages was patched instead.
- PanNuke on Vertex AI: job failures, and repeated torch image builds on a non-CUDA Mac filled the disk. Paused.

## Which outputs are current

- `htan_hms_crc_he_41_dual_mode` is the main line: CoNIC plus ConSeP, 2026-01-27 onward, with the deck and methods. `htan_hms_crc_he_41` is the earlier CoNIC-only exploration (TLS, morphometrics, crops, the confusion-matrix prototype). It is still the only home of the TLS work.
- `classify_immune_phenotypes_v2.py` and `data/immune_phenotype_{metrics,classification}_v2.csv` are current. `methods.md`, `survival_analysis.py`, `plot_immune_compartment_bars.py`, `plot_phenotype_paired_panels.py` and the deck (`immune_phenotype_scatter_v2.png`) use v2. Checked against the working tree: `plot_malignant_mosaic.py`, `plot_region_overlay.py` and `plot_phenotype_spatial_maps.py` still read the v1 classification CSV. `compute_spatial_features.py` and `dbscan_sensitivity.py` import loader functions from the v1 module, so v1 cannot be deleted.
- Figures with v1 and v2 pairs: `immune_phenotype_distribution(_v2).png` and `immune_phenotype_scatter(_v2).png`. The `_v2` files are current.
- `scripts/match_mif.py` (dated 2026-03-12) was rebuilt from old Claude logs in session 00299338, after Adam noticed that the mIF figures had no script in `scripts/`.
- nf-classpose: checked with git and gh, GitHub main (ce699b6, PR #11 merge, 2026-01-27) contains all of claude/issue-3-20260125-2131 including `vips_save`. The local clone is on issue-3-2131 with a clean working tree, but its `origin/main` ref is stale (eab0ffb).

## The cellRingMask CSV

- File: `P37_S29_A24_C59kX_E15__at__20220106_014304_946511_cellRingMask.csv`, 813,143,745 bytes (about 775 MiB), in `htan_hms_crc_he_41_dual_mode/`.
- It is the Orion mIF single-cell table (cell-ring mask quantification, 16 markers, x/y centroids) for the same section as H&E sample "926". `match_mif.py` pairs it with `dg_4DFC_CD2C5DCB-9095-4E90-94F4-1FD37578830D`, which `files_drs.csv` maps to HTA7_926_1013 (biospecimen HTA7_926_8). Matching uses the nearest centroid within 10 px (about 3.25 µm), with the same pixel space assumed. The deck reports 1,211,706 matched cells.
- Origin: Adam moved it from `~/Downloads` on 2026-02-10. macOS quarantine metadata shows it was downloaded with Microsoft Edge on 2026-02-10 at about 10:11 UTC. No source URL was recorded, and no Synapse ID or portal link appears in any prompt, memory file or metadata file. The project memory says the matched mIF exists but the Sorger Lab cell type calls do not.

## nf-classpose claude/issue-3 branches

- Issue #3 "Bundle additional classpose models in Docker image" came from the runner-disk workaround (CoNIC only). The GitHub Claude app created three branches on 2026-01-25.
- claude/issue-3-20260125-2131: ConSeP model and model matrix. PR #11 into main was merged on 2026-01-27 and closed the issue. It also carries the 2026-01-26 work: optional id column, stub blocks, model name in filenames, `vips_save`. It was meant to merge and did.
- claude/issue-3-20260125-2155: an alternative take, merged into 2131 through PR #12, so it is in main.
- claude/issue-3-20260125-2138: PR #10 closed unmerged. One commit is not in main. Superseded.
- All three can be deleted once the 2138 commit is confirmed unwanted.

## Open questions

- Which HTAN or Synapse file the cellRingMask CSV is (likely an HMS Orion Level 4 mIF cell table for HTA7_926). Needs looking up before the file goes into the new Synapse project.
- Whether the outputs under `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/` still exist and should be archived.
- Whether anything finished in `gs://htan-dcc-pannuke/checkpoints`, and whether the bucket and the Vertex jobs should be cleaned up.
- Whether the spatialdata write bug also affects the nf-classpose container, which installs ClassPose from pip without the local patch. The dual-model run did produce zarrs.
- How the GeoJSON-only output, the GPU retry scaling and the runtime EfficientNet download were finally resolved.
- Whether the three plot scripts that still read the v1 phenotype CSV should move to v2.
- Whether the KRAS × immune clustering OS result survives correction for the many tests run. Adam raised p-hacking himself, and it needs an external cohort (TCGA CRC was suggested).
- Whether to get mIF cell type calls from the Sorger Lab or derive them independently.
