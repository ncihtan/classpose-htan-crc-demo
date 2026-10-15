# Pipeline

ClassPose ran on the 41 H&E slides through [nf-classpose](https://github.com/ncihtan/nf-classpose) on Seqera Platform, using AWS Batch GPU spot instances. This folder records what ran, so you can rerun it or check the outputs.

## What ran

| Item | Value |
|---|---|
| Pipeline | `ncihtan/nf-classpose` (moved from `adamjtaylor/nf-classpose` on 2026-10-15; old links redirect) |
| Commit used for the run | `dcb045250a6df043ba0b4eb5912d8b016855aeeb` (branch `claude/issue-3-20260125-2131`) |
| Same code on `main` | merge commit `ce699b6b1038d16eb8d998c916ce8847093514ff` (PR #11, merged 2026-01-27) |
| Models | CoNIC and CoNSeP (`--models conic,consep`) |
| Inputs | 41 HTAN HMS H&E OME-TIFFs, fetched by DRS URI with gen3-client |
| Outputs | `s3://htan-project-tower-bucket/nf-classpose/outputs/htan_hms_crc_he_41_dual_model/` (257 GB: the 460 files in `s3_manifest.tsv` plus 82 SpatialData zarrs, 22.6 GB). Checked on 2026-10-15: every manifest file is present with the same size. |
| Runs | 2026-01-26 16:04, 2026-01-27 10:03, 2026-01-27 11:56 (resumes; the last one is all cached except one task) |

## Files

| File | What it is |
|---|---|
| `samplesheet_htan_hms_crc_he_41.csv` | The samplesheet: one DRS URI per slide. Copied from `nf-classpose/test/`. |
| `params.reconstructed.yaml` | Our best reconstruction of the run params. The original params file was a temporary Seqera upload and is lost. |
| `s3_manifest.tsv` | Every output file in S3 with its HTAN IDs, size and date. Made by `analysis/scripts/generate_s3_manifest.sh`. |
| `reports/execution_trace_*.txt` | Nextflow traces for the three runs. `analysis/scripts/plot_classpose_runtimes.py` reads these. |
| `reports/pipeline_dag_*.svg` | Nextflow DAGs. |
| `reports/execution_report_*.html` | Nextflow execution reports. Copies are also in S3 next to the outputs. |
| `example_run/` | The first end-to-end test on one HTAN slide (2026-01-25, run `fabulous_faggin`): the params, the resolved config, the reports, the density CSV and a draft blog post about the pipeline. Its params are the closest record of the settings used for the main run. |

## Rerun

You need Gen3 credentials for the CRDC (HTAN open-access data), a Nextflow setup with a GPU, and somewhere to write about 250 GB.

```bash
nextflow run ncihtan/nf-classpose -r ce699b6 \
  -params-file params.reconstructed.yaml \
  --gen3_credentials ~/.gen3/credentials.json \
  --outdir <your-output-location> \
  -profile <your-profile>
```

To check the outputs, compare your output listing with `s3_manifest.tsv`.
