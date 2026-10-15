# Automated Cell Classification at Scale: Processing HTAN Whole Slide Images with nf-classpose

*January 25, 2026*

The Human Tumor Atlas Network (HTAN) is generating thousands of whole slide images (WSIs) that need to be analyzed for cell classification. Today, I'm excited to share nf-classpose, a Nextflow pipeline that automates WSI cell classification using [Classpose](https://github.com/sohmandal/classpose)—and how we ran it end-to-end on cloud infrastructure using DRS URIs to pull data directly from the NCI Cancer Research Data Commons.

## The Challenge

Processing whole slide images for cell classification involves several steps:
1. **Data acquisition** - Slides may be stored in federated data repositories with complex access controls
2. **Format conversion** - OME-TIFF files need conversion to OpenSlide-compatible formats
3. **Cell classification** - Running inference with deep learning models (ideally on GPUs)
4. **Output generation** - Producing spatial data in analysis-ready formats

Doing this manually for thousands of slides isn't feasible. We needed automation.

## nf-classpose: The Solution

nf-classpose wraps the Classpose cell classification tool in a Nextflow pipeline that handles the entire workflow:

```
Samplesheet → DRS Download → Format Conversion → Cell Classification → Outputs
```

The pipeline automatically:
- Downloads slides from Gen3/NCI CRDC using DRS URIs
- Converts OME-TIFF files to OpenSlide-compatible pyramidal TIFFs (preserving resolution metadata)
- Runs GPU-accelerated cell classification with the CoNIC model
- Outputs cell contours, centroids, density statistics, and SpatialData objects

## A Real Run: Processing an HTAN Slide

Let's walk through a real execution that processed a single HTAN slide from start to finish.

### Input

The input was simple—a CSV file with a single DRS URI:

```csv
slide_path
drs://nci-crdc.datacommons.io/dg.4DFC/624693b0-7e68-11ee-a75b-033941d3e6da
```

This URI points to a whole slide image stored in the NCI Cancer Research Data Commons.

### Execution

The pipeline ran on AWS Batch via Seqera Platform (Tower), orchestrating three processes:

| Process | Duration | Peak Memory | Description |
|---------|----------|-------------|-------------|
| `GEN3_DOWNLOAD` | 3m 23s | 18 MB | Downloaded the ~12 GB slide using gen3-client |
| `VIPS_CONVERT` | 15.2s | 342 MB | Converted OME-TIFF to pyramidal TIFF |
| `CLASSPOSE_PREDICT_WSI` | 4m 58s | 22.2 GB | GPU inference with CoNIC model |

Total wall-clock time: **~35 minutes** (including AWS Batch job scheduling overhead).

### The Pipeline in Action

Here's what each step accomplished:

**1. GEN3_DOWNLOAD**
The pipeline automatically configured gen3-client with HTAN credentials stored in S3, extracted the GUID from the DRS URI, and downloaded the slide. The container used Seqera Wave for just-in-time container provisioning.

**2. VIPS_CONVERT**
The downloaded OME-TIFF needed conversion because OpenSlide (used by Classpose) doesn't natively support OME-TIFF. VIPS created a pyramidal TIFF with:
- JPEG compression (configurable)
- 256x256 tiles
- Resolution metadata extracted from OME-XML (`PhysicalSizeX`)

This last point is critical—without correct resolution metadata, downstream tissue detection would fail due to incorrect scale calculations.

**3. CLASSPOSE_PREDICT_WSI**
The main inference ran on a GPU-enabled instance (32 GB memory, 1 GPU). The CoNIC model detected and classified cells, producing multiple output formats.

### Configuration

The pipeline ran with these key parameters:

```json
{
  "model_config": "conic",
  "tile_size": 1024,
  "overlap": 64,
  "batch_size": 8,
  "output_type": "csv spatialdata",
  "vips_compression": "jpeg"
}
```

## Results Analysis

The pipeline produced a structured output directory for each sample:

```text
example_run/
├── fabulous_faggin-params.json          # Run parameters
├── resolved.config                       # Full Nextflow config
├── logs/
│   ├── nf-3stNqO29UUEC1M.log            # Detailed execution log
│   ├── nf-3stNqO29UUEC1M.txt            # Summary output
│   └── timeline-3stNqO29UUEC1M.html     # Interactive timeline
├── pipeline_info/
│   ├── execution_report_*.html          # Resource usage report
│   ├── execution_trace_*.txt            # Per-task metrics
│   └── pipeline_dag_*.svg               # DAG visualization
└── dg_4DFC_624693b0-7e68-11ee-a75b-033941d3e6da/
    ├── *_cell_centroids.geojson         # Cell centroid points (615 KB)
    ├── *_cell_contours.geojson          # Cell polygon outlines (1.3 MB)
    ├── *_tissue_contours.geojson        # Tissue regions (108 KB)
    ├── *_cell_densities.csv             # Summary statistics
    ├── *_spatialdata.zarr/              # SpatialData object (600 KB)
    └── versions.yml                      # Software versions
```

### Cell Classification Results

The CoNIC model classified **1,392 cells** across **109 tissue regions**:

| Cell Type   | Count | Percentage | Avg Area (px²) |
| ----------- | ----- | ---------- | -------------- |
| Epithelial  | 636   | 45.7%      | 217.1          |
| Connective  | 433   | 31.1%      | 180.8          |
| Lymphocyte  | 199   | 14.3%      | 104.4          |
| Plasma cell | 120   | 8.6%       | 137.1          |
| Eosinophil  | 4     | 0.3%       | 219.5          |

The cell composition suggests this is likely tumor tissue, with nearly half being epithelial cells and a significant immune infiltrate (lymphocytes + plasma cells = 23%).

### Cell Density Statistics

The `_cell_densities.csv` provides per-region density calculations (cells per mm²):

```csv
region,cell_class,count,density
tissue,Epithelial,636,977.76
tissue,Connective,433,665.68
tissue,Lymphocyte,199,305.93
tissue,Plasma cell,120,184.48
tissue,Eosinophil,4,6.15
tissue,Neutrophil,0,0.0
```

These densities enable quantitative comparisons across samples—critical for studying immune infiltration patterns across the HTAN cohort.

### SpatialData Output

The pipeline outputs a [SpatialData](https://spatialdata.scverse.org/) Zarr store for seamless integration with the scverse ecosystem. The zarr metadata shows:

```json
{
  "mpp_x": 0.325,
  "mpp_y": 0.325,
  "model_config": "conic",
  "n_cells_detected": 1392
}
```

The resolution (0.325 µm/pixel) was correctly preserved from the original OME-TIFF metadata—essential for accurate cell area calculations and cross-sample comparisons.

The SpatialData object contains:

- **points/cell_centroids** - Cell locations as coordinate pairs
- **shapes/cell_contours** - Cell segmentation polygons
- **shapes/tissue_contours** - Tissue region boundaries
- **tables/cellular_densities** - AnnData table with density statistics

This format enables downstream analysis with scanpy, squidpy, and other spatial analysis tools without format conversion.

## Infrastructure

The run leveraged several cloud-native features:

- **AWS Batch** - Serverless compute with automatic scaling
- **Seqera Wave** - Container provisioning without pre-pushing to registries
- **Seqera Platform** - Orchestration, monitoring, and cost tracking
- **Cloud cache** - For resuming failed runs without recomputation

Resource labels tracked costs by project:
```groovy
resourceLabels = [
  CostCenter: '120100',
  Project: 'Human_Tumor_Atlas_Network',
  Department: 'SCCE'
]
```

## Running It Yourself

To process your own slides:

```bash
# Local files
nextflow run adamjtaylor/nf-classpose \
  --input samplesheet.csv \
  --outdir results \
  -profile docker,gpu

# DRS URIs from NCI CRDC
nextflow run adamjtaylor/nf-classpose \
  --input drs_samplesheet.csv \
  --gen3_credentials ~/.gen3/credentials.json \
  -profile docker,gpu
```

The pipeline handles the complexity—you just provide slide paths.

## What's Next

With nf-classpose, we can now process the entire HTAN slide collection systematically. The outputs integrate directly with downstream spatial analysis tools via the SpatialData format.

The code is open source: [github.com/adamjtaylor/nf-classpose](https://github.com/adamjtaylor/nf-classpose)

---

*This work is part of the Human Tumor Atlas Network (HTAN), supported by the National Cancer Institute.*
