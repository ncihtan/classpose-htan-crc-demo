# HTAN Pan-Cancer Cell Phenotyping Study Plan

## Overview

Large-scale analysis of ClassPose cell phenotyping outputs across thousands of whole slide images from the Human Tumor Atlas Network (HTAN), spanning multiple cancer types.

## Data Inventory

### Expected ClassPose Outputs per Slide
- `{slide_id}_cell_contours.geojson` - Cell polygons with classifications
- `{slide_id}_cell_centroids.geojson` - Cell center points
- `{slide_id}_cell_densities.csv` - Per-slide density statistics
- `{slide_id}_tissue_contours.geojson` - Tissue regions
- `{slide_id}_artefact_contours.geojson` - Artifact masks
- `{slide_id}_spatialdata.zarr/` - Integrated SpatialData format

### HTAN Metadata to Integrate
- Cancer type (BRCA, LUAD, CRC, PDAC, etc.)
- Tumor stage and grade
- Patient demographics
- Treatment history
- Molecular subtypes (PAM50, CMS, etc.)
- Survival outcomes
- Spatial region annotations (tumor core, invasive margin, stroma)

---

## Phase 1: Data Aggregation & QC

### 1.1 Build Unified Database
```
htan_cellpose_db/
├── metadata/
│   ├── slide_manifest.parquet      # All slides with HTAN IDs
│   ├── patient_clinical.parquet    # Clinical data
│   └── qc_metrics.parquet          # Per-slide QC scores
├── cell_data/
│   ├── all_cells.parquet           # Aggregated cell-level data
│   └── by_cancer_type/             # Partitioned by cancer
└── spatial_features/
    └── per_slide_features.parquet  # Computed spatial metrics
```

### 1.2 Quality Control Pipeline
- Filter slides with >20% artifact area
- Exclude slides with <500 detected cells
- Flag slides with unusual cell type distributions (outlier detection)
- Validate against pathologist annotations where available

### 1.3 Normalization
- Adjust cell densities for tissue area
- Account for staining batch effects
- Harmonize across imaging centers

---

## Phase 2: Descriptive Analysis

### 2.1 Pan-Cancer Cell Composition
**Question**: How does cellular composition vary across cancer types?

| Metric | Description |
|--------|-------------|
| Cell type proportions | % Epithelial, Lymphocyte, Connective, Plasma, Neutrophil, Eosinophil |
| Absolute densities | Cells/mm² by type |
| Heterogeneity index | Shannon diversity of cell types per slide |

**Visualizations**:
- Stacked bar charts of cell composition by cancer type
- UMAP embedding of slides based on cell composition
- Heatmap of cell type correlations

### 2.2 Spatial Distribution Patterns
**Question**: Do cell types show consistent spatial organization?

| Metric | Description |
|--------|-------------|
| Clustering coefficient | Spatial autocorrelation (Moran's I) |
| Nearest neighbor distances | Mean distance to same/different cell types |
| Ripley's K function | Multi-scale spatial clustering |
| Tumor-stroma interface | Cell composition at boundaries |

**Visualizations**:
- Spatial density maps aggregated by cancer type
- Violin plots of spatial metrics by cancer type
- Representative slide galleries

---

## Phase 3: Hypothesis-Driven Analysis

### Hypothesis 1: Immune Infiltration Signatures
**H1**: Lymphocyte density and spatial distribution correlate with molecular subtypes

- Compare TIL density across PAM50 subtypes (breast)
- Compare TIL density across CMS subtypes (colorectal)
- Correlate with immune gene signatures (if RNA-seq available)

### Hypothesis 2: Tumor Microenvironment Architecture
**H2**: Spatial organization of immune cells predicts patient outcomes

- Define "immune hot" vs "immune cold" vs "immune excluded" phenotypes
- Test association with progression-free survival
- Identify cancer-specific vs universal patterns

### Hypothesis 3: Stromal Composition
**H3**: Connective tissue density reflects tumor aggressiveness

- Correlate with tumor grade
- Compare primary vs metastatic samples
- Association with treatment response

### Hypothesis 4: Inflammatory Signatures
**H4**: Neutrophil/Eosinophil infiltration indicates distinct biology

- Correlate with mutation burden
- Association with specific driver mutations
- Prognostic value by cancer type

### Hypothesis 5: Cell Morphometry
**H5**: Cell size/shape features carry diagnostic/prognostic information

- Compare cell area distributions across grades
- Nuclear-to-cytoplasm ratio proxies
- Morphometric heterogeneity within tumors

---

## Phase 4: Advanced Spatial Analysis

### 4.1 Neighborhood Analysis
- Define local neighborhoods (50μm, 100μm, 200μm radii)
- Compute cell type mixing scores
- Identify recurring spatial motifs (e.g., immune aggregates)

### 4.2 Graph-Based Analysis
- Construct cell proximity graphs (Delaunay triangulation)
- Community detection for cellular neighborhoods
- Graph neural network features

### 4.3 Spatial Statistics
```python
# Key spatial metrics to compute
metrics = {
    'morans_i': per_cell_type_spatial_autocorrelation,
    'ripleys_k': multi_scale_clustering,
    'cross_k': cell_type_co_localization,
    'hotspot': local_density_clusters,
    'interface_score': tumor_stroma_boundary_composition
}
```

---

## Phase 5: Machine Learning

### 5.1 Outcome Prediction
- **Input**: Per-slide feature vectors (composition + spatial metrics)
- **Targets**: Overall survival, treatment response, recurrence
- **Models**: Random Forest, XGBoost, survival models (Cox-PH, DeepSurv)

### 5.2 Cancer Type Classification
- Can cell phenotyping alone distinguish cancer types?
- Which features are most discriminative?
- Transfer learning across cancer types

### 5.3 Subtype Discovery
- Unsupervised clustering of slides by cellular features
- Do data-driven subtypes align with molecular subtypes?
- Novel subtype identification

---

## Phase 6: Visualization & Reporting

### 6.1 Interactive Dashboard
- Slide browser with cell overlays
- Filterable by cancer type, patient characteristics
- Drill-down from cohort to individual slides

### 6.2 Publication Figures
- Pan-cancer overview (n slides, n cells, n patients)
- Cell composition summary by cancer type
- Key hypothesis results with survival curves
- Representative image galleries

### 6.3 Data Release
- Processed cell tables (cell ID, x, y, type, slide_id)
- Per-slide summary statistics
- Spatial feature matrices
- Trained models and code

---

## Technical Implementation

### Compute Requirements
| Task | Resources |
|------|-----------|
| Cell aggregation | 32GB RAM, parallel I/O |
| Spatial metrics | GPU for large slides |
| ML training | GPU cluster |
| Storage | ~100GB for processed data |

### Key Libraries
- `spatialdata` - Data management
- `squidpy` - Spatial analysis
- `scanpy` - Downstream analysis
- `lifelines` - Survival analysis
- `scikit-learn` / `xgboost` - ML models
- `plotly` / `altair` - Interactive viz

### Pipeline Structure
```
01_aggregate_cells.py      # Combine all cell data
02_compute_spatial.py      # Calculate spatial metrics
03_qc_and_filter.py        # Quality control
04_descriptive_stats.py    # Summary statistics
05_hypothesis_tests.py     # Statistical tests
06_ml_models.py            # Predictive modeling
07_generate_figures.py     # Publication figures
```

---

## Timeline & Milestones

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| Data aggregation | 2 weeks | Unified database |
| QC & normalization | 1 week | Clean dataset |
| Descriptive analysis | 2 weeks | Summary report |
| Hypothesis testing | 4 weeks | Statistical results |
| ML modeling | 3 weeks | Predictive models |
| Visualization | 2 weeks | Dashboard + figures |
| **Total** | **14 weeks** | Publication-ready analysis |

---

## Expected Outcomes

1. **Pan-cancer cell atlas**: Comprehensive characterization of cellular composition across HTAN cancer types
2. **Spatial biomarkers**: Novel spatial metrics with prognostic value
3. **Predictive models**: ML models for outcome prediction from histology alone
4. **Resource for community**: Processed data and code for HTAN consortium
