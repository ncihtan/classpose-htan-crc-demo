# Tumor Immune Phenotype Classification: Figures and Methods

## Figure Captions

### Figure 1: Immune Phenotype Distribution
**File:** `immune_phenotype_distribution.png`

Distribution of tumor immune phenotypes across the HTAN CRC cohort (n=41 samples). Samples were classified into four phenotypes based on spatial analysis of cell distributions: Hot (n=9, 22%), Altered-immunosuppressed (n=9, 22%), Altered-excluded (n=18, 44%), and Immune-desert (n=5, 12%). Classification used cohort-relative percentile thresholds for inflammatory cell content and tumor infiltration ratio.

---

### Figure 2: Immune Phenotype Classification Scatter Plot
**File:** `immune_phenotype_scatter.png`

Two-axis classification of immune phenotypes based on overall inflammatory cell percentage (x-axis) and infiltration ratio (y-axis). The infiltration ratio represents the density of inflammatory cells in tumor core regions relative to stromal regions, where values >1 indicate preferential tumor infiltration. Dashed lines indicate cohort-derived percentile thresholds: inflammatory P33 (12.1%) and P67 (18.1%) for low/medium/high inflammatory content, and infiltration P33 (0.39) for poor/good tumor infiltration. Points are colored by assigned phenotype.

---

### Figure 3: Spatial Cell Type Distribution by Immune Phenotype
**File:** `immune_phenotype_spatial_comparison.png`

Representative spatial maps showing individual cell locations for one sample from each immune phenotype. Each point represents a single cell centroid from ConSeP classification, colored by cell type: tumor cells (dark blue), inflammatory/immune cells (magenta), stromal cells (mustard), and other tissue including healthy epithelium and muscle (grey).

**Panels:**
- **(A) Hot** (HTA7_941, 32.9% immune): Dense immune cell infiltration throughout tumor regions, with magenta cells intermixed with blue tumor cells.
- **(B) Altered-immunosuppressed** (HTA7_983, 18.0% immune): Moderate immune presence with cells distributed across the tissue but without dense clustering.
- **(C) Altered-excluded** (HTA7_917, 38.3% immune): High overall immune content but inflammatory cells concentrated in stromal regions rather than infiltrating tumor core.
- **(D) Immune-desert** (HTA7_972, 11.3% immune): Sparse immune cell presence with tumor and stromal cells dominating the tissue architecture.

---

### Figure 4: Tumor Microenvironment Region Assignment
**File:** `region_assignment_comparison.png`

Spatial visualization of tumor microenvironment region assignments overlaid on tissue morphology. Regions were defined by 2D histogram binning of cell densities (75×75 grid) and classified as: tumor core (red, top 25% malignant epithelial density), tumor margin (orange, bins adjacent to core with >5% malignant cells), and stroma (green, <10% malignant and >30% stromal cells). Semi-transparent overlays (alpha=0.5) are shown on greyscale H&E images where available, or on grey cell scatter plots. Uncolored regions represent other tissue types (healthy epithelium, muscle, sparse areas).

**Panels:**
- **(A) Hot** (HTA7_941): Tumor core regions intermixed with stroma, indicating accessible tumor-stroma interface.
- **(B) Altered-immunosuppressed** (HTA7_983): Scattered tumor core regions with extensive stromal areas.
- **(C) Altered-excluded** (HTA7_917): Dense tumor core regions with distinct stromal boundary (shown with H&E background).
- **(D) Immune-desert** (HTA7_972): Large tumor core regions with surrounding stroma.

---

### Figure 5: Individual Cell Type Spatial Maps
**Files:** `spatial_map_hot.png`, `spatial_map_altered_immunosuppressed.png`, `spatial_map_altered_excluded.png`, `spatial_map_immune_desert.png`

High-resolution spatial cell type maps for representative samples of each immune phenotype. Cell centroids from ConSeP classification are plotted with cell type-specific colors: malignant epithelial/tumor (dark blue, #08306b), inflammatory/immune (magenta, #c51b8a), stroma (mustard, #d4a017), and other tissues including healthy epithelium and muscle (grey, #737373). Plot order ensures immune cells are rendered on top for visibility.

---

### Figure 6: Individual Region Overlay Maps
**Files:** `region_overlay_hot.png`, `region_overlay_altered_immunosuppressed.png`, `region_overlay_altered_excluded.png`, `region_overlay_immune_desert.png`

High-resolution tumor microenvironment region assignments for each phenotype (100×100 grid). Overlays show tumor core (red), tumor margin (orange), and stroma (green) with bin counts indicated in legend. These maps illustrate the spatial relationship between tumor and stromal compartments used for calculating infiltration ratios.

---

## Methods

### Cell Classification

Cell segmentation and classification were performed using the ConSeP (Colon Segmentation and Phenotyping) model as part of the nf-classpose pipeline. ConSeP classifies cells into six categories based on H&E morphology:
- Malignant epithelial (tumor)
- Inflammatory (immune cells including lymphocytes, plasma cells, macrophages)
- Stroma (fibroblasts, connective tissue cells)
- Healthy epithelial
- Muscle
- Other

Cell centroids and classifications were stored in SpatialData zarr format.

### Regional Definition

Tumor microenvironment regions were defined using 2D histogram binning of cell densities:

1. **Grid construction**: Each slide was divided into a 50×50 bin grid based on cell centroid coordinates.

2. **Per-bin cell counting**: For each bin, cells were counted by classification type.

3. **Region assignment**:
   - **Tumor core**: Bins in the top 25th percentile of malignant epithelial cell density (among bins containing malignant cells)
   - **Tumor margin**: Bins adjacent to core with >5% malignant cells
   - **Stroma-only**: Bins with <10% malignant cells and >30% stromal cells
   - **Other**: Remaining bins (healthy tissue, muscle, sparse regions)

### Spatial Metrics

For each sample, the following metrics were computed:

| Metric | Definition |
|--------|------------|
| `inflammatory_pct_overall` | Total inflammatory cells / total cells × 100 |
| `inflammatory_density_core` | Inflammatory cells in core bins / core area |
| `inflammatory_density_stroma` | Inflammatory cells in stroma bins / stroma area |
| `infiltration_ratio` | Core density / stroma density |

The infiltration ratio quantifies whether inflammatory cells preferentially localize to tumor regions (ratio >1) or are excluded to stromal regions (ratio <1).

### Phenotype Classification

Samples were classified using a two-axis system based on cohort-relative percentile thresholds:

**Axis 1: Inflammatory Content**
- Low: < P33 of cohort (12.1%)
- Medium: P33-P67 (12.1-18.1%)
- High: > P67 (18.1%)

**Axis 2: Tumor Infiltration**
- Poor: < P33 of infiltration ratio (0.39)
- Good: ≥ P33 (0.39)

**Classification Rules:**

| Inflammatory | Infiltration | Phenotype |
|--------------|--------------|-----------|
| Low | Poor | Immune-desert |
| Low | Good | Altered-excluded |
| High | Good | Hot |
| Medium/High | Poor | Altered-excluded |
| Medium | Good | Altered-immunosuppressed |

This percentile-based approach ensures classifications are relative to the cohort distribution rather than dependent on arbitrary absolute thresholds, improving robustness across datasets.

### Phenotype Definitions

| Phenotype | Characteristics | Interpretation |
|-----------|-----------------|----------------|
| **Hot** | High inflammatory content with good tumor infiltration | Active anti-tumor immune response; potential immunotherapy responders |
| **Altered-excluded** | Inflammatory cells present but excluded from tumor | Physical or functional barrier preventing immune infiltration |
| **Altered-immunosuppressed** | Moderate inflammation with infiltration but no robust response | Immune cells present but functionally suppressed |
| **Immune-desert** | Low inflammatory content and poor infiltration | Lack of immune recognition or recruitment |

### Visualization

Spatial cell type maps were generated by plotting individual cell centroids colored by ConSeP classification:
- Malignant epithelial: dark blue (#08306b)
- Inflammatory: magenta (#c51b8a)
- Stroma: mustard (#d4a017)
- Other tissues: grey (#737373)

Cells were plotted in order (other → stroma → tumor → immune) to ensure immune cells remain visible when overlapping with other cell types.

### Software

- Cell classification: ConSeP via nf-classpose pipeline
- Spatial data handling: SpatialData (v0.6.1)
- Clustering: scikit-learn DBSCAN (supplementary metric)
- Visualization: matplotlib
- Analysis: Python 3.13, pandas, numpy, scipy

### Limitations

1. **No functional markers**: ConSeP classification is based on H&E morphology alone; functional states (e.g., T cell exhaustion, cytokine production) cannot be assessed.

2. **Inflammatory cell subtypes**: The "Inflammatory" class encompasses all immune cell types (T cells, B cells, macrophages, etc.) without subtype discrimination.

3. **Density-based regions**: Tumor boundaries are approximated by cell density rather than pathologist annotation, which may not perfectly capture true tumor margins.

4. **Cohort-specific thresholds**: Percentile-based classification means phenotype assignments are relative to this cohort and may not directly compare to other studies using absolute thresholds.
