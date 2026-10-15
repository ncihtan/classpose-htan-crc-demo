# Tertiary Lymphoid Structure (TLS) Analysis Methods

## Scientific Rationale and Approach

### Background

Tertiary Lymphoid Structures (TLS) are ectopic lymphoid aggregates that form in chronically inflamed tissues, including tumors. In cancer, TLS are associated with favorable prognosis and response to immunotherapy, as they serve as sites for local antigen presentation, B-cell activation, and antibody production.

TLS develop through a maturation process analogous to secondary lymphoid organ development:

1. **Early/Immature TLS**: Loose lymphoid aggregates without clear organization
2. **Primary follicle-like**: Denser B-cell clusters beginning to organize
3. **Secondary follicle-like (Mature)**: Organized structures with germinal centers (GC), featuring follicular dendritic cell networks, active B-cell proliferation, and class-switch recombination
4. **Post-GC/Involuted**: Structures where germinal center activity has ceased, characterized by plasma cell accumulation and reduced lymphocyte density

### Problem with Previous Interpretation

The initial analysis used plasma cell/lymphocyte (P/L) ratio as a direct proxy for TLS maturity, with the assumption that higher P/L ratio = more mature TLS. This interpretation was **inverted** relative to TLS biology:

| Feature | Literature Consensus | Previous (Incorrect) Interpretation |
|---------|---------------------|-------------------------------------|
| **Immature TLS** | Small, loose lymphoid aggregates | Large, actively recruiting |
| **Mature TLS** | Large, organized with active GC | Small, high P/L ratio |

**Key biological insight**: In cancer (chronic antigen exposure), mature TLS are **large, organized structures** with active germinal centers. High P/L ratios in **small** clusters likely represent:

- Post-germinal center remnants where B-cell differentiation to plasma cells has completed
- Involuted structures where the GC reaction has ceased
- Plasma cell-rich zones adjacent to (not within) active TLS

### Data Constraints

This analysis uses cell centroid coordinates and cell type classifications from H&E-based deep learning inference. We do **not** have access to:

- Follicular dendritic cell (FDC) markers (CD21, CD23)
- Germinal center markers (BCL6, Ki67)
- B-cell subset markers (CD20, CD27, IgD)
- T-cell zone markers (CD3, CD4)

Therefore, we must infer TLS maturity from spatial organization patterns rather than molecular markers.

---

## Analytical Approach

### 1. Immune Cell Clustering

TLS are identified as spatially connected clusters of immune cells (lymphocytes + plasma cells) using graph-based clustering:

```
1. Build KD-tree from all immune cell coordinates
2. Connect cells within 50-pixel radius (approximate cell-cell contact distance)
3. Extract connected components as candidate TLS clusters
4. Filter: minimum 15 cells, at least 10 lymphocytes and 3 plasma cells
```

The 50-pixel radius was chosen to capture direct cell-cell interactions while avoiding merger of distinct aggregates.

### 2. Multi-Metric Feature Extraction

For each cluster, we calculate:

#### Size Metrics
- **n_cells**: Total immune cells in cluster
- **n_lymph**: Lymphocyte count
- **n_plasma**: Plasma cell count
- **P/L ratio**: n_plasma / n_lymph

#### Spatial Organization Metrics
- **Area**: Convex hull area (pixels²)
- **Perimeter**: Convex hull perimeter
- **Density**: n_cells / area (cells/pixel²)
- **Compactness**: 4π × area / perimeter² (circularity, range 0-1)

#### Plasma Cell Localization
- **plasma_dist_ratio**: mean(plasma cell distance to centroid) / mean(lymphocyte distance to centroid)

This metric captures whether plasma cells are:
- **Peripheral (ratio > 1.0)**: Consistent with plasma cells emigrating from an active GC
- **Central (ratio < 1.0)**: Plasma cells concentrated in cluster core
- **Mixed (ratio ≈ 1.0)**: No clear spatial organization

### 3. Literature-Aligned TLS Staging

Based on the biological understanding that mature TLS are large and organized:

| Stage | Name | Size Criteria | P/L Ratio | Biological Interpretation |
|-------|------|---------------|-----------|---------------------------|
| 1 | **Early** | < 100 cells | Variable | Nascent lymphoid aggregates, no GC formation |
| 2 | **Developing** | 100-500 cells | Variable | Primary follicle-like, increasing organization |
| 3 | **Mature/Active** | > 500 cells | Low-moderate | Secondary follicle-like with active GC |
| 4 | **Involuted** | < 200 cells | > 0.4 | Post-GC remnants, exhausted structures |

**Stage 4 (Involuted)** is specifically defined by the combination of small size AND high P/L ratio, capturing structures where:
- GC activity has ceased
- Most B-cells have differentiated to plasma cells
- The lymphocyte pool is depleted

### 4. Composite Maturity Score

A continuous maturity score integrates multiple features:

```python
maturity_score = (
    0.50 × normalized_log_size +      # Larger = more mature
    0.25 × organization_score +        # Compact, organized = more mature
    0.25 × plasma_periphery_score      # Peripheral plasma = active GC output
)
```

**Weighting rationale**:
- **Size (50%)**: Dominant factor; mature TLS require sustained immune recruitment
- **Organization (25%)**: Compactness indicates structural organization; penalizes high P/L in small clusters
- **Plasma periphery (25%)**: Only scored for clusters ≥50 cells; peripheral plasma cells suggest active GC with emigrating antibody-secreting cells

**Organization score adjustment**:
```python
# Penalize high P/L ratio in small clusters (involution marker)
size_factor = min(1.0, n_cells / 500)
pl_penalty = (1 - size_factor) × max(0, pl_ratio - 0.3)  # Only if small
organization_score = compactness - pl_penalty
```

This ensures small clusters with high P/L ratios receive lower maturity scores despite potentially high density.

---

## Biological Rationale for Key Decisions

### Why Size is the Primary Maturity Indicator

1. **GC formation requires critical mass**: Germinal center reactions require sufficient B-cells, T follicular helper cells, and FDCs to establish productive interactions
2. **Sustained recruitment indicates active immunity**: Large TLS reflect ongoing immune cell recruitment, suggesting active anti-tumor responses
3. **Clinical correlations**: Studies show larger TLS correlate with better prognosis and immunotherapy response

### Why High P/L Ratio in Small Clusters Indicates Involution

1. **Plasma cells are terminally differentiated**: High plasma cell proportions indicate B-cell differentiation has occurred
2. **Lymphocyte depletion**: Small size with high P/L suggests the lymphocyte pool is not being replenished
3. **Post-GC phenotype**: After GC shutdown, plasma cells persist while naive/memory B-cells disperse

### Why Plasma Cell Localization Matters

In active germinal centers:
1. B-cells undergo somatic hypermutation and selection in the GC dark/light zones
2. Selected B-cells differentiate to plasma cells
3. Plasma cells **emigrate** from the GC toward surrounding tissue or bone marrow

Therefore, peripheral plasma cell localization (ratio > 1.0) is consistent with active GC output, while central/mixed distribution may indicate disorganized or involuted structures.

---

## Limitations and Caveats

### Technical Limitations

1. **No molecular markers**: Cannot directly identify GC (BCL6+, Ki67+), FDC networks, or B/T cell zones
2. **Cell type inference from H&E**: Classification accuracy depends on deep learning model performance
3. **2D analysis**: TLS are 3D structures; single-section analysis may miss organization patterns
4. **Fixed radius for clustering**: 50-pixel threshold may not be optimal for all tissue regions

### Biological Limitations

1. **Plasma cell localization is indirect**: Peripheral plasma cells could also reflect:
   - Pre-existing plasma cells in adjacent tissue
   - Plasma cells from other sources migrating through the region

2. **Size thresholds are approximations**: The 100/200/500 cell thresholds are based on general principles, not validated cut-offs for this specific tissue/cancer type

3. **Temporal dynamics not captured**: A single time-point cannot distinguish:
   - Developing TLS that will become mature
   - Stable early-stage aggregates
   - Recently involuted structures vs. long-quiescent ones

### Recommendations for Validation

1. **Multiplex IHC/IF**: Validate staging with CD20, CD21, BCL6, Ki67 staining
2. **Serial sections**: Assess 3D structure consistency
3. **Clinical correlation**: Compare maturity scores with patient outcomes
4. **Cross-cohort validation**: Test staging criteria in independent datasets

---

## Output Interpretation

### Stage Distribution (Example Results)

```
Stage 1 (Early):        61.5%  - Nascent lymphoid aggregates
Stage 2 (Developing):    9.5%  - Growing, organizing structures
Stage 3 (Mature/Active): 2.4%  - Large, organized TLS with likely GC activity
Stage 4 (Involuted):    26.6%  - Post-GC remnants
```

**Interpretation**: Most immune aggregates are small/early stage. A subset (2.4%) reaches mature/active status with large, organized structures. Notably, ~27% show involuted phenotype (small + high P/L), suggesting prior GC activity that has since ceased.

### Maturity Score Usage

- **Research**: Correlate continuous maturity scores with gene expression, mutation burden, or clinical outcomes
- **Clinical**: Use stage classification for categorical analyses (e.g., "TLS-high" vs "TLS-low" tumors)
- **Spatial analysis**: Map maturity scores to identify tumor regions with active vs. exhausted immune responses

---

## References

1. Sautès-Fridman C, et al. (2019). Tertiary lymphoid structures in the era of cancer immunotherapy. *Nature Reviews Cancer*. 19:307-325.

2. Helmink BA, et al. (2020). B cells and tertiary lymphoid structures promote immunotherapy response. *Nature*. 577:549-555.

3. Cabrita R, et al. (2020). Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature*. 577:561-565.

4. Fridman WH, et al. (2022). B cells and tertiary lymphoid structures as determinants of tumour immune contexture and clinical outcome. *Nature Reviews Clinical Oncology*. 19:441-457.

5. Schumacher TN, Thommen DS. (2022). Tertiary lymphoid structures in cancer. *Science*. 375:eabf9419.
