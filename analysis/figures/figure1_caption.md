# Figure 1: Automated Cell Segmentation Reveals Striking Tumor Microenvironment Heterogeneity in HTAN Colorectal Cancer

## Caption

**Figure 1. Large-scale automated cell segmentation of HTAN colorectal cancer H&E images reveals dramatic heterogeneity in tumor microenvironment composition.**

**(A) Cohort overview.** Clinical and demographic characteristics of the HTAN HMS CRC cohort. Heatmap shows patient age at diagnosis, sex, primary tumor site, AJCC pathologic stage, tumor grade, vital status, and disease recurrence for each H&E image (n = 41 images from 40 patients). Samples are ordered by inflammatory cell percentage (low to high, matching panels C-D). One patient (HTA7_919) contributed two near-serial sections, providing a technical replicate for assessing cell calling reproducibility. The cohort includes patients with a range of disease stages (I-IV), tumor grades, and outcomes.

**(B) Dual-model validation confirms robust cell classification.** Scatter plots comparing cell counts between ConSeP and CoNIC segmentation models for three biologically relevant categories: immune/inflammatory cells (ConSeP: Inflammatory vs. CoNIC: Lymphocyte + Plasma cell + Neutrophil + Eosinophil), epithelial cells (ConSeP: Malignant + Healthy epithelial vs. CoNIC: Epithelial), and stromal cells (ConSeP: Stroma vs. CoNIC: Connective). Strong correlations across all categories (Immune: r = 0.73; Epithelial: r = 0.90; Stromal: r = 0.82; all p < 10⁻⁷) validate that findings are robust to model choice despite different classification schemas. Each point represents one H&E image (n = 41). Dashed lines indicate identity (y = x); solid lines show linear regression fit.

**(C) CRC tumors exhibit striking heterogeneity in cellular composition.** 100% stacked bar chart showing the proportion of each cell type per H&E image, ordered by total inflammatory cell percentage. Inflammatory content varies substantially across samples, revealing a spectrum from "immune-cold" (left) to "immune-hot" (right) tumors. The balance between malignant epithelial (tumor), healthy epithelial, stromal, and muscle compartments also varies dramatically, suggesting distinct tumor microenvironment phenotypes with potential clinical implications for immunotherapy response.

**(D) Hierarchical clustering reveals distinct TME phenotypes.** Heatmap of centered log-ratio (CLR) transformed cell type proportions with hierarchical clustering (Ward's method, Euclidean distance). CLR transformation appropriately handles the compositional nature of proportion data, where values must sum to 100%. Colors indicate enrichment (red) or depletion (blue) relative to the geometric mean across cell types. Samples cluster into groups characterized by distinct cellular compositions, suggesting biologically meaningful TME subtypes rather than random variation.

---

## Methods Summary

Cell segmentation was performed using ClassPose with two pre-trained models: ConSeP (trained on the Colorectal Adenocarcinoma dataset; 6 cell classes) and CoNIC (trained on the CoNSeP + Lizard datasets; 6 cell classes). For the main analysis (Panels A, C, D), we used the ConSeP model as it was specifically trained on colorectal cancer tissue, making it the most appropriate model for this CRC cohort. Panel B demonstrates agreement between models by mapping cell classes to three broad biological categories, validating that findings are robust to model choice.

The ConSeP model distinguishes:
- **Malignant epithelial**: Tumor cells
- **Healthy epithelial**: Normal epithelial cells
- **Stroma**: Connective tissue fibroblasts
- **Muscle**: Smooth muscle cells
- **Inflammatory**: Immune/inflammatory infiltrate
- **Other**: Unclassified cells

CLR transformation in Panel D addresses the compositional data constraint where proportions sum to unity, avoiding spurious correlations inherent in raw proportion analysis.

## Data

- **Cohort**: HTAN HMS CRC (Human Tumor Atlas Network, Harvard Medical School Colorectal Cancer)
- **H&E images**: 41 (from 40 patients; HTA7_919 has 2 near-serial sections)
- **Primary model**: ConSeP (trained on colorectal adenocarcinoma)
- **Total cells analyzed**: 50,191,395 (ConSeP model)
- **Cell classes (ConSeP)**: Malignant epithelial, Healthy epithelial, Stroma, Muscle, Inflammatory, Other
- **Cell classes (CoNIC)**: Epithelial, Connective, Lymphocyte, Plasma cell, Neutrophil, Eosinophil

## Panel A: Clinical Variables

| Variable | Categories | Description |
|----------|------------|-------------|
| Age | <50, 50-59, 60-69, 70+ | Age at diagnosis (years) |
| Sex | Male, Female | Patient sex |
| Site | R. Colon, T. Colon, L. Colon, Sigmoid, Rectum, Appendix | Primary tumor site |
| Stage | I, II, III, IV | AJCC pathologic stage |
| Grade | Low, High | Tumor grade |
| Vital | Alive, Dead | Vital status at last follow-up |
| Recurrence | Yes, No | Disease progression or recurrence |

## Technical Replicate (HTA7_919)

One patient (HTA7_919) has two H&E images from near-serial sections, enabling assessment of cell calling reproducibility:

| Cell Type | Section 1 | Section 2 | Difference |
|-----------|-----------|-----------|------------|
| Inflammatory | 19.2% | 18.1% | 1.1% |
| Stroma | 34.9% | 33.8% | 1.1% |
| Muscle | 17.1% | 21.8% | 4.7% |
| Malignant epithelial | 23.0% | 15.8% | 7.2% |
| Healthy epithelial | 3.7% | 7.0% | 3.3% |

The inflammatory and stromal proportions show excellent concordance between serial sections. Greater variation in epithelial proportions likely reflects true biological heterogeneity within the tumor, as serial sections sample slightly different tissue regions.

## Color Keys

### Panel A: Clinical Variables
| Variable | Color Scheme |
|----------|--------------|
| Age | Yellow → Red (younger → older) |
| Sex | Blue (Male), Pink (Female) |
| Site | Categorical colors by anatomic location |
| Stage | Light → Dark blue (I → IV) |
| Grade | Green (Low), Red (High) |
| Vital | Green (Alive), Red (Dead) |
| Recurrence | Green (No), Red (Yes) |

### Panels C-D: Cell Types
| Cell Type | Color | Hex Code |
|-----------|-------|----------|
| Malignant epithelial | Red | #e41a1c |
| Healthy epithelial | Light red | #ff7f7f |
| Stroma | Blue | #377eb8 |
| Muscle | Teal | #66c2a5 |
| Inflammatory | Green | #4daf4a |
| Other | Gray | #999999 |

## Statistical Notes

- Panel B correlations: Pearson correlation coefficients with two-tailed p-values (n = 41 images)
- Panel D clustering: Ward's minimum variance method with Euclidean distance on CLR-transformed data
- CLR transformation: log(x / geometric_mean(x)), with pseudocount replacement for zeros

## Model Selection Rationale

The ConSeP model was chosen as the primary segmentation model for several reasons:

1. **Training data relevance**: ConSeP was trained on the Colorectal Adenocarcinoma (CoNSeP) dataset, making it directly applicable to this CRC cohort
2. **Tissue-specific classes**: ConSeP distinguishes malignant from healthy epithelial cells and identifies muscle tissue, which are relevant distinctions in colorectal tissue
3. **Validation**: Panel B confirms that despite using different classification schemas, ConSeP and CoNIC models show strong agreement on biologically meaningful cell categories

The tradeoff is that ConSeP groups all immune cells into a single "Inflammatory" class, whereas CoNIC provides finer granularity (Lymphocyte, Plasma cell, Neutrophil, Eosinophil). For studies requiring immune cell subtyping, the CoNIC model results are available in the supplementary materials.

## Files

- `figure1.png` - High-resolution PNG (300 DPI, 16 × 14 inches)
- `figure1.pdf` - Vector PDF for publication
- `figure1_caption.md` - This caption file
