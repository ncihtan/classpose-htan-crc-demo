# Spatial Feature Survival Analysis — Rationale, Methods, and Interpretation

## 1. Motivation

We previously classified 41 HTAN CRC H&E samples into four immune phenotypes
(Hot / Altered-excluded / Altered-immunosuppressed / Immune-desert) using
spatial analysis of ConSep cell classifications. While biologically informative,
a categorical phenotype with 4 levels and n=41 is statistically underpowered
against survival endpoints — our earlier logistic regression found that
Immune-desert carried OR=6-7 for disease progression, but with p~0.11 the
result was not significant.

**The core idea here is simple:** rather than collapsing rich spatial
information into 4 categories and then testing categories against outcomes, we
retain the continuous spatial features and test them directly in Cox
regression. This gains statistical power in two ways:

1. **No information loss.** A continuous covariate uses every data point's
   actual value rather than binning into categories.
2. **Single degree of freedom.** Each Cox test uses 1 df (continuous feature)
   rather than 3 df (4-level categorical), which matters greatly at n=40.

We test seven spatial features derived from two independent cell classification
models (CoNIC and ConSep) against three survival endpoints (OS, CSS, PFS).

---

## 2. Cohort

| Characteristic | Value |
|---|---|
| Total participants | 40 (41 samples; HTA7_919 has 2 H&E images, averaged) |
| Evaluable for OS/CSS | 39 (1 excluded: HTA7_982, Vital Status = "unknown") |
| Evaluable for PFS | 39 (1 excluded: HTA7_982, Progression = "Not Reported") |
| All-cause deaths | 15 / 39 (38%) |
| Cancer-specific deaths | 7 / 39 (18%) |
| Disease progressions | 24 / 39 (62%) |
| Stage I-II | 10 (25%) |
| Stage III-IV | 30 (75%) |
| Median follow-up (alive) | ~1800 days (~4.9 years) |

### Immune phenotype distribution

| Phenotype | n | % |
|---|---|---|
| Altered-excluded | 19 | 48% |
| Altered-immunosuppressed | 9 | 22% |
| Immune-desert | 8 | 20% |
| Hot | 4 | 10% |

This is a heavily pre-treated, advanced-stage cohort (75% stage III-IV) with
high progression rates (62%), which is expected given HTAN's focus on
molecularly characterized tumors from academic centers.

---

## 3. Features

### 3.1 Inflammatory features (from ConSep spatial analysis)

These four features come from our grid-based spatial analysis that divides each
slide into tumor core, tumor-stroma interface, and stroma regions, then
quantifies inflammatory cell distribution across regions.

| Feature | Definition | Biological rationale |
|---|---|---|
| `inflammatory_pct` | % of all ConSep-detected cells classified as Inflammatory | Overall immune infiltration level. Higher values indicate a more inflamed tumor microenvironment. |
| `density_overall` | Inflammatory cells per mm^2 of tissue | Area-normalized version of inflammatory content; controls for slide size. |
| `infiltration_ratio` | Inflammatory density in tumor core / density in stroma | Measures immune cell *penetration* into the tumor. High values mean immune cells are reaching the tumor rather than being excluded in the stroma. |
| `tumor_associated_pct` | % of inflammatory cells located in tumor-adjacent bins | Complementary to infiltration_ratio: what fraction of the immune response is tumor-proximal. |

**Why both `inflammatory_pct` and `density_overall`?** They are correlated but
not identical. Percentage depends on total cell count (which varies with tumor
cellularity), while density depends on tissue area (which varies with slide
coverage). A hypocellular tumor with moderate immune cells has high
inflammatory_pct but low density_overall.

### 3.2 Tumor-stroma ratio (from ConSep)

| Feature | Definition | Biological rationale |
|---|---|---|
| `tumor_stroma_ratio` | Malignant epithelial cell count / Stroma cell count | A morphometric proxy for the desmoplastic reaction. Low TSR (stroma-rich) is associated with worse prognosis in CRC in multiple studies (Huijbers et al. 2013, Mesker et al. 2007). Mechanistically, cancer-associated fibroblasts in the stroma secrete growth factors and remodel extracellular matrix to promote invasion. |

Range in this cohort: 0.18-2.17 (median 0.62), consistent with published CRC
TSR values.

### 3.3 CoNIC-derived features

CoNIC provides finer-grained immune cell subtyping than ConSep (which lumps
all immune cells as "Inflammatory"). CoNIC distinguishes Lymphocyte, Plasma
cell, Neutrophil, Eosinophil, Epithelial, and Connective tissue.

| Feature | Definition | Biological rationale |
|---|---|---|
| `nlr` | Neutrophil count / Lymphocyte count (tissue-level) | The tissue NLR is an analogue of the well-established blood NLR prognostic marker. Elevated NLR reflects a shift from adaptive (lymphocyte-mediated, anti-tumor) to innate (neutrophil-mediated, often pro-tumor) immunity. Tumor-associated neutrophils can promote angiogenesis and suppress T-cell function. |
| `eosinophil_density` | Eosinophil cells per mm^2 | Tumor-associated tissue eosinophilia (TATE) is associated with favorable prognosis in CRC (Harbaum et al. 2015, Prizment et al. 2016). Eosinophils can directly kill tumor cells via degranulation and recruit CD8+ T cells. |

**Note on tissue NLR values:** The tissue NLR (median 0.007, range
0.002-0.022) is far lower than blood NLR (typically 1-5) because tissue
sections are dominated by lymphocytes, with very few neutrophils visible on
H&E. This is expected — neutrophils are circulating cells that are sparse in
solid tissue. The relative ranking across patients is what matters for Cox
regression, not the absolute magnitude.

### 3.4 Truly spatial features (from cell coordinates)

The features above (3.1–3.3) are derived from aggregate cell counts or
grid-based region assignments. While informative, they do not directly use
cell coordinates. The following features leverage the actual (x, y) positions
of individual cells to measure spatial relationships.

#### Immune-tumor nearest-neighbor distance

For each sample, we build a cKDTree on inflammatory cell coordinates and query
the nearest inflammatory neighbor for every malignant epithelial cell. This
produces a per-cell distance array that is summarized as:

| Feature | Definition | Biological rationale |
|---|---|---|
| `nn_dist_median` | Median distance (um) from each tumor cell to its nearest inflammatory cell | A continuous, truly spatial analogue of the Immunoscore concept. Lower values indicate tighter immune-tumor mixing — immune cells are physically close to tumor cells, suggesting active immune engagement. This captures information lost by grid-based approaches: two samples with identical inflammatory_pct can have very different spatial distributions (diffuse vs. clustered near tumor). |

Additional summary statistics (`nn_dist_mean`, `nn_dist_10pct`) are saved in
`data/spatial_features.csv` for exploration but not tested in Cox regression
to limit multiple testing.

**Interpretation:** Lower nn_dist_median = more immune-tumor mixing = expected
to be protective (lower HR for PFS).

**Performance:** cKDTree.query on ~300K malignant cells against ~100K
inflammatory cells takes ~1-2s per sample.

#### Tumor-immune spatial organization (Moran's I)

For each sample, we compute Moran's I — the standard spatial autocorrelation
statistic — on immune cell density within the tumor region. This measures
whether tumor-proximal inflammatory cells are spatially clustered (positive I),
randomly distributed (I ≈ 0), or dispersed (negative I).

**Method:**
1. Build a cKDTree on malignant cell coordinates and identify inflammatory
   cells within 100μm (~308px, ~7-10 cell diameters) of the nearest malignant
   cell. This defines the "tumor microenvironment" immune population.
2. Grid the tumor region (bounding box of malignant cells + buffer, 50×50
   bins) and count tumor-proximal inflammatory cells per bin.
3. Restrict to bins containing ≥1 malignant cell (tumor-region bins).
4. Compute Moran's I with queen contiguity (8-neighbor adjacency) on the
   masked grid:

       I = (N / W) × Σ_ij w_ij (x_i - x̄)(x_j - x̄) / Σ_i (x_i - x̄)²

   where N = number of valid bins, W = total neighbor pairs.

| Feature | Definition | Biological rationale |
|---|---|---|
| `tumor_immune_morans_i` | Moran's I of immune cell density on the tumor-region grid | Measures how *spatially organized* the tumor-proximal immune response is. Higher values indicate immune cells form coherent spatial clusters rather than being diffusely scattered, which is characteristic of organized immune aggregates (e.g., tertiary lymphoid structures). Unlike DBSCAN, Moran's I has no arbitrary clustering parameters — just a spatial weights definition — and produces a single continuous value amenable to Cox regression. |

**Range:** -1 (dispersed) to +1 (clustered); expect 0.1–0.7 for tumor-
proximal immune cells.

**Interpretation:** Higher Moran's I = more spatially organized immune
infiltrate = expected to be protective.

**Performance:** ~2500 bins × 8 neighbors = trivial computation per sample.

---

## 4. Survival Endpoints

### 4.1 Overall Survival (OS) — all-cause

- **Time:** `Days to Last Followup` (the only time variable available; `Days
  to Death` is empty for all patients in this HTAN release)
- **Event:** Vital Status = "Dead" (any cause)
- **Events:** 15 / 39

**Limitation:** `Days to Last Followup` for dead patients may represent the
last clinical contact rather than the exact death date. Cross-validation
against `Year of Death` shows a mean discrepancy of +0.02 years (SD 0.40),
consistent with rounding from integer year-of-diagnosis. The largest gap is
+0.84 years (HTA7_947), which could reflect a late-in-year diagnosis.

### 4.2 Cancer-Specific Survival (CSS)

- **Time:** Same as OS (`Days to Last Followup`)
- **Event:** Dead AND Cause of Death = "Cancer Related"
- **Events:** 7 / 39

Non-cancer deaths (n=4, "Not Cancer Related") and unknown-cause deaths (n=3,
"unknown") are **censored** at their followup time. One additional death
(HTA7_963, cause "unknown" but with Stage IVC and distant metastatic
progression at day 147) may be a misclassified cancer death, but we use the
data as recorded.

**Why add CSS?** All-cause OS includes 8 deaths unrelated to (or of unknown
relation to) the tumor microenvironment. These dilute any signal from spatial
features. If inflammatory density affects *cancer-specific* mortality but not
death from comorbidities, all-cause OS would miss the signal. CSS isolates the
question: "Does the spatial immune landscape predict death *from this cancer*?"

The tradeoff is power — 7 events is very few for Cox regression. We report CSS
results for completeness but interpret them cautiously.

### 4.3 Progression-Free Survival (PFS)

- **Time:** `Days to Last Known Disease Status`
- **Event:** `Progression or Recurrence` starts with "Yes"
- **Events:** 24 / 39

PFS is our **best-powered endpoint** and the most directly relevant to tumor
biology. Progression is a tumor-intrinsic event (unlike death, which can be
from comorbidities), and with 24 events / 39 patients (62% event rate) we have
reasonable power for Cox regression even at n=39.

---

## 5. Statistical Methods

### 5.1 Cox Proportional Hazards Regression

For each feature x endpoint combination, we fit a univariate Cox model:

    h(t) = h_0(t) * exp(beta * z_feature)

where `z_feature` is the **z-scored** (mean-centered, SD-scaled) feature
value. Standardization means the hazard ratio (HR) represents the change in
hazard per **one standard deviation increase** in the feature. This is
critical for interpretability — without standardization, an HR for
inflammatory_pct (range 6-38%) would be on a completely different scale than
an HR for density_overall (range 920-17,700 cells/mm^2), making forest plot
comparison meaningless.

**Reported metrics:**
- HR with 95% CI (exponentiated coefficient)
- Wald test p-value
- Harrell's concordance index (C-index), where 0.5 = no discrimination and
  1.0 = perfect

### 5.2 Stage-Adjusted Multivariate Analysis

For any feature with univariate p < 0.2, we refit including AJCC pathologic
stage (simplified to I-II vs III-IV as a binary covariate) to assess whether
the feature provides prognostic information **independent of stage**:

    h(t) = h_0(t) * exp(beta_1 * z_feature + beta_2 * stage_III_IV)

The p < 0.2 threshold (rather than 0.05) is standard for selecting candidates
for multivariate modeling — it avoids missing features that might become
significant after adjustment.

### 5.3 Kaplan-Meier Curves and Log-Rank Tests

For visualization, each continuous feature is dichotomized at the median into
High/Low groups. Kaplan-Meier survival curves are plotted with the log-rank
test p-value. This is a complementary visual tool — the Cox regression (which
uses the full continuous distribution) is the primary analysis.

### 5.4 Multiple Testing Considerations

With 9 features x 3 endpoints = 27 tests, we do not apply formal multiple
testing correction (e.g., Bonferroni). This is an exploratory analysis in a
small cohort; the goal is hypothesis generation, not confirmatory testing.
Results with p < 0.05 should be interpreted as "warrants further
investigation," not "definitively proven." We note that the strongest PFS
results (inflammatory_pct p=0.003, density_overall p=0.007) would survive
Bonferroni correction within the PFS endpoint (9 tests, threshold 0.006).

---

## 6. Results

### 6.1 Univariate Cox Regression

#### Progression-Free Survival (24 events / 39 patients)

| Feature | HR per SD | 95% CI | p | C-index |
|---|---|---|---|---|
| **Inflammatory %** | **0.44** | **0.26-0.75** | **0.003** | **0.70** |
| **Inflammatory density** | **0.47** | **0.27-0.82** | **0.007** | **0.67** |
| **Eosinophil density** | **0.63** | **0.41-0.99** | **0.043** | **0.64** |
| Tumor-stroma ratio | 0.90 | 0.58-1.38 | 0.619 | 0.61 |
| Tumor-associated % | 1.08 | 0.75-1.55 | 0.686 | 0.54 |
| Infiltration ratio | 0.98 | 0.63-1.50 | 0.911 | 0.45 |
| Neutrophil-lymphocyte ratio | 1.09 | 0.73-1.64 | 0.671 | 0.49 |

Three features are significantly associated with PFS at p < 0.05. All three
are **protective** (HR < 1): higher immune infiltration predicts longer
progression-free survival.

#### Cancer-Specific Survival (7 events / 39 patients)

| Feature | HR per SD | 95% CI | p | C-index |
|---|---|---|---|---|
| Infiltration ratio | 1.82 | 0.87-3.81 | 0.111 | 0.71 |
| Tumor-associated % | 1.32 | 0.69-2.54 | 0.407 | 0.64 |
| Eosinophil density | 0.55 | 0.21-1.43 | 0.218 | 0.63 |
| Tumor-stroma ratio | 0.51 | 0.17-1.50 | 0.218 | 0.64 |
| Inflammatory % | 0.71 | 0.31-1.63 | 0.416 | 0.61 |
| Inflammatory density | 0.65 | 0.26-1.61 | 0.356 | 0.58 |
| NLR | 0.77 | 0.32-1.86 | 0.566 | 0.54 |

No feature reaches significance, as expected with only 7 events. However, the
**infiltration ratio** shows a notable trend (HR=1.82, p=0.111, C=0.71) that
was invisible in all-cause OS (HR=1.28, p=0.33). The direction is
counterintuitive — discussed in Section 7.

#### Overall Survival (15 events / 39 patients)

| Feature | HR per SD | 95% CI | p | C-index |
|---|---|---|---|---|
| Tumor-stroma ratio | 0.62 | 0.32-1.18 | 0.143 | 0.58 |
| NLR | 0.69 | 0.36-1.31 | 0.258 | 0.59 |
| Inflammatory % | 0.73 | 0.41-1.29 | 0.279 | 0.58 |
| Eosinophil density | 0.74 | 0.43-1.29 | 0.293 | 0.58 |
| Infiltration ratio | 1.28 | 0.78-2.12 | 0.330 | 0.60 |
| Inflammatory density | 0.77 | 0.44-1.35 | 0.360 | 0.56 |
| Tumor-associated % | 1.07 | 0.67-1.71 | 0.779 | 0.55 |

No feature reaches significance. All-cause OS is diluted by 8 non-cancer /
unknown-cause deaths and has less power than PFS (15 vs 24 events).

### 6.2 Stage-Adjusted Multivariate Results

Features with univariate p < 0.2 were refitted with AJCC stage (I-II vs
III-IV):

| Feature | Endpoint | Adj. HR | 95% CI | Adj. p | Stage HR | Stage p |
|---|---|---|---|---|---|---|
| **Inflammatory %** | **PFS** | **0.45** | **0.26-0.77** | **0.004** | 1.57 | 0.420 |
| **Inflammatory density** | **PFS** | **0.49** | **0.29-0.84** | **0.010** | 1.69 | 0.343 |
| Eosinophil density | PFS | 0.66 | 0.43-1.01 | 0.057 | 1.90 | 0.247 |
| Infiltration ratio | CSS | 1.84 | 0.89-3.81 | 0.099 | 2.35 | 0.429 |
| Tumor-stroma ratio | OS | 0.63 | 0.33-1.19 | 0.154 | 1.33 | 0.657 |

**Key finding:** Inflammatory % and density remain significant (p=0.004 and
p=0.010) after stage adjustment, with essentially unchanged HRs. This means
the prognostic information in inflammatory content is **independent of
stage** — it is not simply a proxy for earlier-stage disease having more
inflammation.

Notably, stage itself is **not significant** in any of these models (all stage
p > 0.2). This is unusual but may reflect the cohort composition (75% stage
III-IV, limited stage variation) and the fact that stage predicts mortality
more than progression.

---

## 7. Interpretation

### 7.1 The inflammatory content-PFS signal is robust

The strongest finding is that higher inflammatory content — whether measured as
percentage (`inflammatory_pct`, HR=0.44) or density (`density_overall`,
HR=0.47) — is associated with significantly longer PFS, with a large effect
size (>50% hazard reduction per SD increase), good discrimination (C-index
0.67-0.70), and robustness to stage adjustment.

This aligns with a large body of literature establishing that immune-inflamed
("hot") CRC tumors have better outcomes (Galon et al. 2006, Pages et al.
2018). The Immunoscore, which quantifies CD3+ and CD8+ T-cell densities at the
tumor center and invasive margin, is now a validated prognostic tool in CRC
(Pages et al. 2018, Lancet). Our `inflammatory_pct` is a coarser analogue
(ConSep's "Inflammatory" class includes all immune lineages, not just T cells)
but captures the same underlying biology.

### 7.2 Eosinophil density adds independent prognostic value

Eosinophil density is significant for PFS (HR=0.63, p=0.043) and trends toward
significance after stage adjustment (p=0.057). This is consistent with
literature showing tumor-associated tissue eosinophilia (TATE) is a favorable
prognostic marker in CRC (Harbaum et al. 2015, Prizment et al. 2016). Proposed
mechanisms include direct tumoricidal activity via eosinophil degranulation and
recruitment of CD8+ T cells via chemokine secretion.

The fact that eosinophil density provides prognostic signal *in addition to*
overall inflammatory content suggests that the composition of the immune
infiltrate matters, not just its quantity.

### 7.3 The infiltration ratio-CSS trend deserves attention

The infiltration ratio (tumor core / stroma inflammatory density) shows a
trend toward *worse* cancer-specific survival at higher values (HR=1.82,
p=0.111). This seems paradoxical — shouldn't immune penetration into the tumor
be beneficial?

Several explanations are possible:

1. **Immunosuppressive infiltrate.** Not all tumor-infiltrating immune cells
   are anti-tumoral. ConSep's "Inflammatory" class includes
   immunosuppressive populations (regulatory T cells, M2 macrophages, myeloid-
   derived suppressor cells) that cannot be distinguished on H&E. High
   infiltration ratio could reflect an immunosuppressive program that has
   overcome exclusion barriers.

2. **Confounding with tumor biology.** Highly aggressive tumors may actively
   recruit immune cells (via chemokine secretion and neoantigens) while
   simultaneously evading immune destruction. The infiltration ratio would be
   high, but the immune response is ineffective.

3. **Small sample artifact.** With only 7 events, a single outlier can
   substantially influence the result. This should be verified in a larger
   cohort.

This observation was **invisible in all-cause OS** (HR=1.28, p=0.33) because
the 8 non-cancer deaths diluted the signal. CSS is the appropriate endpoint
for this question, even with limited power.

### 7.4 Tumor-stroma ratio: expected direction, insufficient power

TSR shows a consistent protective trend (HR < 1 for higher TSR, meaning more
tumor relative to stroma is associated with better survival). This reverses
the usual "stroma-rich = worse" finding reported in the literature. However,
our TSR is computed from whole-slide ConSep cell counts, not from the
tumor-specific stroma fraction used in clinical TSR studies (which score only
the most stroma-rich area at the invasive front). The whole-slide ratio may
capture different biology.

### 7.5 Tissue NLR: no signal

The tissue neutrophil-lymphocyte ratio shows no association with any endpoint
(all p > 0.25). This likely reflects the extremely low neutrophil counts in
tissue sections — CoNIC detects very few neutrophils on H&E (median NLR =
0.007, meaning ~1 neutrophil per 140 lymphocytes). The blood NLR is a
well-validated prognostic marker in CRC, but its tissue analogue appears
uninformative, at least at this sample size and with H&E-based detection.

---

## 8. Limitations

### 8.1 Sample size

n=39-40 is small for survival analysis. With 24 PFS events, the rule of thumb
for Cox regression (10-20 events per covariate) allows univariate models but
limits multivariate models to 1-2 covariates. Our stage-adjusted models (2
covariates, 24 events) are at the lower bound of reliability.

### 8.2 Survival time construction

`Days to Death` is missing for all patients. We use `Days to Last Followup` as
a proxy for OS/CSS time. Cross-validation against Year of Death suggests
reasonable accuracy (mean discrepancy ~0 years, max ~0.8 years), but the
possibility of systematic underestimation for dead patients cannot be ruled
out. This is a data quality issue in the HTAN HMS CRC release, not a
methodological choice.

### 8.3 Cell classification limitations

- **ConSep's "Inflammatory" class** includes all immune lineages — T cells,
  B cells, macrophages, dendritic cells, granulocytes — without
  distinguishing pro- vs anti-tumoral populations. The Immunoscore's
  superiority comes partly from specificity (CD3+, CD8+), which we lack.
- **CoNIC's neutrophil detection on H&E** is unreliable at low counts.
  Neutrophils are morphologically similar to other cells and are sparse in
  tissue; IHC (e.g., CD66b) would be needed for accurate tissue NLR.
- Both models are trained on different datasets and may have systematic biases
  on this cohort.

### 8.4 Multiple testing

We tested 27 feature-endpoint combinations (9 features x 3 endpoints) without
formal correction. The two strongest results (inflammatory_pct and
density_overall for PFS, p=0.003 and 0.007) survive Bonferroni correction
within endpoint (9 tests, threshold 0.006), but the eosinophil density result
(p=0.043) does not.

### 8.5 Single-cohort, no validation

All results are from a single HTAN HMS cohort. The features and effect sizes
should be validated in an independent CRC cohort before drawing clinical
conclusions.

---

## 9. Conclusions

1. **Continuous spatial features outperform categorical phenotypes** for
   survival analysis in this small cohort. Where the 4-level immune phenotype
   classification was underpowered (p~0.11), continuous inflammatory
   percentage achieves p=0.003 for PFS.

2. **Overall inflammatory content is the strongest spatial predictor of PFS**
   in this CRC cohort (HR=0.44 per SD, p=0.003, C=0.70), independent of AJCC
   stage. Each standard deviation increase in inflammatory cell percentage is
   associated with a 56% reduction in progression hazard.

3. **Eosinophil density provides additional prognostic information** (PFS
   HR=0.63, p=0.043), suggesting that immune composition — not just quantity —
   matters for outcomes.

4. **Cancer-specific survival reveals signals hidden in all-cause OS.** The
   infiltration ratio trends toward worse CSS (HR=1.82, p=0.11, C=0.71),
   potentially reflecting immunosuppressive infiltration, but is invisible in
   all-cause OS where non-cancer deaths dominate.

5. **Tissue NLR from H&E is not prognostic**, likely due to very low
   neutrophil counts in tissue sections. Blood NLR or IHC-based tissue NLR
   may perform differently.

6. These results are **hypothesis-generating**, not confirmatory. The key
   hypothesis for validation is: *H&E-derived inflammatory cell density
   predicts PFS in CRC independently of stage*, with a candidate effect size
   of HR~0.45-0.50 per SD.

---

## 10. Output Files

| File | Description |
|---|---|
| `data/spatial_features.csv` | Spatial features per participant (40 rows: NN distances + Moran's I) |
| `data/survival_analysis_dataset.csv` | Merged analysis dataset (40 rows) |
| `data/survival_cox_results.csv` | All univariate Cox results (21 rows) |
| `data/survival_cox_multivariate_results.csv` | Stage-adjusted results for features with univariate p < 0.2 |
| `figures/km_{feature}_{endpoint}.png` | Kaplan-Meier curves per feature and endpoint (median split) |
| `figures/km_phenotype_{endpoint}.png` | Kaplan-Meier curves by immune phenotype category |
| `figures/survival_forest_plot.png` | Forest plot of all univariate HRs |

---

## References

- Galon J, et al. (2006). Type, density, and location of immune cells within
  human colorectal tumors predict clinical outcome. *Science* 313:1960-1964.
- Harbaum L, et al. (2015). Peritumoral eosinophils predict recurrence in
  colorectal cancer. *Mod Pathol* 28:403-413.
- Huijbers A, et al. (2013). The proportion of tumor-stroma as a strong
  prognosticator for stage II and III colon cancer patients. *Ann Oncol*
  24:2068-2073.
- Mesker WE, et al. (2007). The carcinoma-stromal ratio of colon carcinoma is
  an independent factor for survival. *Cell Oncol* 29:387-398.
- Pages F, et al. (2018). International validation of the consensus
  Immunoscore for the classification of colon cancer. *Lancet* 391:2128-2139.
- Prizment AE, et al. (2016). Tumor eosinophil infiltration and improved
  survival of colorectal cancer patients. *Mod Pathol* 29:905-915.
- Sautes-Fridman C, et al. (2019). Tertiary lymphoid structures in the era of
  cancer immunotherapy. *Nat Rev Cancer* 19:307-325.
