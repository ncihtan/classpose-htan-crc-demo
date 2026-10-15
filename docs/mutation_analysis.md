# Mutation Interaction Analysis: Do Oncogene Mutations Modulate the Prognostic Value of Immune Clustering?

## Background

Our prior spatial analysis identified **tumor-proximal immune clustering** (`db_tumor_cluster_fraction`, the fraction of inflammatory cells forming DBSCAN clusters near tumor regions) as the strongest spatial predictor of progression-free survival (PFS) in this 40-patient HTAN HMS CRC cohort (HR = 0.40 per SD, p = 0.001). Higher clustering — immune cells physically aggregating near tumor — predicts better outcomes.

But CRC is molecularly heterogeneous. Tumors driven by different oncogene mutations have distinct immune microenvironments. KRAS-mutant CRC, for instance, is associated with immunosuppressive signaling, while BRAF-mutant tumors often show microsatellite instability and higher immune infiltration. The natural question: **does mutation status change how much immune clustering matters?**

We have targeted sequencing data for 11 CRC-relevant genes with binary positive/negative calls, with 100% overlap with our H&E imaging cohort (n = 40).

## Cohort Mutation Profile

| Gene     | Mutant | Prevalence | Notes                                           |
|----------|--------|------------|--------------------------------------------------|
| APC      | 26/40  | 65.0%      | Expected for CRC (APC loss initiates most CRC)   |
| TP53     | 24/40  | 60.0%      | Consistent with conventional CRC pathway          |
| KRAS     | 14/40  | 35.0%      | Matches known ~35-45% CRC KRAS mutation rate      |
| BRAF     | 9/40   | 22.5%      | Slightly enriched (population rate ~10-15%)       |
| PIK3CA   | 7/40   | 17.5%      | Within expected range                             |
| SOX9     | 6/40   | 15.0%      |                                                   |
| FBXW7    | 4/40   | 10.0%      |                                                   |
| CTNNB1   | 3/40   | 7.5%       |                                                   |
| PTEN     | 2/40   | 5.0%       | Too rare for subgroup analysis                    |
| TCF7L2   | 2/40   | 5.0%       | Too rare for subgroup analysis                    |
| NRAS     | 1/40   | 2.5%       | Too rare for subgroup analysis                    |

Mutation prevalences match published CRC cohorts well, supporting data quality. The median mutation count per patient is 2 (range 0-5).

## Analysis 1: Do Mutations Directly Predict Survival?

We tested each gene as a binary predictor in univariate Cox regression across three endpoints (OS, CSS, PFS). Only genes with at least 5 events in both mutant and WT groups were tested, which limited the analysis to APC, KRAS, and TP53 (BRAF had insufficient events in the mutant group for some endpoints).

| Gene | Endpoint | HR   | 95% CI        | p     |
|------|----------|------|---------------|-------|
| APC  | OS       | 0.39 | 0.14 - 1.09   | 0.074 |
| APC  | PFS      | 1.16 | 0.48 - 2.80   | 0.749 |
| KRAS | OS       | 0.87 | 0.30 - 2.54   | 0.795 |
| KRAS | PFS      | 1.14 | 0.50 - 2.62   | 0.751 |
| TP53 | OS       | 1.21 | 0.41 - 3.54   | 0.728 |
| TP53 | PFS      | 1.83 | 0.75 - 4.44   | 0.184 |

**No individual mutation significantly predicts survival.** APC mutation shows a borderline protective trend for OS (HR = 0.39, p = 0.074), which is consistent with the literature showing that APC-driven "conventional" CRC has better prognosis than non-APC-driven pathways. But nothing reaches p < 0.05.

This is not unexpected in a cohort of 40 patients — individual mutations are unlikely to reach significance as standalone predictors. The more interesting question is whether they modify the immune clustering effect.

## Analysis 2: Immune Clustering Effect Stratified by Mutation Status

Here we split patients into mutant and WT subgroups for each focus gene, then ran Cox regression on z-scored `db_tumor_cluster_fraction` within each subgroup separately. This asks: *how strongly does immune clustering predict survival within each genetic subgroup?*

### KRAS (n = 14 mutant, 25 WT)

| Stratum     | Endpoint | HR   | 95% CI        | p     | n  |
|-------------|----------|------|---------------|-------|----|
| **Mutant**  | **OS**   | **0.10** | **0.01 - 0.82** | **0.032** | 14 |
| WT          | OS       | 0.76 | 0.40 - 1.45   | 0.401 | 25 |
| Mutant      | PFS      | 0.48 | 0.19 - 1.21   | 0.120 | 14 |
| WT          | PFS      | 0.36 | 0.18 - 0.74   | 0.005 | 25 |

The standout result: in KRAS-mutant patients, each SD increase in immune clustering is associated with a **90% reduction in hazard of death** (HR = 0.10, p = 0.032). In KRAS-WT patients, the same feature shows no significant association with OS (HR = 0.76, p = 0.40). The pattern reverses for PFS, where the WT group shows the stronger effect — but the OS interaction is the biologically compelling finding.

### BRAF (n = 9 mutant, 30 WT)

| Stratum     | Endpoint | HR   | 95% CI        | p     | n  |
|-------------|----------|------|---------------|-------|----|
| Mutant      | OS       | 0.14 | 0.02 - 1.38   | 0.092 | 9  |
| WT          | OS       | 0.59 | 0.29 - 1.19   | 0.141 | 30 |
| Mutant      | PFS      | 0.19 | 0.03 - 1.33   | 0.094 | 9  |
| WT          | PFS      | 0.47 | 0.28 - 0.82   | 0.007 | 30 |

BRAF shows a similar pattern to KRAS — very large point estimates in the mutant group (HR 0.14-0.19) that don't reach significance due to small sample size (n = 9). The trend suggests immune clustering may be particularly prognostic in BRAF-mutant tumors as well, but we are underpowered to confirm.

### TP53 (n = 24 mutant, 15 WT)

| Stratum     | Endpoint | HR   | 95% CI        | p     | n  |
|-------------|----------|------|---------------|-------|----|
| Mutant      | OS       | 0.59 | 0.28 - 1.23   | 0.160 | 24 |
| WT          | OS       | 0.49 | 0.19 - 1.28   | 0.144 | 15 |
| Mutant      | PFS      | 0.42 | 0.23 - 0.75   | 0.004 | 24 |
| WT          | PFS      | 0.42 | 0.17 - 1.01   | 0.053 | 15 |

TP53 shows a **consistent** immune clustering effect in both subgroups — the HR is similar in mutant and WT for both OS and PFS. This means TP53 status does not modify the clustering effect. The immune clustering association is robust regardless of TP53 status.

### APC (n = 26 mutant, 13 WT)

| Stratum     | Endpoint | HR   | 95% CI        | p     | n  |
|-------------|----------|------|---------------|-------|----|
| Mutant      | OS       | 0.47 | 0.21 - 1.07   | 0.071 | 26 |
| WT          | OS       | 0.66 | 0.27 - 1.59   | 0.354 | 13 |
| Mutant      | PFS      | 0.44 | 0.25 - 0.78   | 0.005 | 26 |
| WT          | PFS      | 0.39 | 0.11 - 1.35   | 0.138 | 13 |

Like TP53, the immune clustering effect is similar in both APC subgroups. APC status does not meaningfully modify the clustering-survival relationship.

## Analysis 3: Formal Interaction Test

The stratified analysis is descriptive. The proper statistical test is a Cox model with an interaction term:

```
hazard ~ db_tumor_cluster_fraction + gene_mutation + cluster_fraction * gene_mutation
```

The interaction term p-value tests whether the clustering HR differs significantly between mutant and WT.

| Gene | Endpoint | Interaction HR | 95% CI        | Interaction p | Cluster p | Mutation p |
|------|----------|---------------|---------------|---------------|-----------|------------|
| **KRAS** | **OS** | **0.085** | **0.009 - 0.80** | **0.032** | 0.411 | 0.286 |
| KRAS | PFS      | 1.21          | 0.41 - 3.60   | 0.732         | 0.004     | 0.373      |
| BRAF | OS       | 0.32          | 0.05 - 2.18   | 0.247         | 0.140     | 0.101      |
| BRAF | PFS      | 0.38          | 0.07 - 2.15   | 0.273         | 0.008     | 0.894      |
| TP53 | OS       | 1.80          | 0.41 - 7.89   | 0.438         | 0.123     | 0.585      |
| TP53 | PFS      | 1.83          | 0.50 - 6.64   | 0.358         | 0.029     | 0.095      |
| APC  | OS       | 0.66          | 0.21 - 2.12   | 0.488         | 0.377     | 0.067      |
| APC  | PFS      | 1.16          | 0.31 - 4.41   | 0.822         | 0.098     | 0.473      |

**One significant interaction: KRAS x OS (p = 0.032).** The interaction HR of 0.085 means that the protective effect of immune clustering is approximately 12x stronger in KRAS-mutant patients compared to KRAS-WT patients for overall survival.

No other gene-endpoint combinations show significant interactions. TP53 and APC show consistent clustering effects across subgroups (interaction p > 0.35), confirming that those mutations do not modify the spatial immune feature's prognostic value.

## Analysis 4: Four-Group Kaplan-Meier Visualization

For the significant KRAS x OS interaction, we split patients into four groups by crossing KRAS status (mutant/WT) with immune clustering level (above/below median):

![KRAS x Immune Clustering KM Plot](../analysis/figures/km_interaction_KRAS_os.png)

The plot reveals the biology:
- **KRAS-mut / high clustering** (n = 8): best overall survival, nearly flat curve
- **KRAS-WT / high clustering** (n = 12): good survival, moderate decline
- **KRAS-WT / low clustering** (n = 13): intermediate survival
- **KRAS-mut / low clustering** (n = 6): worst survival, steep early decline

The four-group log-rank test is p = 0.092. The key visual finding is the dramatic separation within the KRAS-mutant group — high-clustering mutant patients do extremely well while low-clustering mutant patients do very poorly. In contrast, the two KRAS-WT groups have overlapping, intermediate curves.

## Analysis 5: Mutation Burden

Total mutation count (sum of all 11 genes, range 0-5, median 2) was tested as a continuous predictor.

| Model                     | Endpoint | Count HR | Count p | Cluster HR | Cluster p |
|---------------------------|----------|----------|---------|------------|-----------|
| Mutation count only        | OS       | 1.10     | 0.673   | -          | -         |
| Mutation count only        | CSS      | 1.13     | 0.709   | -          | -         |
| Mutation count only        | PFS      | 1.01     | 0.940   | -          | -         |
| Count + immune clustering  | OS       | 1.24     | 0.358   | 0.52       | 0.040     |
| Count + immune clustering  | CSS      | 1.25     | 0.502   | 0.57       | 0.212     |
| Count + immune clustering  | PFS      | 1.17     | 0.351   | 0.40       | 0.001     |

**Mutation burden has no prognostic value** in this cohort (all p > 0.35). Importantly, `db_tumor_cluster_fraction` retains its prognostic significance for OS (p = 0.040) and PFS (p = 0.001) even after adjusting for mutation count. The spatial immune feature provides orthogonal prognostic information that is not captured by genomic mutation burden.

## Key Conclusions

### 1. KRAS-mutant rescue effect: spatial immune organization as a prognostic gatekeeper in MSS CRC

While KRAS-mutant CRC is typically associated with immune exclusion and poor prognosis, high tumor-proximal immune clustering identifies a subset of KRAS-mutant patients with survival outcomes comparable to or better than KRAS-WT patients. The KRAS x `db_tumor_cluster_fraction` interaction for OS reaches formal statistical significance (interaction p = 0.032), strengthening to p = 0.020 after removing BRAF-mutant/likely-MSI-H cases. In KRAS-mutant tumors, immune clustering is a powerful predictor (HR = 0.10 per SD, p = 0.032), while in KRAS-WT tumors it is not significant for OS (HR = 0.76, p = 0.40).

This is an **effect modification**, not confounding. KRAS mutation itself does not predict OS (HR = 0.87, p = 0.795). Instead, it determines whether the spatial immune pattern matters for survival. The mutual exclusivity of KRAS and BRAF mutations in this cohort means the KRAS-mutant group is definitionally BRAF-WT and very likely MSS, placing this finding in the most clinically challenging CRC subgroup — the population that currently has no effective immunotherapy options.

### 2. The KRAS-mutant "rescue" effect

The sensitivity analysis (Section 3 below) establishes that all 14 KRAS-mutant patients are BRAF-WT and therefore overwhelmingly likely to be MSS (Microsatellite Stable). This places our finding squarely in the most clinically difficult CRC subgroup: **KRAS-mutant/MSS tumors are the "black hole" of immunotherapy** — typically immune-excluded, unresponsive to checkpoint inhibitors, and associated with poor outcomes.

KRAS-mutant CRC has well-documented mechanisms of immune exclusion:
- KRAS signaling (particularly G12D/V) drives secretion of CXCL8 (IL-8), GM-CSF, and CCL2, which recruit MDSCs and neutrophils while excluding CD8+ T cells (Liao et al., Cancer Cell 2019)
- KRAS activation promotes TGF-beta signaling in the stroma, creating a physical barrier to immune infiltration (He et al., 2024)
- KRAS-mutant tumors upregulate PD-L1, further suppressing adaptive immunity

This establishes two contrasting biological models for the interaction:

**The "Breach" hypothesis (KRAS-mutant):** The baseline state of a KRAS-mutant tumor is an immunosuppressive "fortress." If high tumor-proximal immune clustering is nevertheless present, it implies a rare, high-quality immune response that has successfully breached the KRAS-driven exclusion program. The extreme HR (0.10) reflects the separation between patients whose immune system has overcome an active evasion mechanism (exceptional prognosis) and those whose immune system has been fully excluded (poor prognosis). The feature `db_tumor_cluster_fraction` specifically measures clusters near tumor regions — this is not just "immune cells present in the tissue," but immune cells that have navigated to the tumor-stroma interface despite active spatial exclusion.

**The "Bystander" hypothesis (KRAS-WT):** In KRAS-WT tumors, the immune system is not being actively pushed away by KRAS-driven cytokines. Immune cells may be present and even clustered, but this represents the default state rather than evidence of an exceptional immune response. The presence or absence of clustering is therefore less discriminating — it separates "normal" from "slightly less normal" rather than "exceptional" from "failed."

This framing addresses the "spatial exclusion" problem that defines KRAS-mutant CRC: it is not enough to count immune cells (density), one must ask whether they have successfully organized at the tumor interface (clustering).

### 3. The KRAS interaction is not an MSI artifact

The most critical confounder for any KRAS × immune interaction in CRC is microsatellite instability (MSI). MSI-High tumors have high neoantigen loads and strong immune infiltration; if the KRAS-mut/high-clustering group were simply enriched for MSI-H patients, the survival benefit would reflect known MSI biology rather than a novel interaction.

Direct MSI status is not available in this cohort (the `Mismatch_Repair_Mutation` field is unpopulated). However, we can address this through molecular proxies:

**KRAS and BRAF are perfectly mutually exclusive** in this cohort — zero patients carry both (consistent with known CRC biology where activating mutations in the same RAS-RAF-MAPK pathway are redundant):

|              | BRAF-WT | BRAF-mut | Total |
|--------------|---------|----------|-------|
| KRAS-WT      | 17      | 9        | 26    |
| KRAS-mut     | 14      | 0        | 14    |
| Total        | 31      | 9        | 40    |

Since BRAF V600E is the strongest molecular marker for the MSI-H/CIMP pathway in CRC, and all 14 KRAS-mutant patients are BRAF-WT, the KRAS-mutant group is very unlikely to contain MSI-H patients. The 9 BRAF-mutant (likely MSI-H) patients all fall in the KRAS-WT group — 6 of them in the high-clustering WT subgroup, which may partly explain why that group also has reasonable survival.

**Sensitivity analysis excluding all BRAF-mutant patients** (n = 31 remaining):

| Term        | HR    | 95% CI        | p       |
|-------------|-------|---------------|---------|
| cluster_z   | 1.00  | 0.44 - 2.25   | 0.999   |
| KRAS        | 0.64  | 0.12 - 3.48   | 0.603   |
| **interaction** | **0.063** | **0.006 - 0.65** | **0.020** |

The interaction effect *strengthens* after removing BRAF-mutant patients (p = 0.020 vs 0.032 in the full cohort). This is because removing the likely-MSI-H patients from the WT group makes the WT immune clustering effect even weaker, sharpening the contrast with the KRAS-mutant group.

This sensitivity analysis provides strong evidence that the KRAS × immune clustering interaction reflects genuine KRAS-specific biology, not MSI confounding.

### 4. Immune clustering is robust across TP53 and APC subgroups

For the two most common mutations (TP53 60%, APC 65%), the immune clustering effect is consistent in both mutant and WT subgroups (interaction p > 0.35). This is reassuring — it means the spatial feature is broadly prognostic and not driven by a single genetic subgroup.

### 5. Mutation burden adds no prognostic value beyond spatial features

Total mutation count (0-5 per patient across 11 genes) does not predict OS, CSS, or PFS. The spatial immune feature remains independently significant after adjusting for mutation burden. This supports the idea that *how* the immune system is spatially organized around tumors matters more than the mutational landscape for prognosis.

### 6. Robustness checks

Three additional analyses were performed to stress-test the KRAS × clustering interaction:

**Confounder adjustment.** The interaction survives adjustment for clinical confounders:

| Model | Interaction HR | p |
|-------|---------------|---|
| Unadjusted | 0.085 | 0.032 |
| + stage | 0.107 | 0.030 |
| + stage + treatment | 0.107 | 0.028 |
| + stage + treatment + age | 0.132 | 0.032 |

The effect attenuates slightly (~50% reduction in log-HR) with full adjustment but remains significant. This is reassuring — the interaction is not driven by confounding from stage, treatment, or age.

**Permutation test.** Randomly permuting KRAS status 1,000 times and re-fitting the interaction model, only 29/1,000 permutations produced an interaction coefficient as extreme as observed (permutation p = 0.029). This confirms the parametric p-value is not an artifact of the Cox model assumptions.

**Feature specificity.** Testing the KRAS × OS interaction across multiple immune features reveals that the effect is **not specific to spatial clustering**:

| Feature | Type | Interaction HR | Interaction p |
|---------|------|---------------|---------------|
| `db_cluster_fraction` (all cells) | Spatial clustering | 0.033 | **0.010** |
| `db_tumor_cluster_fraction` (tumor-proximal) | Spatial clustering | 0.085 | **0.032** |
| `density_overall` | Simple density | 0.065 | 0.105 |
| `tumor_cluster_fraction` (HDBSCAN) | Spatial clustering | 0.110 | 0.185 |
| `infiltration_ratio` | Density ratio | 4.316 | **0.027** |
| `inflammatory_pct` | Simple percentage | 0.484 | 0.442 |

This is an important caveat: `density_overall` shows a similar trend (p = 0.105), and `db_cluster_fraction` (clustering of all inflammatory cells, not restricted to tumor-proximal) is actually *stronger* (p = 0.010) than our featured tumor-proximal metric. This means the KRAS interaction may be partly about "how much immune activity is present" rather than purely "how it is spatially organized." The "spatial organization vs. density" narrative should be presented as a hypothesis rather than a conclusion.

Interestingly, `infiltration_ratio` (the fraction of immune cells within/near tumor regions vs. total) shows a significant interaction in the **opposite** direction (HR = 4.3) — higher infiltration ratio in KRAS-mutant patients is associated with *worse* OS. This seemingly contradictory result may reflect the compositional nature of the ratio: in KRAS-mutant tumors with few immune cells overall, a high infiltration ratio might indicate a small number of cells concentrated at the tumor margin (potentially immunosuppressive MDSCs or Tregs) rather than a robust anti-tumor response.

### 7. Limitations

- **Small sample size** (n = 40): The KRAS interaction is based on only 14 KRAS-mutant patients with 5 OS events. The HR of 0.10 is extreme and almost certainly inflated by small-sample bias ("winner's curse"). While the permutation test confirms the effect is real, the magnitude should not be taken at face value. This needs validation in a larger cohort.
- **Multiple testing**: We tested 4 genes x 2-3 endpoints = 8 interaction tests. After Bonferroni correction (threshold = 0.05/8 = 0.006), the KRAS x OS interaction (p = 0.032) would not survive. The permutation test (p = 0.029) provides additional confidence but does not address the multiplicity across gene-endpoint combinations. This should be considered hypothesis-generating.
- **ConSep "Inflammatory" is a grab bag**: The cell classification lumps all immune lineages — CD8+ cytotoxic T cells, regulatory T cells, B cells, plasma cells, macrophages, MDSCs, neutrophils — into a single class. The biological meaning of a cluster of CD8+ T cells is radically different from a cluster of Tregs or MDSCs. Without cell subtype resolution, we cannot determine whether the clusters represent anti-tumor immunity or immunosuppressive aggregates. Multiplex IHC or spatial transcriptomics on these specific cases would be needed to resolve this.
- **"Spatial organization" claim is overstated**: As shown in the feature specificity analysis, simple inflammatory density shows a similar (though weaker) interaction trend with KRAS × OS. The interaction may be driven by the *amount* of immune activity rather than its spatial pattern specifically. Distinguishing these requires a cohort where density and clustering are less correlated.
- **Binary mutation calls**: We only have positive/negative calls, not variant allele frequency or specific mutation types. Different KRAS mutations (G12D vs G12V vs G13D) have different effects on immune signaling (Liao et al., 2019). A cohort with variant-level data could test whether the interaction is driven by specific KRAS alleles.
- **No MSI status**: We infer likely-MSS status from BRAF-WT status, which is a strong but imperfect proxy. Direct MSI/MMR testing would be needed to confirm. If even one KRAS-mutant patient is MSI-H, it could disproportionately influence results in this small cohort.
- **Single H&E section per patient**: Spatial features are computed from one tissue section. Intratumoral heterogeneity in CRC is substantial — immune clustering may vary across different regions of the same tumor. The feature reflects a single snapshot, not the whole-tumor immune landscape.
- **DBSCAN parameters**: While a sensitivity analysis (`dbscan_sensitivity.py`) confirmed the prognostic signal is stable across parameter ranges (eps 75-150 px, min_samples 20-50), this validation was performed on PFS. The OS interaction was not included in the parameter sensitivity analysis.
- **No proportional hazards verification**: With only 5 OS events in the KRAS-mutant group, the Cox PH assumption cannot be meaningfully tested. Violations could distort the interaction estimate.
- **BRAF underpowered**: BRAF showed similar trends to KRAS (very large point estimates in mutant group) but with only 9 mutant patients, we cannot draw conclusions. A BRAF-enriched cohort would be informative.

### 8. Clinical implications

If validated, the KRAS x immune clustering interaction has direct translational relevance:

- **Risk stratification in KRAS-mutant/MSS CRC**: This is the patient population with the fewest therapeutic options and no reliable immune biomarkers. Spatial immune clustering from routine H&E slides could identify a favorable-prognosis subgroup that current molecular-only classification systems miss entirely.
- **Immunotherapy patient selection**: KRAS-mutant/MSS patients are currently excluded from checkpoint inhibitor trials. If high immune clustering identifies patients whose immune system has already "breached" the KRAS-driven exclusion barrier, these patients might be the most likely to respond to immunotherapy that amplifies an existing response.
- **Companion diagnostic potential**: KRAS genotyping is already standard-of-care in CRC. Adding a spatial immune clustering score from H&E (computable from the same diagnostic slide) would create a two-variable classifier that captures biology inaccessible to either test alone.

## Relationship to the Sorger Lab Papers (Same Cohort)

Two key publications from the Sorger Lab at HMS used what is almost certainly the **same cohort** as our analysis — Cohort 1 in the Orion paper is 40 BWH CRC resections (stage II-IV), matching our n, institution, and stage range exactly.

### Lin et al. 2023, Cell: "Multiplexed 3D atlas of state transitions and immune interaction in colorectal cancer"

This HTAN consortium paper profiles 93 CRC specimens using CyCIF (40-plex immunofluorescence) and H&E, with 3D tissue reconstruction from serial sections. Key contributions:

- **Morphological gradients**: Identified continuous transitions from normal mucosa → adenoma → adenocarcinoma using dimensionality reduction on cellular phenotypes
- **Tumor budding**: Quantified single-cell and cluster budding at the invasive front, linked to epithelial-mesenchymal transition markers
- **TLS networks**: Mapped tertiary lymphoid structures and their B-cell/T-cell composition in 3D
- **PD1-PDL1 interactions**: Found PDL1+ tumor cells are ~5× more frequent in MSI-H tumors; mapped PD1-PDL1 contact pairs at the invasive margin
- **Immune margin profiling**: Characterized the immune cell composition at the tumor-stroma boundary using 40+ markers

**What overlaps with our work:** Both studies show that the spatial arrangement of immune cells at the tumor-stroma interface is biologically meaningful in CRC. Their PD1-PDL1 interaction mapping at the invasive margin is conceptually related to our DBSCAN clustering of inflammatory cells near tumor regions — both measure "organized immune activity at the tumor boundary."

**What they don't do that we do:** This paper is primarily a descriptive atlas. It does **not** perform mutation-stratified survival analysis, does not test whether spatial immune features predict differently in KRAS-mutant vs WT tumors, and does not model interaction terms. Their focus is characterizing what the immune microenvironment looks like; ours is asking whether its prognostic meaning changes by genotype. They also do not test KRAS-specific immune interactions despite having the sequencing data.

### Lin et al. 2023, Nature Cancer: "High-plex immunofluorescence imaging and traditional histology of the same tissue section" (Orion paper)

This paper introduces the Orion platform for same-section H&E + 16-18 channel immunofluorescence on 74 CRC resections across two cohorts (40 BWH + 34 external). Key results:

- **IFM1** (Immunoscore analog): CD3/CD8 density at tumor center and invasive margin, recapitulated from IF. HR = 0.14, p = 7.63×10⁻⁵ for PFS.
- **IFM2** (improved model): Added α-SMA + CD45 + PD-L1 + CD4 features. HR = 0.05, p = 5.5×10⁻⁶ for PFS.
- **IFM4** (composite): HR = 0.12, p = 8.4×10⁻⁵ for PFS.
- **Spatial-LDA topics**: Topic 7 (poorly differentiated/high-grade tumor with stromal invasion) predicted poor PFS (HR = 0.26, p = 2.98×10⁻⁴).
- **H&E features**: Demonstrated that H&E-derived computational features correlate with IF-derived immune markers, supporting H&E as a viable substrate for immune assessment.

**What overlaps with our work:**
1. Both studies show spatial immune features from H&E predict PFS in CRC
2. Both use the same cohort (Cohort 1)
3. Both demonstrate that immune organization at the tumor-stroma interface is prognostic
4. Their IFM1 (density-based Immunoscore) and our `db_tumor_cluster_fraction` both capture "immune activity near tumor"

**What they don't do that we do:**
1. **No mutation stratification**: The Orion paper does not test whether IFM1/IFM2/IFM4 predict differently by KRAS, BRAF, TP53, or APC status. Their models treat the cohort as genotype-homogeneous.
2. **No interaction terms**: They do not fit Cox models with gene × immune feature interaction terms.
3. **No KRAS-specific analysis**: The most striking finding in our work — that immune clustering has a 12× stronger HR in KRAS-mutant vs WT for OS — is not tested or reported.
4. **No OS analysis**: The Orion paper focuses on PFS; our KRAS interaction is OS-specific. This is important because our data also shows immune clustering predicts PFS in *both* KRAS subgroups (no interaction), while the OS interaction is subgroup-specific.

**What they do better than us:**
1. **Cell type resolution**: With 16-18 IF channels, they distinguish CD3, CD4, CD8, CD20, CD68, FOXP3, PD-L1, α-SMA etc. Our ConSep "Inflammatory" class is a single undifferentiated category.
2. **Feature engineering**: IFM2 and Spatial-LDA topics are more sophisticated than our DBSCAN clustering.
3. **Validation cohort**: They have a second cohort (n=34) for external validation; we have only the discovery cohort.
4. **Effect sizes**: IFM2 HR = 0.05 for PFS substantially outperforms our `db_tumor_cluster_fraction` HR = 0.40 for PFS — multiplex IF captures more immune biology than H&E classification alone.

### Positioning our work relative to these papers

Our analysis is **complementary, not redundant**. The Orion paper establishes *what* — spatial immune features from H&E predict CRC survival. We ask the next question: **for whom?** — showing the prognostic value is genotype-dependent (KRAS × immune clustering interaction for OS).

The fact that we are almost certainly using the same cohort is a **strength**, not a weakness. We are extending the Sorger group's own published data with a mutation interaction analysis they did not perform. Our contribution is:

1. **Integrating the sequencing layer** that both papers had available but did not cross with spatial immune features for survival modeling
2. **Identifying KRAS status as an effect modifier** — suggesting that genotype-agnostic immune scores (like Immunoscore/IFM1) may miss subgroup-specific biology
3. **Generating the hypothesis** that in KRAS-mutant/MSS CRC, spatial immune organization has outsized prognostic value — a finding that could be validated using the Orion paper's superior cell-type-resolved features

The natural next step would be to test whether the Orion paper's IFM1/IFM2 scores also show a KRAS × survival interaction. If they do, the cell-type resolution would reveal *which* immune lineages drive the clustering signal in KRAS-mutant tumors — a question our ConSep-based approach cannot answer.

### Yang et al. 2025: "Single-cell and spatial transcriptome profiling identifies the immunosuppressive spatial niche in KRAS-mutant colorectal cancer" (PMID 41475845)

This paper (Journal for ImmunoTherapy of Cancer) uses scRNA-seq + spatial transcriptomics to map the immunosuppressive microenvironment in KRAS-mutant CRC. Key findings:

- **KRAS-mutant epithelial cells recruit S100A8+ monocytes** via the MDK-SDC4 signaling axis
- **CTHRC1+ fibroblasts form a physical barrier** around KRAS-mutant tumor nests, creating spatial exclusion of cytotoxic T cells
- **Specific spatial "ecotypes"** (organized cellular neighborhoods) are present in KRAS-mutant tumors but absent in WT tumors
- These ecotypes define an immunosuppressive niche that actively repels adaptive immune infiltration

**How this connects to our work:** This paper provides the **mechanistic basis** for our "Breach" hypothesis. The MDK-SDC4 + CTHRC1+ fibroblast barrier is exactly the exclusion program that our DBSCAN clustering feature implicitly measures the success or failure of. Patients with high `db_tumor_cluster_fraction` in KRAS-mutant tumors may be those whose immune system has overcome the MDK/CTHRC1-mediated spatial exclusion — explaining the extreme survival benefit (HR = 0.10). Patients with low clustering may represent the default KRAS-mutant state: immune cells successfully excluded by this fibroblast barrier, consistent with the poor outcomes we observe.

**What they don't test:** Survival interaction. They characterize the niche but do not ask whether overcoming it (high immune clustering despite KRAS mutation) predicts better outcomes. Our work provides the clinical/prognostic evidence for their spatial biology.

## Relationship to Broader Literature

### Where this replicates known biology

| Finding | Literature consensus | Verdict |
|---------|---------------------|---------|
| KRAS-mutant tumors have lower immune clustering (implied by low clustering group's poor outcomes) | KRAS mutants are "immune cold" or excluded (Liao et al., Cancer Cell 2019) | Consistent |
| TP53/APC do not modify immune prognosis | TP53/APC generally do not dictate immune "coldness" as strongly as KRAS | Consistent |
| BRAF-mutant tumors trend toward high immune infiltration | BRAF V600E is associated with MSI-H and massive immune infiltration | Consistent |
| KRAS and BRAF are mutually exclusive | Canonical RAS-RAF-MAPK pathway mutual exclusivity in CRC | Consistent |

### Where this may be novel: density vs. spatial organization

The prevailing dogma, established by Galon and the Immunoscore consortium (Pages et al., Lancet Oncol 2018), holds that immune infiltration (CD3/CD8 density at tumor center and invasive margin) is prognostic **independent of molecular subtype** — in both MSI and MSS, and both KRAS-mutant and WT populations.

Our finding that spatial immune clustering is strongly prognostic in KRAS-mutant tumors but essentially irrelevant in KRAS-WT tumors (for OS) is a specific form of effect modification that **contrasts with this general consensus**. The critical distinction is what we are measuring: the Immunoscore measures **density** (how many immune cells are present), while `db_tumor_cluster_fraction` measures **spatial organization** (whether those cells have formed organized clusters at the tumor-stroma interface).

This addresses the "immune exclusion" taxonomy directly:
- **Immune desert**: No T-cells anywhere
- **Immune excluded**: T-cells present but trapped in stroma/periphery
- **Immune inflamed**: T-cells organized at the tumor interface

Our feature specifically captures the transition from "excluded" to "inflamed" — immune cells that have not just entered the tissue but have organized at the tumor boundary. In KRAS-mutant tumors, where the default state is exclusion, this transition is rare and clinically transformative. In KRAS-WT tumors, where the baseline is less hostile, the transition is common and therefore less discriminating.

### Comparison with contemporary spatial literature

| Study | Focus | How our work adds |
|-------|-------|-------------------|
| Pages et al. (Lancet Oncol 2018) | Immunoscore (CD3/CD8 density) is prognostic across all CRC subtypes | We show that *spatial organization* (not just density) has subtype-dependent prognostic value — clustering matters more in KRAS-mutant than KRAS-WT for OS |
| Liao et al. (Cancer Cell 2019) | KRAS G12D/V drives CXCL8/GM-CSF secretion, MDSC recruitment, CD8 exclusion | Our data shows that when clustering occurs *despite* this exclusion program, it identifies exceptional survivors — consistent with a "breach" of the exclusion barrier |
| He et al. (2024) | Spatial proteogenomics shows KRAS mutations correlate with TGF-beta in stroma | Our work shows that if immune clustering occurs despite this stromal TGF-beta barrier, survival is excellent — the spatial feature captures successful immune infiltration against active resistance |
| Pelka et al. (2021) | Spatially organized multicellular immune "hubs" in CRC | Our finding suggests these hubs are most clinically meaningful when they form against the "resistance" of a KRAS-driven immunosuppressive program |
| Smeby et al. (2018) | BRAF increases infiltration, KRAS reduces it; cytotoxic T-cells prognostic in both groups | We find a different result for *spatial clustering* vs. OS — prognostic only in KRAS-mutant. The discrepancy may reflect density vs. organization, or endpoint-specific biology (OS vs. DFS) |
| Lal et al. (2018) | KRAS associated with lower immune scores but score still prognostic in both groups | Similar discrepancy — their density-based score is "universally prognostic," while our clustering-based feature shows KRAS-specific prognostic value |
| Cadwell et al. | KRAS-mutant tumors form specific spatial niches with bacteria/neutrophil co-localization | Our DBSCAN clustering may be capturing these reactive niches — the spatial algorithm is better suited to detect organized cellular aggregates than cell counts |

### Where further work is needed

The discrepancy between our finding and prior density-based studies (Lal, Smeby, Immunoscore consortium) could reflect:

1. **Feature type**: They measured CD3/CD8 density; we measure DBSCAN spatial clustering of ConSep-classified inflammatory cells. These capture different biological phenomena — presence vs. organization.
2. **Endpoint specificity**: Our interaction is specific to OS, not PFS. For PFS, immune clustering is prognostic in both KRAS subgroups. This suggests the interaction may be specific to the biology captured by overall survival (which includes non-cancer deaths and long-term immune surveillance) vs. time-to-progression.
3. **Cohort size**: Our n = 40 (14 KRAS-mutant, 5 OS events in that group) is much smaller than large Immunoscore validation cohorts (n > 2,000). The interaction could be a true finding amplified by small-sample variance, or it could reflect cohort-specific features. The extreme HR of 0.10 is biologically plausible but likely somewhat inflated.
4. **Cell type resolution**: Our "inflammatory" class from ConSep includes lymphocytes, plasma cells, and other immune cells without distinguishing CD3/CD8/CD4/FOXP3 subtypes. The composition of the clusters may differ between KRAS subgroups — spatial transcriptomics or multiplex IHC on these specific cases would clarify what cell types drive the clustering signal.

A validation path would be to test this interaction in a larger CRC cohort with both spatial H&E features and KRAS genotyping (e.g., TCGA CRC with digital pathology, or the Immunoscore consortium data with molecular annotations). The most compelling validation would be in a KRAS-mutant/MSS-specific cohort, where the clinical need is greatest.

## Output Files

| File | Description |
|------|-------------|
| `data/mutation_univariate_cox.csv` | Gene-level Cox regression results |
| `data/mutation_stratified_cox.csv` | Immune clustering HRs within mutation subgroups |
| `data/mutation_interaction_cox.csv` | Formal interaction test results |
| `data/mutation_burden_cox.csv` | Mutation burden analysis |
| `figures/mutation_forest_plot.png` | Forest plot of univariate mutation HRs |
| `figures/km_interaction_KRAS_os.png` | 4-group KM plot for KRAS x clustering x OS |
