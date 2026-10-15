#!/usr/bin/env python3
"""
Spatial feature survival analysis for HTAN CRC cohort.

Tests continuous spatial features (inflammatory density, tumor-stroma ratio, NLR)
in Cox regression against OS and PFS endpoints. Features are standardized (z-scored)
so hazard ratios represent the change per standard deviation increase.

Usage:
    python scripts/survival_analysis.py
"""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*delta method.*")

sys.path.insert(0, str(Path(__file__).parent))
from load_cell_densities import get_project_root, load_cell_densities

FEATURES = [
    # Clinical
    "age_at_diagnosis",
    "sex_male",
    "stage_binary",
    "tumor_grade_high",
    # Inflammatory (ConSep spatial)
    "inflammatory_pct",
    "infiltration_ratio",
    "tumor_associated_pct",
    "density_overall",
    "tumor_stroma_ratio",
    # CoNIC-derived
    "nlr",
    "eosinophil_density",
    # Truly spatial
    "nn_dist_median",
    "db_cluster_fraction",
    "db_tumor_cluster_fraction",
    "cluster_fraction",
    "tumor_cluster_fraction",
    "tumor_immune_morans_i",
]

FEATURE_LABELS = {
    "age_at_diagnosis": "Age at diagnosis",
    "sex_male": "Sex (male)",
    "stage_binary": "Stage III-IV",
    "tumor_grade_high": "Tumor grade (high)",
    "inflammatory_pct": "Inflammatory %",
    "infiltration_ratio": "Infiltration ratio",
    "tumor_associated_pct": "Tumor-associated %",
    "density_overall": "Inflammatory density",
    "tumor_stroma_ratio": "Tumor-stroma ratio",
    "nlr": "Neutrophil-lymphocyte ratio",
    "eosinophil_density": "Eosinophil density",
    "nn_dist_median": "Immune-tumor distance (median)",
    "db_cluster_fraction": "DBSCAN cluster frac (all)",
    "db_tumor_cluster_fraction": "DBSCAN cluster frac (tumor)",
    "cluster_fraction": "HDBSCAN cluster frac (all)",
    "tumor_cluster_fraction": "HDBSCAN cluster frac (tumor)",
    "tumor_immune_morans_i": "Tumor immune Moran's I",
}

ENDPOINTS = {
    "os": {"time": "os_time", "event": "os_event", "label": "Overall Survival"},
    "css": {
        "time": "css_time",
        "event": "css_event",
        "label": "Cancer-Specific Survival",
    },
    "pfs": {
        "time": "pfs_time",
        "event": "pfs_event",
        "label": "Progression-Free Survival",
    },
}


def build_survival_dataset() -> pd.DataFrame:
    """Build merged analysis-ready DataFrame with all features and survival endpoints."""
    project_root = get_project_root()

    # 1. Load inflammatory metrics (per-sample) → aggregate to participant level
    print("  Loading inflammatory metrics...")
    metrics = pd.read_csv(project_root / "data" / "immune_phenotype_metrics_v2.csv")
    feature_cols = [
        "inflammatory_pct",
        "infiltration_ratio",
        "tumor_associated_pct",
        "density_overall",
    ]
    metrics_agg = metrics.groupby("participant_id")[feature_cols].mean().reset_index()

    # 2. Load phenotype classifications (take first per participant)
    print("  Loading phenotype classifications...")
    phenotypes = pd.read_csv(
        project_root / "data" / "immune_phenotype_classification_v2.csv"
    )
    phenotypes_agg = (
        phenotypes.groupby("participant_id")["phenotype"].first().reset_index()
    )

    # 3. CoNIC cell densities → NLR and eosinophil density
    print("  Loading CoNIC cell densities...")
    conic = load_cell_densities("conic")
    conic_counts = (
        conic.groupby(["participant_id", "cell_class"])["count"].sum().reset_index()
    )
    conic_wide = conic_counts.pivot_table(
        index="participant_id", columns="cell_class", values="count", aggfunc="sum"
    ).fillna(0)

    nlr_series = conic_wide["Neutrophil"] / conic_wide["Lymphocyte"].replace(
        0, np.nan
    )
    nlr_df = nlr_series.reset_index()
    nlr_df.columns = ["participant_id", "nlr"]

    # Eosinophil density: average density across samples per participant
    eos = conic[conic["cell_class"] == "Eosinophil"]
    eos_density = eos.groupby("participant_id")["density"].mean().reset_index()
    eos_density.columns = ["participant_id", "eosinophil_density"]

    # 4. ConSep cell densities → tumor-stroma ratio
    print("  Loading ConSep cell densities...")
    consep = load_cell_densities("consep")
    consep_counts = (
        consep.groupby(["participant_id", "cell_class"])["count"].sum().reset_index()
    )
    consep_wide = consep_counts.pivot_table(
        index="participant_id", columns="cell_class", values="count", aggfunc="sum"
    ).fillna(0)

    tsr_series = consep_wide["Malignant epithelial"] / consep_wide["Stroma"].replace(
        0, np.nan
    )
    tsr_df = tsr_series.reset_index()
    tsr_df.columns = ["participant_id", "tumor_stroma_ratio"]

    # 5. Spatial features (NN distances + Moran's I)
    print("  Loading spatial features...")
    spatial_path = project_root / "data" / "spatial_features.csv"
    if spatial_path.exists():
        spatial_features = pd.read_csv(spatial_path)
    else:
        print("    WARNING: spatial_features.csv not found — run compute_spatial_features.py first")
        spatial_features = None

    # 6. Clinical metadata
    print("  Loading clinical metadata...")
    clinical = pd.read_csv(project_root / "metadata" / "participant.tsv", sep="\t")
    clinical = clinical.rename(columns={"HTAN Participant ID": "participant_id"})
    survival_cols = [
        "participant_id",
        "Vital Status",
        "Days to Last Followup",
        "Days to Last Known Disease Status",
        "Progression or Recurrence",
        "AJCC Pathologic Stage",
        "Cause of Death",
        "Age at Diagnosis (years)",
        "Gender",
        "Tumor Grade",
    ]
    clinical = clinical[survival_cols].copy()

    # 7. Merge everything
    print("  Merging datasets...")
    df = metrics_agg.merge(phenotypes_agg, on="participant_id", how="outer")
    df = df.merge(nlr_df, on="participant_id", how="left")
    df = df.merge(eos_density, on="participant_id", how="left")
    df = df.merge(tsr_df, on="participant_id", how="left")
    if spatial_features is not None:
        df = df.merge(spatial_features, on="participant_id", how="left")
    df = df.merge(clinical, on="participant_id", how="left")

    # 8. Construct survival endpoints
    print("  Constructing survival endpoints...")

    # OS (all-cause): time = Days to Last Followup, event = Dead
    df["os_time"] = pd.to_numeric(df["Days to Last Followup"], errors="coerce")
    df["os_event"] = (df["Vital Status"] == "Dead").astype(int)
    # Exclude unknown vital status from OS
    unknown_mask = df["Vital Status"] == "unknown"
    df.loc[unknown_mask, ["os_time", "os_event"]] = np.nan

    # CSS (cancer-specific): same time, but only cancer deaths count as events
    # Non-cancer and unknown-cause deaths are censored at their followup time
    df["css_time"] = df["os_time"].copy()
    df["css_event"] = (
        (df["Vital Status"] == "Dead")
        & (df["Cause of Death"] == "Cancer Related")
    ).astype(int)
    df.loc[unknown_mask, ["css_time", "css_event"]] = np.nan

    # PFS: time = Days to Last Known Disease Status, event = progression
    df["pfs_time"] = pd.to_numeric(
        df["Days to Last Known Disease Status"], errors="coerce"
    )
    df["pfs_event"] = (
        df["Progression or Recurrence"]
        .fillna("")
        .str.startswith("Yes")
        .astype(int)
    )
    # Exclude Not Reported from PFS
    not_reported_mask = df["Progression or Recurrence"] == "Not Reported"
    df.loc[not_reported_mask, ["pfs_time", "pfs_event"]] = np.nan

    # Stage simplification: I-II vs III-IV
    def simplify_stage(stage):
        if pd.isna(stage):
            return np.nan
        if "IV" in stage or "III" in stage:
            return "III-IV"
        return "I-II"

    df["stage_group"] = df["AJCC Pathologic Stage"].apply(simplify_stage)
    df["stage_binary"] = df["stage_group"].map({"I-II": 0, "III-IV": 1}).astype(float)

    # Clinical features
    df["age_at_diagnosis"] = pd.to_numeric(df["Age at Diagnosis (years)"], errors="coerce")
    df["sex_male"] = (df["Gender"] == "male").astype(float)
    df["tumor_grade_high"] = df["Tumor Grade"].map(
        {"High Grade": 1.0, "Low Grade": 0.0}
    )  # "Not Reported" → NaN

    # Save
    output_path = project_root / "data" / "survival_analysis_dataset.csv"
    df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path} ({len(df)} rows)")

    return df


def run_cox_univariate(df: pd.DataFrame) -> pd.DataFrame:
    """Run univariate Cox regression for each feature x endpoint.

    Features are z-scored so HRs represent the effect per SD increase.
    """
    results = []

    for feature in FEATURES:
        for ep_name, ep in ENDPOINTS.items():
            cols = [feature, ep["time"], ep["event"]]
            subset = df[cols].dropna().copy()

            if len(subset) < 10:
                print(f"    Skipping {feature} x {ep_name}: n={len(subset)}")
                continue

            # Standardize for interpretable per-SD HRs
            mean_val = subset[feature].mean()
            std_val = subset[feature].std()
            if std_val == 0:
                continue
            subset[feature] = (subset[feature] - mean_val) / std_val

            try:
                cph = CoxPHFitter()
                cph.fit(subset, duration_col=ep["time"], event_col=ep["event"])
                summary = cph.summary

                results.append(
                    {
                        "feature": feature,
                        "feature_label": FEATURE_LABELS[feature],
                        "endpoint": ep_name,
                        "endpoint_label": ep["label"],
                        "n": len(subset),
                        "n_events": int(subset[ep["event"]].sum()),
                        "hr": summary.loc[feature, "exp(coef)"],
                        "ci_lower": summary.loc[feature, "exp(coef) lower 95%"],
                        "ci_upper": summary.loc[feature, "exp(coef) upper 95%"],
                        "p_value": summary.loc[feature, "p"],
                        "concordance": cph.concordance_index_,
                    }
                )
            except Exception as e:
                print(f"    Error fitting {feature} x {ep_name}: {e}")

    return pd.DataFrame(results)


def plot_kaplan_meier(df: pd.DataFrame, output_dir: Path):
    """Plot KM curves for each feature (median split) and by phenotype."""
    for feature in FEATURES:
        for ep_name, ep in ENDPOINTS.items():
            cols = [feature, ep["time"], ep["event"]]
            subset = df[cols].dropna()

            if len(subset) < 10:
                continue

            median_val = subset[feature].median()
            high = subset[subset[feature] >= median_val]
            low = subset[subset[feature] < median_val]

            if len(high) < 2 or len(low) < 2:
                continue

            fig, ax = plt.subplots(figsize=(8, 6))

            kmf_high = KaplanMeierFitter()
            kmf_low = KaplanMeierFitter()

            kmf_high.fit(
                high[ep["time"]], high[ep["event"]], label=f"High (n={len(high)})"
            )
            kmf_low.fit(
                low[ep["time"]], low[ep["event"]], label=f"Low (n={len(low)})"
            )

            kmf_high.plot_survival_function(ax=ax, ci_show=True, color="#b2182b")
            kmf_low.plot_survival_function(ax=ax, ci_show=True, color="#2166ac")

            lr = logrank_test(
                high[ep["time"]],
                low[ep["time"]],
                high[ep["event"]],
                low[ep["event"]],
            )

            ax.set_title(
                f"{FEATURE_LABELS[feature]} — {ep['label']}\n"
                f"Median split = {median_val:.4g}, log-rank p = {lr.p_value:.3f}",
                fontsize=12,
            )
            ax.set_xlabel("Days", fontsize=11)
            ax.set_ylabel("Survival probability", fontsize=11)
            ax.legend(loc="lower left", fontsize=10)

            plt.tight_layout()
            plt.savefig(
                output_dir / f"km_{feature}_{ep_name}.png",
                dpi=150,
                bbox_inches="tight",
            )
            plt.close()

    # KM by immune phenotype category
    phenotype_colors = {
        "Immune-desert": "#2166ac",
        "Altered-excluded": "#92c5de",
        "Altered-immunosuppressed": "#f4a582",
        "Hot": "#b2182b",
    }

    for ep_name, ep in ENDPOINTS.items():
        cols = ["phenotype", ep["time"], ep["event"]]
        subset = df[cols].dropna()

        if len(subset) < 10:
            continue

        fig, ax = plt.subplots(figsize=(8, 6))

        for phenotype, color in phenotype_colors.items():
            group = subset[subset["phenotype"] == phenotype]
            if len(group) < 2:
                continue
            kmf = KaplanMeierFitter()
            kmf.fit(
                group[ep["time"]],
                group[ep["event"]],
                label=f"{phenotype} (n={len(group)})",
            )
            kmf.plot_survival_function(ax=ax, ci_show=False, color=color)

        ax.set_title(f"Immune Phenotype — {ep['label']}", fontsize=12)
        ax.set_xlabel("Days", fontsize=11)
        ax.set_ylabel("Survival probability", fontsize=11)
        ax.legend(loc="lower left", fontsize=10)

        plt.tight_layout()
        plt.savefig(
            output_dir / f"km_phenotype_{ep_name}.png",
            dpi=150,
            bbox_inches="tight",
        )
        plt.close()


def run_cox_multivariate(
    df: pd.DataFrame, univariate_results: pd.DataFrame
) -> pd.DataFrame:
    """Run stage-adjusted Cox for features with univariate p < 0.2."""
    # Skip clinical variables from stage adjustment (can't adjust stage for stage,
    # and other clinical vars are confounders not spatial features)
    clinical_features = {"age_at_diagnosis", "sex_male", "stage_binary", "tumor_grade_high"}
    sig = univariate_results[
        (univariate_results["p_value"] < 0.2)
        & (~univariate_results["feature"].isin(clinical_features))
    ][["feature", "endpoint"]].drop_duplicates()

    if len(sig) == 0:
        print("  No features with p < 0.2 for multivariate analysis")
        return pd.DataFrame()

    results = []

    for _, row in sig.iterrows():
        feature = row["feature"]
        ep_name = row["endpoint"]
        ep = ENDPOINTS[ep_name]

        cols = [feature, "stage_binary", ep["time"], ep["event"]]
        subset = df[cols].dropna().copy()

        if len(subset) < 10:
            continue

        # Standardize feature (same as univariate)
        mean_val = subset[feature].mean()
        std_val = subset[feature].std()
        if std_val == 0:
            continue
        subset[feature] = (subset[feature] - mean_val) / std_val

        try:
            cph = CoxPHFitter()
            cph.fit(subset, duration_col=ep["time"], event_col=ep["event"])
            summary = cph.summary

            results.append(
                {
                    "feature": feature,
                    "feature_label": FEATURE_LABELS[feature],
                    "endpoint": ep_name,
                    "endpoint_label": ep["label"],
                    "n": len(subset),
                    "adjusted_hr": summary.loc[feature, "exp(coef)"],
                    "adjusted_ci_lower": summary.loc[
                        feature, "exp(coef) lower 95%"
                    ],
                    "adjusted_ci_upper": summary.loc[
                        feature, "exp(coef) upper 95%"
                    ],
                    "adjusted_p": summary.loc[feature, "p"],
                    "stage_hr": summary.loc["stage_binary", "exp(coef)"],
                    "stage_p": summary.loc["stage_binary", "p"],
                }
            )
        except Exception as e:
            print(f"    Error fitting {feature} + stage x {ep_name}: {e}")

    return pd.DataFrame(results)


def plot_forest(univariate_results: pd.DataFrame, output_path: Path):
    """Forest plot of all univariate HRs (per SD) with 95% CIs."""
    if len(univariate_results) == 0:
        return

    n_endpoints = len(ENDPOINTS)
    fig, axes = plt.subplots(1, n_endpoints, figsize=(7 * n_endpoints, 6), sharey=True)

    for idx, (ep_name, ep) in enumerate(ENDPOINTS.items()):
        ax = axes[idx]
        ep_data = univariate_results[
            univariate_results["endpoint"] == ep_name
        ].copy()
        ep_data = ep_data.sort_values("feature").reset_index(drop=True)

        if len(ep_data) == 0:
            ax.set_visible(False)
            continue

        y_pos = np.arange(len(ep_data))

        for i, (_, row) in enumerate(ep_data.iterrows()):
            color = "#b2182b" if row["p_value"] < 0.05 else "black"
            marker = "D" if row["p_value"] < 0.05 else "o"
            ax.errorbar(
                row["hr"],
                i,
                xerr=[[row["hr"] - row["ci_lower"]], [row["ci_upper"] - row["hr"]]],
                fmt=marker,
                color=color,
                capsize=4,
                markersize=7,
                linewidth=1.5,
            )

        ax.axvline(x=1, color="gray", linestyle="--", alpha=0.5)
        ax.set_xscale("log")

        labels = [
            f"{row['feature_label']} (p={row['p_value']:.3f})"
            for _, row in ep_data.iterrows()
        ]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Hazard Ratio per SD (95% CI)", fontsize=11)
        ax.set_title(ep["label"], fontsize=12)

    plt.suptitle("Univariate Cox Regression — Hazard Ratios", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    project_root = get_project_root()
    figures_dir = project_root / "figures"
    data_dir = project_root / "data"
    figures_dir.mkdir(exist_ok=True)
    data_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Spatial Feature Survival Analysis")
    print("=" * 70)

    # 1. Build dataset
    print("\n1. Building survival dataset...")
    df = build_survival_dataset()

    print(f"\n   Dataset: {len(df)} participants, {len(df.columns)} columns")
    for ep_name, ep in ENDPOINTS.items():
        n_events = df[ep["event"]].sum()
        n_eval = df[ep["time"]].notna().sum()
        print(f"   {ep['label']:30s}: {n_events:.0f} events / {n_eval} evaluable")

    print("\n   Feature ranges:")
    for f in FEATURES:
        if f in df.columns:
            vals = df[f].dropna()
            print(
                f"     {FEATURE_LABELS[f]:30s}: "
                f"{vals.min():.4f} — {vals.max():.4f} "
                f"(median {vals.median():.4f})"
            )

    # 2. Univariate Cox regression
    print("\n2. Univariate Cox regression...")
    cox_results = run_cox_univariate(df)
    cox_path = data_dir / "survival_cox_results.csv"
    cox_results.to_csv(cox_path, index=False)
    print(f"   Saved: {cox_path}")

    print("\n   Results (HR per SD increase):")
    for _, row in cox_results.iterrows():
        sig = "*" if row["p_value"] < 0.05 else " "
        print(
            f"   {sig} {row['feature_label']:30s} {row['endpoint']:3s}  "
            f"HR={row['hr']:.3f} ({row['ci_lower']:.3f}—{row['ci_upper']:.3f})  "
            f"p={row['p_value']:.3f}  C={row['concordance']:.3f}"
        )

    # 3. Kaplan-Meier curves
    print("\n3. Plotting Kaplan-Meier curves...")
    plot_kaplan_meier(df, figures_dir)
    print("   Saved KM plots to figures/")

    # 4. Stage-adjusted multivariate
    print("\n4. Stage-adjusted multivariate Cox regression...")
    mv_results = run_cox_multivariate(df, cox_results)
    if len(mv_results) > 0:
        mv_path = data_dir / "survival_cox_multivariate_results.csv"
        mv_results.to_csv(mv_path, index=False)
        print(f"   Saved: {mv_path}")

        print("\n   Stage-adjusted results (features with univariate p < 0.2):")
        for _, row in mv_results.iterrows():
            sig = "*" if row["adjusted_p"] < 0.05 else " "
            print(
                f"   {sig} {row['feature_label']:30s} {row['endpoint']:3s}  "
                f"adj.HR={row['adjusted_hr']:.3f} "
                f"({row['adjusted_ci_lower']:.3f}—{row['adjusted_ci_upper']:.3f})  "
                f"p={row['adjusted_p']:.3f}  "
                f"| stage HR={row['stage_hr']:.3f} p={row['stage_p']:.3f}"
            )

    # 5. Forest plot
    print("\n5. Plotting forest plot...")
    forest_path = figures_dir / "survival_forest_plot.png"
    plot_forest(cox_results, forest_path)
    print(f"   Saved: {forest_path}")

    print("\n" + "=" * 70)
    print("Done!")


if __name__ == "__main__":
    main()
