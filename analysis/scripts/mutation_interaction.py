#!/usr/bin/env python3
"""
Mutation × immune clustering interaction analysis for HTAN CRC cohort.

Tests whether mutation status in key CRC genes modulates the prognostic value
of tumor-proximal immune clustering (db_tumor_cluster_fraction).

Usage:
    python scripts/mutation_interaction.py
"""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*delta method.*")

sys.path.insert(0, str(Path(__file__).parent))
from load_cell_densities import get_project_root

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

FOCUS_GENES = ["KRAS", "BRAF", "TP53", "APC"]
IMMUNE_FEATURE = "db_tumor_cluster_fraction"


def load_mutation_data(project_root: Path) -> pd.DataFrame:
    """Load targeted sequencing data and pivot to wide binary format."""
    seq = pd.read_csv(project_root / "metadata" / "bquxjob_6358239c_19c48c9d3ac.csv")

    # Fix gene symbol artifacts
    seq.loc[seq["Gene_Symbol"].str.startswith("TCF7L", na=False), "Gene_Symbol"] = "TCF7L2"
    seq.loc[seq["Gene_Symbol"] == "p53", "Gene_Symbol"] = "TP53"

    # Filter to targeted sequencing results with positive/negative calls
    seq = seq[seq["Test_Result"].isin(["positive", "negative"])].copy()

    # Pivot: participant × gene, values = 1/0
    seq["mutation"] = (seq["Test_Result"] == "positive").astype(int)
    mut_wide = seq.pivot_table(
        index="HTAN_Participant_ID",
        columns="Gene_Symbol",
        values="mutation",
        aggfunc="max",
    ).reset_index()
    mut_wide = mut_wide.rename(columns={"HTAN_Participant_ID": "participant_id"})

    return mut_wide


def load_therapy_data(project_root: Path) -> pd.DataFrame:
    """Load therapy metadata and create treated binary."""
    therapy = pd.read_csv(project_root / "metadata" / "bquxjob_4136c6be_19c48c7c3e6.csv")
    therapy["treated"] = (therapy["Treatment_or_Therapy"] == "yes").astype(int)
    therapy = therapy[["HTAN_Participant_ID", "treated"]].rename(
        columns={"HTAN_Participant_ID": "participant_id"}
    )
    return therapy


def run_univariate_mutation_cox(df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Analysis 1: Cox PH for each gene (binary) × endpoint."""
    results = []

    for gene in genes:
        if gene not in df.columns:
            continue
        for ep_name, ep in ENDPOINTS.items():
            cols = [gene, ep["time"], ep["event"]]
            subset = df[cols].dropna().copy()

            if len(subset) < 10:
                continue

            # Check ≥5 events in both mutant and WT groups
            events_mut = subset.loc[subset[gene] == 1, ep["event"]].sum()
            events_wt = subset.loc[subset[gene] == 0, ep["event"]].sum()
            if events_mut < 5 or events_wt < 5:
                continue

            try:
                cph = CoxPHFitter()
                cph.fit(subset, duration_col=ep["time"], event_col=ep["event"])
                s = cph.summary

                results.append({
                    "gene": gene,
                    "endpoint": ep_name,
                    "endpoint_label": ep["label"],
                    "n": len(subset),
                    "n_mutant": int(subset[gene].sum()),
                    "n_wt": int((subset[gene] == 0).sum()),
                    "n_events": int(subset[ep["event"]].sum()),
                    "events_mutant": int(events_mut),
                    "events_wt": int(events_wt),
                    "hr": s.loc[gene, "exp(coef)"],
                    "ci_lower": s.loc[gene, "exp(coef) lower 95%"],
                    "ci_upper": s.loc[gene, "exp(coef) upper 95%"],
                    "p_value": s.loc[gene, "p"],
                    "concordance": cph.concordance_index_,
                })
            except Exception as e:
                print(f"    Error fitting {gene} x {ep_name}: {e}")

    return pd.DataFrame(results)


def run_stratified_cox(df: pd.DataFrame, focus_genes: list[str]) -> pd.DataFrame:
    """Analysis 2: Cox on z-scored immune clustering within mutation subgroups."""
    results = []

    for gene in focus_genes:
        if gene not in df.columns:
            continue
        for ep_name, ep in ENDPOINTS.items():
            cols = [IMMUNE_FEATURE, gene, ep["time"], ep["event"]]
            subset = df[cols].dropna().copy()

            for mut_val, stratum_label in [(1, "mutant"), (0, "WT")]:
                stratum = subset[subset[gene] == mut_val].copy()

                if len(stratum) < 5 or stratum[ep["event"]].sum() < 3:
                    continue

                # Z-score within stratum
                mean_val = stratum[IMMUNE_FEATURE].mean()
                std_val = stratum[IMMUNE_FEATURE].std()
                if std_val == 0:
                    continue
                stratum[IMMUNE_FEATURE] = (stratum[IMMUNE_FEATURE] - mean_val) / std_val

                try:
                    cph = CoxPHFitter()
                    cph.fit(
                        stratum[[IMMUNE_FEATURE, ep["time"], ep["event"]]],
                        duration_col=ep["time"],
                        event_col=ep["event"],
                    )
                    s = cph.summary

                    results.append({
                        "gene": gene,
                        "stratum": stratum_label,
                        "endpoint": ep_name,
                        "endpoint_label": ep["label"],
                        "n": len(stratum),
                        "n_events": int(stratum[ep["event"]].sum()),
                        "hr": s.loc[IMMUNE_FEATURE, "exp(coef)"],
                        "ci_lower": s.loc[IMMUNE_FEATURE, "exp(coef) lower 95%"],
                        "ci_upper": s.loc[IMMUNE_FEATURE, "exp(coef) upper 95%"],
                        "p_value": s.loc[IMMUNE_FEATURE, "p"],
                    })
                except Exception as e:
                    print(f"    Error fitting stratified {gene}/{stratum_label} x {ep_name}: {e}")

    return pd.DataFrame(results)


def run_interaction_test(df: pd.DataFrame, focus_genes: list[str]) -> pd.DataFrame:
    """Analysis 3: Formal interaction test (cluster × mutation)."""
    results = []

    for gene in focus_genes:
        if gene not in df.columns:
            continue
        for ep_name, ep in ENDPOINTS.items():
            cols = [IMMUNE_FEATURE, gene, ep["time"], ep["event"]]
            subset = df[cols].dropna().copy()

            if len(subset) < 15 or subset[ep["event"]].sum() < 8:
                continue

            # Z-score the immune feature
            mean_val = subset[IMMUNE_FEATURE].mean()
            std_val = subset[IMMUNE_FEATURE].std()
            if std_val == 0:
                continue
            subset["cluster_z"] = (subset[IMMUNE_FEATURE] - mean_val) / std_val

            # Create interaction term
            subset["interaction"] = subset["cluster_z"] * subset[gene]

            try:
                cph = CoxPHFitter()
                fit_cols = ["cluster_z", gene, "interaction", ep["time"], ep["event"]]
                cph.fit(subset[fit_cols], duration_col=ep["time"], event_col=ep["event"])
                s = cph.summary

                results.append({
                    "gene": gene,
                    "endpoint": ep_name,
                    "endpoint_label": ep["label"],
                    "n": len(subset),
                    "n_events": int(subset[ep["event"]].sum()),
                    "cluster_hr": s.loc["cluster_z", "exp(coef)"],
                    "cluster_p": s.loc["cluster_z", "p"],
                    "mutation_hr": s.loc[gene, "exp(coef)"],
                    "mutation_p": s.loc[gene, "p"],
                    "interaction_hr": s.loc["interaction", "exp(coef)"],
                    "interaction_ci_lower": s.loc["interaction", "exp(coef) lower 95%"],
                    "interaction_ci_upper": s.loc["interaction", "exp(coef) upper 95%"],
                    "interaction_p": s.loc["interaction", "p"],
                })
            except Exception as e:
                print(f"    Error fitting interaction {gene} x {ep_name}: {e}")

    return pd.DataFrame(results)


def plot_km_interaction(
    df: pd.DataFrame,
    gene: str,
    ep_name: str,
    ep: dict,
    output_dir: Path,
):
    """Analysis 4: 4-group KM plot for {mutant/WT} × {high/low cluster}."""
    cols = [IMMUNE_FEATURE, gene, ep["time"], ep["event"]]
    subset = df[cols].dropna().copy()

    if len(subset) < 10:
        return

    median_cluster = subset[IMMUNE_FEATURE].median()
    subset["cluster_high"] = (subset[IMMUNE_FEATURE] >= median_cluster).astype(int)

    groups = {
        f"{gene}-mut / high cluster": (subset[gene] == 1) & (subset["cluster_high"] == 1),
        f"{gene}-mut / low cluster": (subset[gene] == 1) & (subset["cluster_high"] == 0),
        f"{gene}-WT / high cluster": (subset[gene] == 0) & (subset["cluster_high"] == 1),
        f"{gene}-WT / low cluster": (subset[gene] == 0) & (subset["cluster_high"] == 0),
    }

    colors = ["#b2182b", "#ef8a62", "#2166ac", "#67a9cf"]

    fig, ax = plt.subplots(figsize=(8, 6))

    group_data = []
    for (label, mask), color in zip(groups.items(), colors):
        g = subset[mask]
        if len(g) < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(g[ep["time"]], g[ep["event"]], label=f"{label} (n={len(g)})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color, linewidth=2)
        group_data.append(g)

    if len(group_data) < 2:
        plt.close()
        return

    # Multi-group log-rank test
    subset["group"] = np.nan
    for i, (label, mask) in enumerate(groups.items()):
        subset.loc[mask, "group"] = i
    subset_valid = subset.dropna(subset=["group"])

    try:
        lr = multivariate_logrank_test(
            subset_valid[ep["time"]],
            subset_valid["group"],
            subset_valid[ep["event"]],
        )
        p_val = lr.p_value
    except Exception:
        p_val = np.nan

    ax.set_title(
        f"{gene} mutation × immune clustering — {ep['label']}\n"
        f"Cluster median split = {median_cluster:.3f}, 4-group log-rank p = {p_val:.3f}",
        fontsize=12,
    )
    ax.set_xlabel("Days", fontsize=11)
    ax.set_ylabel("Survival probability", fontsize=11)
    ax.legend(loc="lower left", fontsize=9)

    plt.tight_layout()
    plt.savefig(
        output_dir / f"km_interaction_{gene}_{ep_name}.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close()
    print(f"    Saved: km_interaction_{gene}_{ep_name}.png")


def run_mutation_burden(df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Analysis 5: Mutation burden (count) as predictor."""
    available_genes = [g for g in genes if g in df.columns]
    df = df.copy()
    df["mutation_count"] = df[available_genes].sum(axis=1)

    results = []

    # Test mutation_count alone
    for ep_name, ep in ENDPOINTS.items():
        cols = ["mutation_count", ep["time"], ep["event"]]
        subset = df[cols].dropna().copy()

        if len(subset) < 10:
            continue

        try:
            cph = CoxPHFitter()
            cph.fit(subset, duration_col=ep["time"], event_col=ep["event"])
            s = cph.summary

            results.append({
                "model": "mutation_count_only",
                "endpoint": ep_name,
                "endpoint_label": ep["label"],
                "n": len(subset),
                "n_events": int(subset[ep["event"]].sum()),
                "mutation_count_hr": s.loc["mutation_count", "exp(coef)"],
                "mutation_count_ci_lower": s.loc["mutation_count", "exp(coef) lower 95%"],
                "mutation_count_ci_upper": s.loc["mutation_count", "exp(coef) upper 95%"],
                "mutation_count_p": s.loc["mutation_count", "p"],
                "cluster_hr": np.nan,
                "cluster_p": np.nan,
            })
        except Exception as e:
            print(f"    Error fitting mutation_count x {ep_name}: {e}")

    # Test mutation_count + db_tumor_cluster_fraction
    for ep_name, ep in ENDPOINTS.items():
        cols = ["mutation_count", IMMUNE_FEATURE, ep["time"], ep["event"]]
        subset = df[cols].dropna().copy()

        if len(subset) < 10:
            continue

        # Z-score the immune feature
        mean_val = subset[IMMUNE_FEATURE].mean()
        std_val = subset[IMMUNE_FEATURE].std()
        if std_val == 0:
            continue
        subset["cluster_z"] = (subset[IMMUNE_FEATURE] - mean_val) / std_val

        try:
            cph = CoxPHFitter()
            fit_cols = ["mutation_count", "cluster_z", ep["time"], ep["event"]]
            cph.fit(subset[fit_cols], duration_col=ep["time"], event_col=ep["event"])
            s = cph.summary

            results.append({
                "model": "mutation_count_plus_cluster",
                "endpoint": ep_name,
                "endpoint_label": ep["label"],
                "n": len(subset),
                "n_events": int(subset[ep["event"]].sum()),
                "mutation_count_hr": s.loc["mutation_count", "exp(coef)"],
                "mutation_count_ci_lower": s.loc["mutation_count", "exp(coef) lower 95%"],
                "mutation_count_ci_upper": s.loc["mutation_count", "exp(coef) upper 95%"],
                "mutation_count_p": s.loc["mutation_count", "p"],
                "cluster_hr": s.loc["cluster_z", "exp(coef)"],
                "cluster_p": s.loc["cluster_z", "p"],
            })
        except Exception as e:
            print(f"    Error fitting mutation_count + cluster x {ep_name}: {e}")

    return pd.DataFrame(results)


def plot_mutation_forest(uni_results: pd.DataFrame, output_path: Path):
    """Forest plot of univariate mutation HRs across endpoints."""
    if len(uni_results) == 0:
        return

    # Only plot endpoints that have data
    active_eps = [
        (ep_name, ep) for ep_name, ep in ENDPOINTS.items()
        if len(uni_results[uni_results["endpoint"] == ep_name]) > 0
    ]
    if not active_eps:
        return

    fig, axes = plt.subplots(1, len(active_eps), figsize=(7 * len(active_eps), 5))
    if len(active_eps) == 1:
        axes = [axes]

    for idx, (ep_name, ep) in enumerate(active_eps):
        ax = axes[idx]
        ep_data = uni_results[uni_results["endpoint"] == ep_name].copy()
        ep_data = ep_data.sort_values("gene").reset_index(drop=True)

        for i, (_, row) in enumerate(ep_data.iterrows()):
            color = "#b2182b" if row["p_value"] < 0.05 else "black"
            marker = "D" if row["p_value"] < 0.05 else "o"
            ax.errorbar(
                row["hr"], i,
                xerr=[[row["hr"] - row["ci_lower"]], [row["ci_upper"] - row["hr"]]],
                fmt=marker, color=color, capsize=4, markersize=7, linewidth=1.5,
            )

        ax.axvline(x=1, color="gray", linestyle="--", alpha=0.5)
        ax.set_xscale("log")

        labels = [
            f"{row['gene']} ({row['n_mutant']}/{row['n']}, p={row['p_value']:.3f})"
            for _, row in ep_data.iterrows()
        ]
        ax.set_yticks(range(len(ep_data)))
        ax.set_yticklabels(labels)
        ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=11)
        ax.set_title(ep["label"], fontsize=12)

    plt.suptitle("Univariate Cox — Mutation Status", fontsize=14)
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
    print("Mutation × Immune Clustering Interaction Analysis")
    print("=" * 70)

    # --- Data loading ---
    print("\nLoading data...")
    surv = pd.read_csv(data_dir / "survival_analysis_dataset.csv")
    print(f"  Survival dataset: {len(surv)} participants")

    mut_wide = load_mutation_data(project_root)
    all_genes = [c for c in mut_wide.columns if c != "participant_id"]
    print(f"  Sequencing data: {len(mut_wide)} participants, {len(all_genes)} genes")
    print(f"  Genes: {all_genes}")

    therapy = load_therapy_data(project_root)
    print(f"  Therapy data: {len(therapy)} participants")

    # Merge
    df = surv.merge(mut_wide, on="participant_id", how="inner")
    df = df.merge(therapy, on="participant_id", how="left")
    print(f"  Merged dataset: {len(df)} participants")

    # Mutation prevalences
    print("\n  Mutation prevalences:")
    for gene in all_genes:
        if gene in df.columns:
            n_mut = df[gene].sum()
            pct = 100 * n_mut / len(df)
            print(f"    {gene:10s}: {n_mut:.0f}/{len(df)} ({pct:.1f}%)")

    # --- Analysis 1: Univariate mutation → survival ---
    print("\n" + "-" * 70)
    print("Analysis 1: Univariate mutation → survival")
    print("-" * 70)
    uni_results = run_univariate_mutation_cox(df, all_genes)
    uni_path = data_dir / "mutation_univariate_cox.csv"
    uni_results.to_csv(uni_path, index=False)
    print(f"  Saved: {uni_path}")

    if len(uni_results) > 0:
        print("\n  Results:")
        for _, row in uni_results.iterrows():
            sig = "*" if row["p_value"] < 0.05 else " "
            print(
                f"  {sig} {row['gene']:10s} {row['endpoint']:3s}  "
                f"HR={row['hr']:.3f} ({row['ci_lower']:.3f}–{row['ci_upper']:.3f})  "
                f"p={row['p_value']:.3f}  "
                f"mut={row['n_mutant']}/{row['n']}"
            )
    else:
        print("  No genes met the ≥5 events criterion.")

    # Forest plot
    print("\n  Plotting mutation forest plot...")
    plot_mutation_forest(uni_results, figures_dir / "mutation_forest_plot.png")

    # --- Analysis 2: Stratified immune clustering ---
    print("\n" + "-" * 70)
    print("Analysis 2: Stratified immune clustering within mutation subgroups")
    print("-" * 70)
    strat_results = run_stratified_cox(df, FOCUS_GENES)
    strat_path = data_dir / "mutation_stratified_cox.csv"
    strat_results.to_csv(strat_path, index=False)
    print(f"  Saved: {strat_path}")

    if len(strat_results) > 0:
        print(f"\n  {IMMUNE_FEATURE} HR (per SD) within mutation strata:")
        for _, row in strat_results.iterrows():
            sig = "*" if row["p_value"] < 0.05 else " "
            print(
                f"  {sig} {row['gene']:6s} {row['stratum']:7s} {row['endpoint']:3s}  "
                f"HR={row['hr']:.3f} ({row['ci_lower']:.3f}–{row['ci_upper']:.3f})  "
                f"p={row['p_value']:.3f}  n={row['n']}"
            )

    # --- Analysis 3: Formal interaction test ---
    print("\n" + "-" * 70)
    print("Analysis 3: Formal interaction test (cluster × mutation)")
    print("-" * 70)
    interaction_results = run_interaction_test(df, FOCUS_GENES)
    interaction_path = data_dir / "mutation_interaction_cox.csv"
    interaction_results.to_csv(interaction_path, index=False)
    print(f"  Saved: {interaction_path}")

    if len(interaction_results) > 0:
        print("\n  Interaction term results:")
        for _, row in interaction_results.iterrows():
            sig = "*" if row["interaction_p"] < 0.05 else "~" if row["interaction_p"] < 0.2 else " "
            print(
                f"  {sig} {row['gene']:6s} {row['endpoint']:3s}  "
                f"interaction HR={row['interaction_hr']:.3f} "
                f"({row['interaction_ci_lower']:.3f}–{row['interaction_ci_upper']:.3f})  "
                f"p={row['interaction_p']:.3f}  "
                f"| cluster p={row['cluster_p']:.3f}, mutation p={row['mutation_p']:.3f}"
            )

    # --- Analysis 4: 4-group KM plots ---
    print("\n" + "-" * 70)
    print("Analysis 4: 4-group Kaplan-Meier plots (interaction p < 0.2)")
    print("-" * 70)

    if len(interaction_results) > 0:
        sig_interactions = interaction_results[interaction_results["interaction_p"] < 0.2]
        if len(sig_interactions) > 0:
            for _, row in sig_interactions.iterrows():
                plot_km_interaction(df, row["gene"], row["endpoint"], ENDPOINTS[row["endpoint"]], figures_dir)
        else:
            print("  No interactions with p < 0.2; plotting all focus gene × PFS combinations.")
            for gene in FOCUS_GENES:
                if gene in df.columns:
                    plot_km_interaction(df, gene, "pfs", ENDPOINTS["pfs"], figures_dir)
    else:
        print("  No interaction results to plot.")

    # --- Analysis 5: Mutation burden ---
    print("\n" + "-" * 70)
    print("Analysis 5: Mutation burden analysis")
    print("-" * 70)

    burden_results = run_mutation_burden(df, all_genes)
    if len(burden_results) > 0:
        burden_path = data_dir / "mutation_burden_cox.csv"
        burden_results.to_csv(burden_path, index=False)
        print(f"  Saved: {burden_path}")

        available_genes = [g for g in all_genes if g in df.columns]
        df_tmp = df.copy()
        df_tmp["mutation_count"] = df_tmp[available_genes].sum(axis=1)
        print(f"\n  Mutation count distribution: "
              f"mean={df_tmp['mutation_count'].mean():.1f}, "
              f"median={df_tmp['mutation_count'].median():.0f}, "
              f"range={df_tmp['mutation_count'].min():.0f}–{df_tmp['mutation_count'].max():.0f}")

        print("\n  Results:")
        for _, row in burden_results.iterrows():
            sig = "*" if row["mutation_count_p"] < 0.05 else " "
            cluster_str = ""
            if pd.notna(row["cluster_p"]):
                cluster_str = f"  | cluster HR={row['cluster_hr']:.3f} p={row['cluster_p']:.3f}"
            print(
                f"  {sig} {row['model']:30s} {row['endpoint']:3s}  "
                f"count HR={row['mutation_count_hr']:.3f} "
                f"({row['mutation_count_ci_lower']:.3f}–{row['mutation_count_ci_upper']:.3f})  "
                f"p={row['mutation_count_p']:.3f}"
                f"{cluster_str}"
            )

    # --- Summary ---
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    if len(interaction_results) > 0:
        sig_any = interaction_results[interaction_results["interaction_p"] < 0.05]
        trend_any = interaction_results[
            (interaction_results["interaction_p"] >= 0.05)
            & (interaction_results["interaction_p"] < 0.2)
        ]
        print(f"  Significant interactions (p < 0.05): {len(sig_any)}")
        if len(sig_any) > 0:
            for _, row in sig_any.iterrows():
                print(f"    {row['gene']} × {row['endpoint']}: interaction p = {row['interaction_p']:.4f}")
        print(f"  Trending interactions (0.05 ≤ p < 0.2): {len(trend_any)}")
        if len(trend_any) > 0:
            for _, row in trend_any.iterrows():
                print(f"    {row['gene']} × {row['endpoint']}: interaction p = {row['interaction_p']:.3f}")

    print("\nDone!")


if __name__ == "__main__":
    main()
