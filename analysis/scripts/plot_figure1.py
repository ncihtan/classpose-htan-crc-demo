#!/usr/bin/env python3
"""
Figure 1: Multi-panel overview of HTAN CRC ClassPose cell segmentation analysis.

Outputs separate files for each panel:
- figure1a_cohort.png - Clinical metadata heatmap
- figure1b_validation.png - Model validation scatter plots (including total cells)
- figure1c_heterogeneity.png - 100% stacked bar charts (CoNIC and ConSeP side by side)
- figure1d_serial_sections.png - Technical replicate comparison (HTA7_919)
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
from scipy import stats

from load_cell_densities import (
    get_project_root,
    load_cell_densities,
    pivot_for_stacked_bar,
    get_total_cells_per_sample,
    get_sample_metadata,
    CONSEP_IMMUNE_CLASSES,
    CONIC_IMMUNE_CLASSES,
)

# =============================================================================
# Color Palettes - ColorBrewer qualitative palettes (colorblind-friendly)
# =============================================================================
# CoNIC: Using "Dark2" palette from ColorBrewer (8-class qualitative)
# Good for categorical data, distinguishable, colorblind-safe
CONIC_COLORS = {
    "Epithelial": "#1b9e77",    # Teal
    "Connective": "#d95f02",    # Orange
    "Lymphocyte": "#7570b3",    # Purple
    "Plasma cell": "#e7298a",   # Pink
    "Neutrophil": "#66a61e",    # Green
    "Eosinophil": "#e6ab02",    # Gold
}

# ConSeP: Using "Set1" palette from ColorBrewer (9-class qualitative)
# Vibrant, highly distinguishable colors
CONSEP_COLORS = {
    "Malignant epithelial": "#e41a1c",  # Red
    "Healthy epithelial": "#377eb8",    # Blue
    "Stroma": "#4daf4a",                # Green
    "Muscle": "#984ea3",                # Purple
    "Inflammatory": "#ff7f00",          # Orange
    "Other": "#999999",                 # Gray
}

# Cell class mappings for model comparison
CONIC_IMMUNE = ["Lymphocyte", "Plasma cell", "Neutrophil", "Eosinophil"]
CONSEP_IMMUNE = ["Inflammatory"]
CONIC_EPITHELIAL = ["Epithelial"]
CONSEP_EPITHELIAL = ["Malignant epithelial", "Healthy epithelial"]
CONIC_STROMAL = ["Connective"]
CONSEP_STROMAL = ["Stroma"]


def get_category_counts(df, categories, by="sample_id"):
    """Sum counts for a list of cell class categories."""
    pivoted = df.pivot_table(
        index=by, columns="cell_class", values="count", aggfunc="sum"
    ).fillna(0)
    valid_cats = [c for c in categories if c in pivoted.columns]
    if valid_cats:
        return pivoted[valid_cats].sum(axis=1)
    return pd.Series(0, index=pivoted.index)


def load_clinical_metadata():
    """Load and process clinical metadata."""
    project_root = get_project_root()
    part_df = pd.read_csv(project_root / "metadata" / "participant.tsv", sep="\t")
    part_df = part_df.rename(columns={
        "HTAN Participant ID": "participant_id",
        "Age at Diagnosis (years)": "age",
        "Gender": "gender",
        "Site of Resection or Biopsy": "site",
        "AJCC Pathologic Stage": "stage",
        "Tumor Grade": "grade",
        "Vital Status": "vital_status",
        "Progression or Recurrence": "progression",
    })
    return part_df


def get_sample_order_by_participant(df):
    """Get sample order sorted by participant ID."""
    meta = get_sample_metadata(df)
    # Sort by participant_id, then by sample_id for consistent ordering
    meta_sorted = meta.sort_values(["participant_id", "sample_id"])
    return meta_sorted["sample_id"].tolist()


def load_drs_mapping():
    """Load mapping from sample_id (dg_4DFC_UUID) to Data File ID."""
    project_root = get_project_root()
    drs_df = pd.read_csv(project_root / "metadata" / "files_drs.csv")

    # Extract UUID from drs_uri and convert to sample_id format
    # drs://nci-crdc.datacommons.io/dg.4DFC/UUID -> dg_4DFC_UUID
    mapping = {}
    for _, row in drs_df.iterrows():
        drs_uri = row["drs_uri"]
        # Extract UUID from URI
        uuid = drs_uri.split("/")[-1]
        sample_id = f"dg_4DFC_{uuid}"
        mapping[sample_id] = row["data_file_id"]

    return mapping


def create_sample_display_labels(df):
    """Create display labels using Data File ID from metadata."""
    meta = get_sample_metadata(df)
    drs_mapping = load_drs_mapping()

    labels = {}
    for _, row in meta.iterrows():
        sample_id = row["sample_id"]
        # Use Data File ID from metadata, fall back to sample_id if not found
        labels[sample_id] = drs_mapping.get(sample_id, sample_id)

    return labels


def plot_figure1a_cohort(consep_df, sample_order, sample_labels, output_dir):
    """Panel A: Clinical cohort overview as colored annotation bars."""

    clinical_df = load_clinical_metadata()
    meta = get_sample_metadata(consep_df)
    sample_to_participant = dict(zip(meta["sample_id"], meta["participant_id"]))

    # Map sample order to participants
    participants_in_order = [sample_to_participant[s] for s in sample_order]

    # Merge clinical data
    clinical = clinical_df.set_index("participant_id")

    # Define variables and their color mappings
    variables = {
        "Age": {
            "func": lambda x: "<50" if x < 50 else ("50-59" if x < 60 else ("60-69" if x < 70 else "70+")),
            "colors": {"<50": "#ffffb2", "50-59": "#fecc5c", "60-69": "#fd8d3c", "70+": "#e31a1c"},
            "col": "age"
        },
        "Sex": {
            "func": lambda x: str(x).capitalize() if pd.notna(x) else "Unknown",
            "colors": {"Male": "#6baed6", "Female": "#fc9272", "Unknown": "#d9d9d9"},
            "col": "gender"
        },
        "Site": {
            "func": lambda x: "Rectum" if "rect" in str(x).lower() else (
                "R.Colon" if any(s in str(x).lower() for s in ["ascending", "cecum"]) else (
                "L.Colon" if "descending" in str(x).lower() else (
                "Sigmoid" if "sigmoid" in str(x).lower() else (
                "T.Colon" if "transverse" in str(x).lower() else "Other")))),
            "colors": {"R.Colon": "#1b9e77", "T.Colon": "#d95f02", "L.Colon": "#7570b3",
                      "Sigmoid": "#e7298a", "Rectum": "#66a61e", "Other": "#999999"},
            "col": "site"
        },
        "Stage": {
            "func": lambda x: "I" if str(x).startswith("Stage I") and "II" not in str(x) and "IV" not in str(x) else (
                "II" if "II" in str(x) and "IV" not in str(x) else (
                "III" if "III" in str(x) else (
                "IV" if "IV" in str(x) else "Unk"))),
            "colors": {"I": "#eff3ff", "II": "#bdd7e7", "III": "#6baed6", "IV": "#2171b5", "Unk": "#d9d9d9"},
            "col": "stage"
        },
        "Grade": {
            "func": lambda x: "Low" if "Low" in str(x) else ("High" if "High" in str(x) else "Unk"),
            "colors": {"Low": "#74c476", "High": "#e31a1c", "Unk": "#d9d9d9"},
            "col": "grade"
        },
        "Vital": {
            "func": lambda x: "Alive" if str(x) == "Alive" else ("Dead" if str(x) == "Dead" else "Unk"),
            "colors": {"Alive": "#31a354", "Dead": "#de2d26", "Unk": "#d9d9d9"},
            "col": "vital_status"
        },
        "Recur": {
            "func": lambda x: "Yes" if "Yes" in str(x) else ("No" if str(x) == "no" else "Unk"),
            "colors": {"No": "#74c476", "Yes": "#e31a1c", "Unk": "#d9d9d9"},
            "col": "progression"
        },
    }

    fig, ax = plt.subplots(figsize=(16, 4))

    n_samples = len(sample_order)
    n_vars = len(variables)
    bar_height = 0.8

    for var_idx, (var_name, var_info) in enumerate(variables.items()):
        y_pos = n_vars - var_idx - 1

        for sample_idx, sample_id in enumerate(sample_order):
            participant = sample_to_participant[sample_id]

            if participant in clinical.index:
                raw_val = clinical.loc[participant, var_info["col"]]
                cat_val = var_info["func"](raw_val)
            else:
                cat_val = "Unk" if "Unk" in var_info["colors"] else "Unknown"

            color = var_info["colors"].get(cat_val, "#d9d9d9")
            ax.barh(y_pos, 1, left=sample_idx, height=bar_height, color=color, edgecolor="white", linewidth=0.5)

    # Y-axis labels
    ax.set_yticks(range(n_vars))
    ax.set_yticklabels(list(variables.keys())[::-1], fontsize=10)
    ax.set_ylim(-0.5, n_vars - 0.5)

    # X-axis - show participant IDs (shorter than Data File IDs)
    ax.set_xticks([i + 0.5 for i in range(n_samples)])
    ax.set_xticklabels([sample_to_participant[s] for s in sample_order], fontsize=7, rotation=90, ha="center")
    ax.set_xlim(0, n_samples)
    ax.set_xlabel("Participant ID", fontsize=10)

    ax.set_title("A. Cohort Overview", fontsize=14, fontweight="bold", loc="left")

    # Create legend
    legend_elements = []
    for var_name, var_info in list(variables.items())[:4]:  # First 4 variables
        for cat, color in var_info["colors"].items():
            if cat not in ["Unk", "Unknown"]:
                legend_elements.append(Patch(facecolor=color, label=f"{var_name}: {cat}"))

    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.01, 1),
              fontsize=7, ncol=2, title="Clinical Variables")

    n_participants = len(set(participants_in_order))
    ax.annotate(f"{n_samples} images, {n_participants} patients",
                xy=(0.99, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8))

    plt.tight_layout()

    output_path = output_dir / "figure1a_cohort.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(output_dir / "figure1a_cohort.pdf", bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_path}")
    plt.close()


def plot_figure1b_validation(conic_df, consep_df, output_dir):
    """Panel B: Model validation - 4 scatter plots including total cells."""

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    # Total cells
    conic_total = get_total_cells_per_sample(conic_df, by="sample_id")
    consep_total = get_total_cells_per_sample(consep_df, by="sample_id")
    common = conic_total.index.intersection(consep_total.index)

    # Category counts
    conic_immune = get_category_counts(conic_df, CONIC_IMMUNE, by="sample_id")
    consep_immune = get_category_counts(consep_df, CONSEP_IMMUNE, by="sample_id")
    conic_epithelial = get_category_counts(conic_df, CONIC_EPITHELIAL, by="sample_id")
    consep_epithelial = get_category_counts(consep_df, CONSEP_EPITHELIAL, by="sample_id")
    conic_stromal = get_category_counts(conic_df, CONIC_STROMAL, by="sample_id")
    consep_stromal = get_category_counts(consep_df, CONSEP_STROMAL, by="sample_id")

    comparisons = [
        ("Total Cells", conic_total.loc[common], consep_total.loc[common], "#666666"),
        ("Immune Cells", conic_immune.loc[common], consep_immune.loc[common], "#4daf4a"),
        ("Epithelial Cells", conic_epithelial.loc[common], consep_epithelial.loc[common], "#e41a1c"),
        ("Stromal Cells", conic_stromal.loc[common], consep_stromal.loc[common], "#377eb8"),
    ]

    for ax, (title, conic_vals, consep_vals, color) in zip(axes, comparisons):
        ax.scatter(conic_vals, consep_vals, s=50, alpha=0.7, c=color, edgecolors="white", linewidth=0.5)

        r, p_value = stats.pearsonr(conic_vals, consep_vals)

        # Find max for this comparison (same for x and y)
        plot_max = max(conic_vals.max(), consep_vals.max()) * 1.05

        # Identity line
        ax.plot([0, plot_max], [0, plot_max], "k--", alpha=0.4, linewidth=1, label="y=x")

        # Regression line
        slope, intercept, _, _, _ = stats.linregress(conic_vals, consep_vals)
        x_line = np.array([0, plot_max])
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, color=color, alpha=0.8, linewidth=2)

        ax.set_xlabel("CoNIC", fontsize=10)
        ax.set_ylabel("ConSeP", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")

        ax.annotate(f"r = {r:.2f}", xy=(0.05, 0.95), xycoords="axes fraction",
                   ha="left", va="top", fontsize=10, fontweight="bold",
                   bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="gray", alpha=0.8))

        # Set identical x and y limits for this plot (square aspect)
        ax.set_xlim(0, plot_max)
        ax.set_ylim(0, plot_max)

        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}K"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}K"))
        ax.tick_params(axis="both", labelsize=9)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle("B. Model Validation: CoNIC vs ConSeP Agreement", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    output_path = output_dir / "figure1b_validation.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(output_dir / "figure1b_validation.pdf", bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_path}")
    plt.close()


def plot_figure1c_heterogeneity(conic_df, consep_df, sample_order, sample_labels, output_dir):
    """Panel C: 100% stacked bar charts for both models side by side."""

    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    # CoNIC plot
    ax_conic = axes[0]
    pivoted_conic = pivot_for_stacked_bar(conic_df, by="sample_id")
    pivoted_conic = pivoted_conic.reindex(sample_order)
    pivoted_conic_pct = pivoted_conic.div(pivoted_conic.sum(axis=1), axis=0) * 100

    conic_col_order = ["Connective", "Epithelial", "Eosinophil", "Neutrophil", "Plasma cell", "Lymphocyte"]
    conic_col_order = [c for c in conic_col_order if c in pivoted_conic_pct.columns]
    pivoted_conic_pct = pivoted_conic_pct[conic_col_order]

    bottom = np.zeros(len(pivoted_conic_pct))
    x = np.arange(len(pivoted_conic_pct))

    for col in conic_col_order:
        color = CONIC_COLORS.get(col, "#999999")
        ax_conic.bar(x, pivoted_conic_pct[col].values, bottom=bottom, color=color,
                     edgecolor="white", linewidth=0.3, width=0.9, label=col)
        bottom += pivoted_conic_pct[col].values

    ax_conic.set_ylabel("Proportion (%)", fontsize=11)
    ax_conic.set_title("CoNIC Model", fontsize=12, fontweight="bold")
    ax_conic.set_ylim(0, 100)
    ax_conic.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax_conic.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax_conic.set_axisbelow(True)

    # ConSeP plot
    ax_consep = axes[1]
    pivoted_consep = pivot_for_stacked_bar(consep_df, by="sample_id")
    pivoted_consep = pivoted_consep.reindex(sample_order)
    pivoted_consep_pct = pivoted_consep.div(pivoted_consep.sum(axis=1), axis=0) * 100

    consep_col_order = ["Other", "Muscle", "Stroma", "Healthy epithelial", "Malignant epithelial", "Inflammatory"]
    consep_col_order = [c for c in consep_col_order if c in pivoted_consep_pct.columns]
    pivoted_consep_pct = pivoted_consep_pct[consep_col_order]

    bottom = np.zeros(len(pivoted_consep_pct))

    for col in consep_col_order:
        color = CONSEP_COLORS.get(col, "#999999")
        ax_consep.bar(x, pivoted_consep_pct[col].values, bottom=bottom, color=color,
                      edgecolor="white", linewidth=0.3, width=0.9, label=col)
        bottom += pivoted_consep_pct[col].values

    ax_consep.set_ylabel("Proportion (%)", fontsize=11)
    ax_consep.set_title("ConSeP Model", fontsize=12, fontweight="bold")
    ax_consep.set_ylim(0, 100)
    ax_consep.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax_consep.yaxis.grid(True, linestyle="--", alpha=0.3)
    ax_consep.set_axisbelow(True)

    # X-axis labels (file IDs) - show all
    ax_consep.set_xticks(x)
    ax_consep.set_xticklabels([sample_labels[s] for s in sample_order], fontsize=8, rotation=90, ha="center")
    ax_consep.set_xlim(-0.5, len(sample_order) - 0.5)
    ax_consep.set_xlabel("Sample (ordered by participant ID)", fontsize=11)

    fig.suptitle("C. Tumor Microenvironment Heterogeneity", fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout()

    output_path = output_dir / "figure1c_heterogeneity.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(output_dir / "figure1c_heterogeneity.pdf", bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_path}")
    plt.close()


def plot_figure1d_serial_sections(conic_df, consep_df, output_dir):
    """Panel D: Technical replicate comparison for HTA7_919 serial sections."""

    meta = get_sample_metadata(consep_df)
    hta919_samples = meta[meta["participant_id"] == "HTA7_919"]["sample_id"].tolist()

    if len(hta919_samples) != 2:
        print(f"Warning: Expected 2 samples for HTA7_919, found {len(hta919_samples)}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    sample1, sample2 = sorted(hta919_samples)

    # ConSeP comparison
    consep_data = {}
    for sample_id in hta919_samples:
        sample_df = consep_df[consep_df["sample_id"] == sample_id]
        total = sample_df["count"].sum()
        props = sample_df.set_index("cell_class")["count"] / total * 100
        consep_data[sample_id] = props

    # CoNIC comparison
    conic_data = {}
    for sample_id in hta919_samples:
        sample_df = conic_df[conic_df["sample_id"] == sample_id]
        total = sample_df["count"].sum()
        props = sample_df.set_index("cell_class")["count"] / total * 100
        conic_data[sample_id] = props

    # ConSeP bar chart
    ax1 = axes[0, 0]
    cell_types_consep = list(CONSEP_COLORS.keys())
    x = np.arange(len(cell_types_consep))
    width = 0.35

    vals1 = [consep_data[sample1].get(ct, 0) for ct in cell_types_consep]
    vals2 = [consep_data[sample2].get(ct, 0) for ct in cell_types_consep]

    ax1.bar(x - width/2, vals1, width, label="Section 1", color="#1f77b4", alpha=0.8)
    ax1.bar(x + width/2, vals2, width, label="Section 2", color="#ff7f0e", alpha=0.8)

    ax1.set_ylabel("Proportion (%)", fontsize=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels([ct.split()[0] for ct in cell_types_consep], rotation=45, ha="right", fontsize=9)
    ax1.legend(fontsize=9)
    ax1.set_title("ConSeP: Cell Proportions", fontsize=11, fontweight="bold")
    ax1.yaxis.grid(True, linestyle="--", alpha=0.3)

    # ConSeP scatter
    ax2 = axes[0, 1]
    ax2.scatter(vals1, vals2, s=100, c=[CONSEP_COLORS[ct] for ct in cell_types_consep],
                edgecolors="black", linewidth=1)
    for i, ct in enumerate(cell_types_consep):
        ax2.annotate(ct.split()[0], (vals1[i], vals2[i]), fontsize=7, xytext=(3, 3), textcoords="offset points")

    r_consep, _ = stats.pearsonr(vals1, vals2)
    max_val = max(max(vals1), max(vals2)) * 1.1
    ax2.plot([0, max_val], [0, max_val], "k--", alpha=0.5)
    ax2.set_xlim(0, max_val)
    ax2.set_ylim(0, max_val)
    ax2.set_xlabel("Section 1 (%)", fontsize=10)
    ax2.set_ylabel("Section 2 (%)", fontsize=10)
    ax2.set_title(f"ConSeP: r = {r_consep:.3f}", fontsize=11, fontweight="bold")
    ax2.set_aspect("equal")
    ax2.grid(True, linestyle="--", alpha=0.3)

    # CoNIC bar chart
    ax3 = axes[1, 0]
    cell_types_conic = list(CONIC_COLORS.keys())
    x = np.arange(len(cell_types_conic))

    vals1_conic = [conic_data[sample1].get(ct, 0) for ct in cell_types_conic]
    vals2_conic = [conic_data[sample2].get(ct, 0) for ct in cell_types_conic]

    ax3.bar(x - width/2, vals1_conic, width, label="Section 1", color="#1f77b4", alpha=0.8)
    ax3.bar(x + width/2, vals2_conic, width, label="Section 2", color="#ff7f0e", alpha=0.8)

    ax3.set_ylabel("Proportion (%)", fontsize=10)
    ax3.set_xticks(x)
    ax3.set_xticklabels(cell_types_conic, rotation=45, ha="right", fontsize=9)
    ax3.legend(fontsize=9)
    ax3.set_title("CoNIC: Cell Proportions", fontsize=11, fontweight="bold")
    ax3.yaxis.grid(True, linestyle="--", alpha=0.3)

    # CoNIC scatter
    ax4 = axes[1, 1]
    ax4.scatter(vals1_conic, vals2_conic, s=100, c=[CONIC_COLORS[ct] for ct in cell_types_conic],
                edgecolors="black", linewidth=1)
    for i, ct in enumerate(cell_types_conic):
        ax4.annotate(ct, (vals1_conic[i], vals2_conic[i]), fontsize=7, xytext=(3, 3), textcoords="offset points")

    r_conic, _ = stats.pearsonr(vals1_conic, vals2_conic)
    max_val = max(max(vals1_conic), max(vals2_conic)) * 1.1
    ax4.plot([0, max_val], [0, max_val], "k--", alpha=0.5)
    ax4.set_xlim(0, max_val)
    ax4.set_ylim(0, max_val)
    ax4.set_xlabel("Section 1 (%)", fontsize=10)
    ax4.set_ylabel("Section 2 (%)", fontsize=10)
    ax4.set_title(f"CoNIC: r = {r_conic:.3f}", fontsize=11, fontweight="bold")
    ax4.set_aspect("equal")
    ax4.grid(True, linestyle="--", alpha=0.3)

    fig.suptitle("D. Technical Replicate: HTA7_919 Serial Sections", fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout()

    output_path = output_dir / "figure1d_serial_sections.png"
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(output_dir / "figure1d_serial_sections.pdf", bbox_inches="tight", facecolor="white")
    print(f"Saved: {output_path}")
    plt.close()


def create_figure1():
    """Create all Figure 1 panels as separate files."""
    print("Loading data...")
    conic_df = load_cell_densities("conic")
    consep_df = load_cell_densities("consep")

    meta = get_sample_metadata(consep_df)
    n_files = meta["sample_id"].nunique()
    n_participants = meta["participant_id"].nunique()

    # Order by participant ID
    sample_order = get_sample_order_by_participant(consep_df)
    sample_labels = create_sample_display_labels(consep_df)

    print(f"Loaded {n_files} H&E images from {n_participants} participants")
    print(f"Total ConSeP cells: {consep_df['count'].sum():,}")
    print(f"Total CoNIC cells: {conic_df['count'].sum():,}")

    output_dir = get_project_root() / "figures"
    output_dir.mkdir(exist_ok=True)

    print("\nGenerating Panel A (Cohort Overview)...")
    plot_figure1a_cohort(consep_df, sample_order, sample_labels, output_dir)

    print("Generating Panel B (Model Validation)...")
    plot_figure1b_validation(conic_df, consep_df, output_dir)

    print("Generating Panel C (TME Heterogeneity - both models)...")
    plot_figure1c_heterogeneity(conic_df, consep_df, sample_order, sample_labels, output_dir)

    print("Generating Panel D (Serial Sections)...")
    plot_figure1d_serial_sections(conic_df, consep_df, output_dir)

    print("\n" + "=" * 60)
    print("Figure 1 Summary")
    print("=" * 60)
    print(f"H&E images: {n_files} (from {n_participants} participants)")
    print(f"Total cells (ConSeP): {consep_df['count'].sum():,}")
    print(f"Total cells (CoNIC): {conic_df['count'].sum():,}")
    print("\nOutput files:")
    print("  - figure1a_cohort.png/pdf")
    print("  - figure1b_validation.png/pdf")
    print("  - figure1c_heterogeneity.png/pdf")
    print("  - figure1d_serial_sections.png/pdf")


if __name__ == "__main__":
    create_figure1()
