#!/usr/bin/env python3
"""
Parse Nextflow trace files and plot ClassPose (CLASSPOSE_PREDICT_WSI) runtimes.

Three panels:
  1. Raw wall-clock duration (hours)
  2. Duration normalised by total cell count (seconds / 1000 cells)
  3. Duration normalised by input data size from trace rchar (seconds / GB read)

Cell counts come from classpose_outputs/data/<patient>/<slide>/{conic,consep}/*_cell_densities.csv.
Input size comes from the rchar column in the trace (bytes read during the task).

Usage:
  uv run python scripts/plot_classpose_runtimes.py
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRACE_DIR = PROJECT_ROOT.parent / "pipeline" / "reports"
OUTPUTS_DIR = PROJECT_ROOT / "classpose_outputs" / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_duration(s: str) -> float:
    """Convert Nextflow duration string to seconds. E.g. '1h 27m 30s' → 5250.0"""
    s = s.strip()
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)\s*([hms])", s):
        v = float(value)
        if unit == "h":
            total += v * 3600
        elif unit == "m":
            total += v * 60
        elif unit == "s":
            total += v
    return total


def parse_bytes(s: str) -> float:
    """Convert Nextflow size string (e.g. '44.9 GB', '512 MB') to bytes."""
    s = s.strip()
    m = re.match(r"([\d.]+)\s*([KMGT]?B)", s)
    if not m:
        return float("nan")
    value, unit = float(m.group(1)), m.group(2)
    multipliers = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}
    return value * multipliers.get(unit, 1)


# ── Parse trace ───────────────────────────────────────────────────────────────

def load_trace() -> pd.DataFrame:
    """Load the most recent trace file and return CLASSPOSE_PREDICT_WSI rows."""
    trace_files = sorted(TRACE_DIR.glob("execution_trace_*.txt"))
    if not trace_files:
        raise FileNotFoundError(f"No trace files found in {TRACE_DIR}")
    trace_path = trace_files[-1]
    print(f"Using trace: {trace_path.name}")

    df = pd.read_csv(trace_path, sep="\t")

    cp = df[df["name"].str.startswith("CLASSPOSE_PREDICT_WSI")].copy()
    print(f"  {len(cp)} CLASSPOSE_PREDICT_WSI tasks")

    # Extract sample_id and model from name like "CLASSPOSE_PREDICT_WSI (dg_4DFC_...-conic)"
    def parse_name(name):
        m = re.search(r"\((.+)-(conic|consep)\)", name)
        if m:
            return m.group(1), m.group(2)
        return None, None

    cp[["sample_id", "model"]] = cp["name"].apply(
        lambda n: pd.Series(parse_name(n))
    )

    cp["duration_s"] = cp["duration"].apply(parse_duration)
    cp["realtime_s"] = cp["realtime"].apply(parse_duration)
    cp["rchar_bytes"] = cp["rchar"].apply(parse_bytes)

    return cp.reset_index(drop=True)


# ── Cell counts ───────────────────────────────────────────────────────────────

def load_cell_counts() -> dict[tuple[str, str], int]:
    """Return {(sample_id, model): total_cell_count} from density CSVs."""
    counts = {}
    for patient_dir in OUTPUTS_DIR.iterdir():
        if not patient_dir.is_dir():
            continue
        for slide_dir in patient_dir.iterdir():
            for model in ("conic", "consep"):
                model_dir = slide_dir / model
                if not model_dir.exists():
                    continue
                csvs = list(model_dir.glob("*_cell_densities.csv"))
                if not csvs:
                    continue
                # sample_id is the dg_4DFC_... part of the filename
                fname = csvs[0].stem  # e.g. dg_4DFC_..._consep_cell_densities
                sample_id = re.sub(f"_{model}_cell_densities$", "", fname)
                df = pd.read_csv(csvs[0])
                total = df.loc[df["region"] == "tissue", "count"].sum()
                counts[(sample_id, model)] = int(total)
    print(f"Loaded cell counts for {len(counts)} (sample, model) pairs")
    return counts


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_runtime_histograms(cp: pd.DataFrame, counts: dict) -> None:
    # Attach cell counts
    cp["cell_count"] = cp.apply(
        lambda r: counts.get((r["sample_id"], r["model"]), np.nan), axis=1
    )

    # Normalised metrics
    cp["s_per_1k_cells"] = cp["realtime_s"] / (cp["cell_count"] / 1000)
    cp["s_per_gb_rchar"] = cp["realtime_s"] / (cp["rchar_bytes"] / 1e9)

    models = ["conic", "consep"]
    colors = {"conic": "#3498db", "consep": "#e74c3c"}
    positions = {"conic": 0, "consep": 1}

    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    panels = [
        ("realtime_s", "Realtime duration (h)", lambda x: x / 3600, "Raw runtime (realtime)"),
        ("s_per_1k_cells", "Realtime per 1,000 cells (s)", lambda x: x, "Normalised by cell count"),
        ("s_per_gb_rchar", "Realtime per GB read (s)", lambda x: x, "Normalised by input data (rchar)"),
    ]

    # x positions: box_left, dots_conic, dots_consep, box_right
    xpos = {"box_conic": 0, "conic": 1, "consep": 2, "box_consep": 3}
    box_width = 0.35

    for ax, (col, ylabel, transform, title) in zip(axes, panels):
        wide = cp.pivot(index="sample_id", columns="model", values=col).dropna()

        # Connecting lines
        for _, row in wide.iterrows():
            ax.plot([xpos["conic"], xpos["consep"]],
                    [transform(pd.Series([row["conic"]])).iloc[0],
                     transform(pd.Series([row["consep"]])).iloc[0]],
                    color="grey", alpha=0.3, linewidth=0.8, zorder=1)

        # Dots
        for model in models:
            vals = transform(wide[model])
            ax.scatter([xpos[model]] * len(vals), vals,
                       color=colors[model], s=30, zorder=2, alpha=0.8)

        # Boxes on the outside
        for model, box_x in [("conic", xpos["box_conic"]), ("consep", xpos["box_consep"])]:
            vals = transform(wide[model]).values
            bp = ax.boxplot(vals, positions=[box_x], widths=box_width,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", linewidth=2))
            bp["boxes"][0].set_facecolor(colors[model])
            bp["boxes"][0].set_alpha(0.75)

        ax.set_xticks([xpos["box_conic"], xpos["conic"], xpos["consep"], xpos["box_consep"]])
        ax.set_xticklabels(["", "conic", "consep", ""], fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlim(-0.5, 3.5)

    plt.suptitle("CLASSPOSE_PREDICT_WSI runtimes", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = FIGURES_DIR / "classpose_runtimes.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out.name}")

    # Summary table
    summary = cp.groupby("model").agg(
        n=("realtime_s", "count"),
        median_h=("realtime_s", lambda x: np.median(x) / 3600),
        mean_h=("realtime_s", lambda x: x.mean() / 3600),
        max_h=("realtime_s", lambda x: x.max() / 3600),
        median_s_per_1k=("s_per_1k_cells", "median"),
        median_s_per_gb=("s_per_gb_rchar", "median"),
    ).round(2)
    print("\nSummary:")
    print(summary.to_string())


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    cp = load_trace()
    counts = load_cell_counts()
    plot_runtime_histograms(cp, counts)


if __name__ == "__main__":
    main()
