#!/usr/bin/env python3
"""
Boxplot of Pearson r across conditions for each method-pair.

Each box = one method pair
Points inside each box = Pearson r from:
  0_vs_0, 0_vs_1, 1_vs_0, 1_vs_1

Updated for 4 methods:
- concatenated
- consensus_avg
- ancestral_avg
- PAO1PA14_ref
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ==========================
# HARD-CODED PATHS
# ==========================
SUMMARY_CSV = (
    "/storage/home/hcoda1/2/bli629/scratch/three_method_comparison/regression/out_log_v3/regression_summary.csv"
)

OUTDIR = (
    "/storage/home/hcoda1/2/bli629/scratch/three_method_comparison/regression/out_log_v3"
)

METRIC = "pearson_r"   # change to "spearman_rho" if needed
# ==========================


def clean_method_name(s):
    """Make method names shorter/cleaner for plotting."""
    s = s.replace("_avg", "")
    s = s.replace("PAO1PA14_ref", "PAO1PA14")
    return s


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    df = pd.read_csv(SUMMARY_CSV)

    # Ensure numeric
    df[METRIC] = pd.to_numeric(df[METRIC], errors="coerce")

    # Drop skipped rows
    use = df.dropna(subset=[METRIC]).copy()
    if use.empty:
        raise RuntimeError(f"No valid {METRIC} values found.")

    # Clean pair labels
    use["x_clean"] = use["x_method"].map(clean_method_name)
    use["y_clean"] = use["y_method"].map(clean_method_name)
    use["pair"] = use["x_clean"] + " vs " + use["y_clean"]

    # Canonical order for 4 methods = 6 pairs
    pair_order = [
        "concatenated vs consensus",
        "concatenated vs ancestral",
        "concatenated vs PAO1PA14",
        "consensus vs ancestral",
        "consensus vs PAO1PA14",
        "ancestral vs PAO1PA14",
    ]

    data = []
    labels = []

    for p in pair_order:
        vals = use.loc[use["pair"] == p, METRIC].to_numpy()
        if len(vals) > 0:
            data.append(vals)
            labels.append(p)

    if len(data) == 0:
        raise RuntimeError("No data collected for boxplot.")

    # Save raw values used
    points_csv = os.path.join(OUTDIR, f"{METRIC}_points.csv")
    use.sort_values(["pair", "condition"]).to_csv(points_csv, index=False)

    # ---- Plot ----
    colors = [
        "#4C72B0",  # blue
        "#DD8452",  # orange
        "#55A868",  # green
        "#C44E52",  # red
        "#8172B3",  # purple
        "#937860",  # brown
    ]

    plt.figure(figsize=(10, 5))

    bp = plt.boxplot(
        data,
        patch_artist=True,
        showmeans=True,
    )

    # Color boxes
    for patch, color in zip(bp["boxes"], colors[:len(data)]):
        patch.set_facecolor(color)

    # Style medians and means
    for median in bp["medians"]:
        median.set_color("black")

    for mean in bp["means"]:
        mean.set_markerfacecolor("black")
        mean.set_markeredgecolor("black")

    # Remove x ticks (legend explains)
    plt.xticks([])

    ylabel = "Pearson r" if METRIC == "pearson_r" else "Spearman rho"
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} across conditions for each method pair")

    # Legend
    legend_handles = [
        Patch(facecolor=c, label=l)
        for c, l in zip(colors[:len(labels)], labels)
    ]
    plt.legend(
        handles=legend_handles,
        loc="best",
        frameon=False,
        fontsize=9,
        handlelength=1.0,
        handleheight=0.7,
        labelspacing=0.3,
    )

    plt.tight_layout()

    # SAVE FIGURE
    out_png = os.path.join(OUTDIR, f"{METRIC}_boxplot.png")
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"[OK] Wrote boxplot: {out_png}")
    print(f"[OK] Wrote values : {points_csv}")


if __name__ == "__main__":
    main()