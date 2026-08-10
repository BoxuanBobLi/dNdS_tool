#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Input

base_dir = "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output"
out_png = "average_dnds_ref_comparison.png"
out_csv = "average_dnds_ref_comparison_summary.csv"

# 0 = Clade B, 1 = Clade A
conditions = {
    "1_vs_1": "Clade A vs. A reference",
    "1_vs_0": "Clade A vs. B reference",
    "0_vs_1": "Clade B vs. A reference",
    "0_vs_0": "Clade B vs. B reference",
}

plot_order = [
    "Clade A vs. A reference",
    "Clade A vs. B reference",
    "Clade B vs. A reference",
    "Clade B vs. B reference",
]


# Read dN/dS values gene-first

results = []
gene_level_rows = []

for folder, label in conditions.items():
    folder_path = os.path.join(base_dir, folder)

    files = (
        sorted(glob.glob(os.path.join(folder_path, "*.csv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.tsv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.txt")))
    )

    gene_means = []

    for file in files:
        try:
            df = pd.read_csv(file, sep=None, engine="python")
        except Exception:
            continue

        if "dN/dS" not in df.columns:
            continue

        dnds = pd.to_numeric(df["dN/dS"], errors="coerce")
        dnds = dnds.replace([np.inf, -np.inf], np.nan).dropna()

        if len(dnds) == 0:
            continue

        gene_mean = dnds.mean()
        gene_means.append(gene_mean)

        gene_level_rows.append({
            "Condition": label,
            "Gene_file": os.path.basename(file),
            "Gene_mean_dnds": gene_mean,
            "N_pairwise_values": len(dnds),
        })

    gene_means = np.array(gene_means, dtype=float)

    mean_dnds = np.mean(gene_means)
    se_dnds = np.std(gene_means, ddof=1) / np.sqrt(len(gene_means))

    results.append({
        "Condition": label,
        "Mean_dnds": mean_dnds,
        "SE": se_dnds,
        "N_genes": len(gene_means),
    })

res_df = pd.DataFrame(results)
res_df = res_df.set_index("Condition").loc[plot_order].reset_index()

gene_level_df = pd.DataFrame(gene_level_rows)
gene_level_df.to_csv(out_csv, index=False)

print("\n=== Gene-first average dN/dS summary ===")
print(res_df.to_string(index=False))

print(f"\nSaved gene-level summary: {out_csv}")


# =========================
# Plot
# =========================
from matplotlib.patches import Patch

# Internal = same clade reference (AA, BB)
internal_color = "#4C78A8"   # blue

# External = opposite clade reference (AB, BA)
external_color = "#E45756"   # red

bar_colors = [
    internal_color,  # Clade A vs. A reference
    external_color,  # Clade A vs. B reference
    external_color,  # Clade B vs. A reference
    internal_color,  # Clade B vs. B reference
]

plt.figure(figsize=(8, 6))

plt.bar(
    res_df["Condition"],
    res_df["Mean_dnds"],
    yerr=res_df["SE"],
    capsize=6,
    color=bar_colors,
    edgecolor="black"
)

# Legend
legend_elements = [
    Patch(facecolor=internal_color, edgecolor="black", label="Internal reference"),
    Patch(facecolor=external_color, edgecolor="black", label="External reference"),
]

plt.legend(
    handles=legend_elements,
    frameon=False,
    loc="upper right"
)

plt.title("Average dN/dS Across Conditions")
plt.xlabel("Reference type")
plt.ylabel("Average gene dN/dS")

plt.xticks(rotation=25, ha="right")
plt.grid(False)

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()

print(f"Saved plot: {out_png}")