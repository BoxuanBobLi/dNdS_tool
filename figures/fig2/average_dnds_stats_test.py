#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu, ttest_rel, wilcoxon

# =========================
# Input
# =========================
base_dir = "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output"
out_csv = "gene_first_dnds_stats.csv"

conditions = {
    "1_vs_1": "AA_avg_dnds",
    "1_vs_0": "AB_avg_dnds",
    "0_vs_1": "BA_avg_dnds",
    "0_vs_0": "BB_avg_dnds",
}

# =========================
# Helper: extract gene name from file
# =========================
def get_gene_name(file):
    name = os.path.basename(file)
    return name.split(".")[0]

# =========================
# Read gene-level mean dN/dS
# =========================
gene_rows = []

for folder, label in conditions.items():
    folder_path = os.path.join(base_dir, folder)

    files = (
        sorted(glob.glob(os.path.join(folder_path, "*.csv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.tsv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.txt")))
    )

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

        gene_rows.append({
            "Gene": get_gene_name(file),
            "Condition": label,
            "mean_dnds": dnds.mean(),
            "n_values": len(dnds),
        })

gene_df = pd.DataFrame(gene_rows)

# =========================
# Convert to wide format
# =========================
wide = gene_df.pivot_table(
    index="Gene",
    columns="Condition",
    values="mean_dnds",
    aggfunc="mean"
).reset_index()

needed_cols = [
    "AA_avg_dnds",
    "AB_avg_dnds",
    "BA_avg_dnds",
    "BB_avg_dnds",
]

wide = wide.dropna(subset=needed_cols).copy()

wide["inter_clade_avg"] = (
    wide["AB_avg_dnds"] + wide["BA_avg_dnds"]
) / 2

wide["intra_clade_avg"] = (
    wide["AA_avg_dnds"] + wide["BB_avg_dnds"]
) / 2

wide.to_csv(out_csv, index=False)

print(f"\nSaved gene-level table: {out_csv}")
print(f"Genes used in paired tests: {len(wide)}")

# =========================
# Test helper
# =========================
def paired_test(name, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    t_stat, t_p = ttest_rel(x, y, nan_policy="omit")
    w_stat, w_p = wilcoxon(x, y)

    print(f"\n=== {name} ===")
    print(f"Mean x = {np.mean(x):.6f}")
    print(f"Mean y = {np.mean(y):.6f}")
    print(f"Mean difference x-y = {np.mean(x - y):.6f}")
    print(f"Paired t-test p = {t_p:.3e}")
    print(f"Wilcoxon signed-rank p = {w_p:.3e}")

def unpaired_test(name, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    t_stat, t_p = ttest_ind(x, y, equal_var=False, nan_policy="omit")
    u_stat, u_p = mannwhitneyu(x, y, alternative="two-sided")

    print(f"\n=== {name} ===")
    print(f"Mean x = {np.mean(x):.6f}")
    print(f"Mean y = {np.mean(y):.6f}")
    print(f"Mean difference x-y = {np.mean(x) - np.mean(y):.6f}")
    print(f"Welch t-test p = {t_p:.3e}")
    print(f"Mann-Whitney U p = {u_p:.3e}")

# =========================
# Paired tests, recommended
# =========================
paired_test(
    "AB vs AA",
    wide["AB_avg_dnds"],
    wide["AA_avg_dnds"]
)

paired_test(
    "BA vs BB",
    wide["BA_avg_dnds"],
    wide["BB_avg_dnds"]
)

paired_test(
    "Inter-clade average vs intra-clade average",
    wide["inter_clade_avg"],
    wide["intra_clade_avg"]
)

# =========================
# Optional unpaired pooled test
# =========================
inter_pooled = pd.concat([
    wide["AB_avg_dnds"],
    wide["BA_avg_dnds"]
])

intra_pooled = pd.concat([
    wide["AA_avg_dnds"],
    wide["BB_avg_dnds"]
])

unpaired_test(
    "Pooled inter-clade vs pooled intra-clade",
    inter_pooled,
    intra_pooled
)