#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.patches import Patch

# =========================
# Input
# =========================
base_dir = "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output"
top_n = 50

out_png_combined = "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/visualization/venn/venn_top50_combined_cladeA_cladeB.png"

conditions = {
    "1_vs_1": "AA",
    "1_vs_0": "AB",
    "0_vs_1": "BA",
    "0_vs_0": "BB",
}

# =========================
# Colors
# =========================
REF_A_COLOR = "#9EC3F5"      # blue for ref A
REF_B_COLOR = "#F4A261"      # orange for ref B
OVERLAP_COLOR = "#E8D9B5"    # pale yellow overlap

# =========================
# Helpers
# =========================
def get_gene_name(file):
    name = os.path.basename(file)
    return name.split(".")[0]

def read_gene_means(folder):
    folder_path = os.path.join(base_dir, folder)

    files = (
        sorted(glob.glob(os.path.join(folder_path, "*.csv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.tsv"))) +
        sorted(glob.glob(os.path.join(folder_path, "*.txt")))
    )

    rows = []

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

        rows.append({
            "Gene": get_gene_name(file),
            "mean_dnds": dnds.mean(),
            "n_values": len(dnds),
        })

    return pd.DataFrame(rows)

def top_genes(gene_df, n=50):
    return set(
        gene_df.sort_values("mean_dnds", ascending=False)
               .head(n)["Gene"]
    )

def draw_venn_on_ax(ax, set_left, set_right, left_label, right_label, title):
    # Left circle = ref A = blue
    # Right circle = ref B = orange
    v = venn2(
        [set_left, set_right],
        set_labels=(left_label, right_label),
        set_colors=(REF_A_COLOR, REF_B_COLOR),
        alpha=1.0,
        ax=ax
    )

    # Left-only region: ref A
    patch = v.get_patch_by_id("10")
    if patch is not None:
        patch.set_color(REF_A_COLOR)
        patch.set_alpha(0.75)

    # Right-only region: ref B
    patch = v.get_patch_by_id("01")
    if patch is not None:
        patch.set_color(REF_B_COLOR)
        patch.set_alpha(0.75)

    # Overlap region: yellow
    patch = v.get_patch_by_id("11")
    if patch is not None:
        patch.set_color(OVERLAP_COLOR)
        patch.set_alpha(0.85)

    # Style number labels
    for subset_id in ("10", "01", "11"):
        label = v.get_label_by_id(subset_id)
        if label is not None:
            label.set_fontsize(18)
            label.set_color("#222222")

    # Style set labels
    for label in v.set_labels:
        if label is not None:
            label.set_fontsize(18)
            label.set_color("#222222")

    # Move bottom labels slightly toward center
    if v.set_labels[0] is not None:
        v.set_labels[0].set_position((-0.42, -0.58))
        v.set_labels[0].set_ha("center")

    if v.set_labels[1] is not None:
        v.set_labels[1].set_position((0.42, -0.58))
        v.set_labels[1].set_ha("center")

    ax.set_title(title, fontsize=20, pad=22, color="#222222")

    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.75, 0.65)
    ax.set_aspect("equal")
    ax.axis("off")

    left_only = len(set_left - set_right)
    right_only = len(set_right - set_left)
    overlap = len(set_left & set_right)

    print(f"{title}: left only={left_only}, overlap={overlap}, right only={right_only}")

# =========================
# Read gene-level means
# =========================
gene_means = {}

for folder, label in conditions.items():
    gene_means[label] = read_gene_means(folder)

# =========================
# Top genes
# =========================
top_AA = top_genes(gene_means["AA"], top_n)
top_AB = top_genes(gene_means["AB"], top_n)
top_BA = top_genes(gene_means["BA"], top_n)
top_BB = top_genes(gene_means["BB"], top_n)

# =========================
# Combined figure
# =========================
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Clade A query: ref A vs ref B
draw_venn_on_ax(
    ax=axes[0],
    set_left=top_AA,
    set_right=top_AB,
    left_label="Top 50 genes\nwith ref A",
    right_label="Top 50 genes\nwith ref B",
    title="Clade A"
)

# Clade B query: ref A vs ref B
draw_venn_on_ax(
    ax=axes[1],
    set_left=top_BA,
    set_right=top_BB,
    left_label="Top 50 genes\nwith ref A",
    right_label="Top 50 genes\nwith ref B",
    title="Clade B"
)

# =========================
# Legend
# =========================
legend_handles = [
    Patch(facecolor=REF_A_COLOR, edgecolor="none", alpha=0.75, label="Reference A"),
    Patch(facecolor=REF_B_COLOR, edgecolor="none", alpha=0.75, label="Reference B"),
    Patch(facecolor=OVERLAP_COLOR, edgecolor="none", alpha=0.85, label="Overlap"),
]

fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=3,
    frameon=False,
    fontsize=16,
    bbox_to_anchor=(0.5, 0.02)
)

plt.subplots_adjust(wspace=0.25, bottom=0.16, top=0.88)

plt.savefig(out_png_combined, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()

print(f"Saved: {out_png_combined}")