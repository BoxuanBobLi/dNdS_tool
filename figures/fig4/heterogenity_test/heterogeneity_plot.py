#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# INPUT
# ==========================
CSV = "/data2/B_Li/vfdb/vgrG/hypothesis_test/vgrG_acute_vs_chronic_heterogeneity_summary.csv"   # replace with your csv
OUTPNG = "vgrG_clean_compare_v2.png"
OUTPDF = "vgrG_clean_compare_v2.pdf"

df = pd.read_csv(CSV)

keep = [
    "VFG041012.nogap.dedup.nt_ali.fasta",
    "VFG051487.nogap.dedup.nt_ali.fasta",
    "VFG051489.nogap.dedup.nt_ali.fasta",
    "VFG051491.nogap.dedup.nt_ali.fasta",
    "VFG051593.nogap.dedup.nt_ali.fasta",
    "VFG051596.nogap.dedup.nt_ali.fasta",
]
df = df[df["gene"].isin(keep)].copy()

df["gene_short"] = df["gene"].str.replace(".nogap.dedup.nt_ali.fasta", "", regex=False)

# order by chronic - acute diversity
df["diff_pw"] = df["chronic_mean_pairwise_dist"] - df["acute_mean_pairwise_dist"]
df = df.sort_values("diff_pw")

genes = df["gene_short"].tolist()
y = np.arange(len(df))

# consistent state colors
acute_color = "#1f77b4"
chronic_color = "#d62728"

fig, axes = plt.subplots(
    3, 1, figsize=(8, 10),
    sharey=False,
    constrained_layout=True
)

fig.suptitle(
    "vgrG heterogeneity test",
    fontsize=18,
    y=1.05
)

fig.set_constrained_layout_pads(h_pad=0.08, w_pad=0.04, hspace=0.08)

# --------------------------
# common style
# --------------------------
for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.2, linewidth=0.8)

# ==========================
# A. within-group diversity
# ==========================
ax = axes[0]

for i, row in enumerate(df.itertuples(index=False)):
    ax.plot(
        [row.acute_mean_pairwise_dist, row.chronic_mean_pairwise_dist],
        [i, i],
        color="0.75",
        linewidth=1.2,
        zorder=1
    )

ax.scatter(
    df["acute_mean_pairwise_dist"], y,
    s=55, facecolor=acute_color, edgecolor=acute_color,
    label="Acute", zorder=3
)
ax.scatter(
    df["chronic_mean_pairwise_dist"], y,
    s=55, facecolor=chronic_color, edgecolor=chronic_color,
    label="Chronic", zorder=3
)

ax.set_yticks(y)
ax.set_yticklabels(genes)
ax.set_xlabel("Mean pairwise distance")
ax.set_title("A. Within-group sequence diversity")
ax.legend(frameon=False, loc="best")

# ==========================
# B. acute sequences
# ==========================
ax = axes[1]

for i, row in enumerate(df.itertuples(index=False)):
    ax.plot(
        [row.acute_to_acute_consensus, row.acute_to_chronic_consensus],
        [i, i],
        color="0.75",
        linewidth=1.2,
        zorder=1
    )

# filled = own ref
ax.scatter(
    df["acute_to_acute_consensus"], y,
    s=55, facecolor=acute_color, edgecolor=acute_color,
    label="Own ref", zorder=3
)

# open = other ref
ax.scatter(
    df["acute_to_chronic_consensus"], y,
    s=55, facecolor="white", edgecolor=acute_color, linewidth=1.5,
    label="Other ref", zorder=3
)

ax.set_yticks(y)
ax.set_yticklabels(genes)
ax.set_xlabel("Distance to reference")
ax.set_title("B. Acute sequences")
ax.legend(frameon=False, loc="best")

# ==========================
# C. chronic sequences
# ==========================
ax = axes[2]

for i, row in enumerate(df.itertuples(index=False)):
    ax.plot(
        [row.chronic_to_chronic_consensus, row.chronic_to_acute_consensus],
        [i, i],
        color="0.75",
        linewidth=1.2,
        zorder=1
    )

ax.scatter(
    df["chronic_to_chronic_consensus"], y,
    s=55, facecolor=chronic_color, edgecolor=chronic_color,
    label="Own ref", zorder=3
)

ax.scatter(
    df["chronic_to_acute_consensus"], y,
    s=55, facecolor="white", edgecolor=chronic_color, linewidth=1.5,
    label="Other ref", zorder=3
)

ax.set_yticks(y)
ax.set_yticklabels(genes)
ax.set_xlabel("Distance to reference")
ax.set_title("C. Chronic sequences")
ax.legend(frameon=False, loc="best")

plt.savefig(OUTPNG, dpi=300, bbox_inches="tight")
plt.savefig(OUTPDF, bbox_inches="tight")
print(f"Saved: {OUTPNG}")
print(f"Saved: {OUTPDF}")