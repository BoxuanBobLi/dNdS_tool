#!/usr/bin/env python3

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT_CSV = "VF_presence_absence_enrichment_Clade_A_vs_B_all.csv"
OUTPUT_PNG = "VF_enrichment_effect_size_FDR_volcano.png"

FDR_THRESHOLD = 0.05
EFFECT_THRESHOLD = 5.0       # percentage points
Y_CAP = 300                  # cap extremely small FDR values for plotting
N_LABELS = 6                 # labels on each side


# ============================================================
# HELPERS
# ============================================================

def extract_vfg(accession):
    """Extract VFGxxxxxx from the Accession column."""
    if pd.isna(accession):
        return ""
    match = re.search(r"(VFG\d+)", str(accession))
    return match.group(1) if match else str(accession)


def make_label(row):
    """
    Use PA locus tag when available.
    Otherwise use the VFG ID.
    """
    if pd.notna(row["locus_tag"]) and str(row["locus_tag"]).strip() != "":
        return str(row["locus_tag"])

    return extract_vfg(row["Accession"])


# ============================================================
# READ DATA
# ============================================================

df = pd.read_csv(INPUT_CSV)

# x-axis: prevalence difference in percentage points
df["x"] = (
    pd.to_numeric(
        df["freq_difference_CladeA_minus_CladeB"],
        errors="coerce"
    ) * 100
)

# y-axis: -log10(BH-adjusted FDR)
fdr = pd.to_numeric(df["fdr_bh"], errors="coerce")

# Some FDR values are exactly 0 because they are extremely small.
# Replace 0 with the smallest positive floating-point value.
fdr_for_log = fdr.copy()
fdr_for_log.loc[fdr_for_log <= 0] = np.nextafter(0, 1)

df["minus_log10_fdr"] = -np.log10(fdr_for_log)

# Cap the displayed y-value
df["y_plot"] = df["minus_log10_fdr"].clip(upper=Y_CAP)


# ============================================================
# CLASSIFY POINTS
# ============================================================

df["category"] = "Not significant"

df.loc[
    (df["significant_fdr"] == True)
    & (df["significant_and_meaningful"] == False),
    "category"
] = "Significant, <5% difference"

df.loc[
    df["enriched_clade"] == "Clade_A",
    "category"
] = "Meaningful: Clade A"

df.loc[
    df["enriched_clade"] == "Clade_B",
    "category"
] = "Meaningful: Clade B"


# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))

categories = [
    ("Not significant", "tab:blue"),
    ("Significant, <5% difference", "tab:orange"),
    ("Meaningful: Clade A", "tab:green"),
    ("Meaningful: Clade B", "tab:red"),
]

for category, color in categories:

    sub = df[df["category"] == category]

    ax.scatter(
        sub["x"],
        sub["y_plot"],
        s=25,
        alpha=0.65,
        color=color,
        label=f"{category} (n={len(sub)})"
    )


# ============================================================
# THRESHOLD LINES
# ============================================================

# ±5 percentage-point effect-size thresholds
ax.axvline(
    -EFFECT_THRESHOLD,
    linestyle=":",
    linewidth=1.2,
    color="tab:blue"
)

ax.axvline(
    EFFECT_THRESHOLD,
    linestyle=":",
    linewidth=1.2,
    color="tab:blue"
)

# FDR = 0.05 threshold
ax.axhline(
    -np.log10(FDR_THRESHOLD),
    linestyle="--",
    linewidth=1.2,
    color="tab:blue"
)


# ============================================================
# SELECT EXTREME GENES TO LABEL
# ============================================================

# Most Clade-B-enriched = most negative prevalence differences
clade_b_labels = (
    df[df["enriched_clade"] == "Clade_B"]
    .sort_values("x", ascending=True)
    .head(N_LABELS)
    .copy()
)

# Most Clade-A-enriched = most positive prevalence differences
clade_a_labels = (
    df[df["enriched_clade"] == "Clade_A"]
    .sort_values("x", ascending=False)
    .head(N_LABELS)
    .copy()
)


# ============================================================
# LEFT-SIDE LABELS
# ============================================================

left_y_positions = np.linspace(
    285,
    125,
    len(clade_b_labels)
)

for (_, row), label_y in zip(
    clade_b_labels.iterrows(),
    left_y_positions
):

    ax.annotate(
        make_label(row),
        xy=(row["x"], row["y_plot"]),
        xytext=(-130, label_y),
        textcoords="data",
        ha="right",
        va="center",
        fontsize=10,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            linewidth=1
        ),
        annotation_clip=False
    )


# ============================================================
# RIGHT-SIDE LABELS
# ============================================================

right_y_positions = np.linspace(
    285,
    125,
    len(clade_a_labels)
)

for (_, row), label_y in zip(
    clade_a_labels.iterrows(),
    right_y_positions
):

    ax.annotate(
        make_label(row),
        xy=(row["x"], row["y_plot"]),
        xytext=(116, label_y),
        textcoords="data",
        ha="left",
        va="center",
        fontsize=10,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            linewidth=1
        ),
        annotation_clip=False
    )


# ============================================================
# FORMATTING
# ============================================================

ax.set_xlim(-102, 105)
ax.set_ylim(0, 315)

ax.set_xlabel(
    "Clade A − Clade B prevalence difference (percentage points)",
    fontsize=11
)

ax.set_ylabel(
    "−log10(BH-adjusted FDR)",
    fontsize=11
)

ax.set_title(
    "VF presence/absence enrichment by effect size and FDR",
    fontsize=13
)

ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.13),
    ncol=2,
    frameon=False
)

# Leave extra space for labels outside the axes
plt.subplots_adjust(
    left=0.18,
    right=0.82,
    bottom=0.20,
    top=0.92
)

plt.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Saved: {OUTPUT_PNG}")