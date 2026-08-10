#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt


input_csv = "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_4groups_clade_seperated/metadata/allGenomesAllNichesMetadata_BronchMapped.csv"
out_png = "clade_niche_distributions_v2.png"


CF = {"CF", "early.CF"}

CHRONIC = {"CF", "early.CF", "Bronchiectasis"}

NON_HUMAN = {
    "Animal", "Aquatic", "Built.Environment", "Hospital.Environment",
    "Ocean", "Plant", "Terrestrial", "Waste.Water",
}

HUMAN = {
    "Blood", "Bronchiectasis", "Burn", "Ear", "Eye", "Lung", "Pneumonia",
    "Rectal/Feces", "Skin", "Sputum", "Throat", "Upper.Respiratory",
    "Urinary", "Wound", "CF", "early.CF"
}


df = pd.read_csv(input_csv)

df["Clade"] = df["Clade"].astype(str).str.strip()
df["Niche"] = df["Niche"].astype(str).str.strip()

df = df[df["Clade"].isin(["CladeA", "CladeB"])].copy()


# Classification
df["cf_noncf"] = df["Niche"].apply(
    lambda x: "CF" if x in CF else "Non-CF"
)

def human_nonhuman(x):
    if x in HUMAN:
        return "Human"
    elif x in NON_HUMAN:
        return "Environmental"
    else:
        return "Other"

df["human_nonhuman"] = df["Niche"].apply(human_nonhuman)

def chronic_other_human(x):
    if x in CHRONIC:
        return "Chronic"
    elif x in HUMAN:
        return "Other human"
    else:
        return "Environmental/Other"

df["chronic_other_human"] = df["Niche"].apply(chronic_other_human)

# =========================
# Plot helper
# =========================
def stacked_percent_plot(ax, data, class_col, categories, title, colors, drop_other=False):
    d = data.copy()

    if drop_other:
        d = d[d[class_col].isin(categories)].copy()

    counts = (
        d.groupby(["Clade", class_col])
         .size()
         .unstack(fill_value=0)
         .reindex(index=["CladeA", "CladeB"], fill_value=0)
    )

    for cat in categories:
        if cat not in counts.columns:
            counts[cat] = 0

    counts = counts[categories]

    totals = counts.sum(axis=1)
    percents = counts.div(totals.replace(0, pd.NA), axis=0) * 100

    percents.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        edgecolor="black",
        width=0.7
    )

    ax.set_title(title, fontsize=11, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("Percentage of samples")
    ax.set_ylim(0, 112)
    ax.tick_params(axis="x", rotation=0)

    ax.legend(
        frameon=False,
        fontsize=9,
        title="",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=len(categories)
    )

    # raw counts inside the percentage bars
    for cat_idx, container in enumerate(ax.containers):
        cat = categories[cat_idx]
        labels = []

        for bar_idx, height in enumerate(container.datavalues):
            clade = counts.index[bar_idx]
            raw_count = counts.loc[clade, cat]

            if raw_count > 0 and pd.notna(height):
                labels.append(str(int(raw_count)))
            else:
                labels.append("")

        ax.bar_label(
            container,
            labels=labels,
            label_type="center",
            fontsize=8
        )

    # total n above each bar
    for i, clade in enumerate(counts.index):
        ax.text(
            i,
            104,
            f"n={int(totals.loc[clade])}",
            ha="center",
            va="bottom",
            fontsize=9
        )


# Create 1x3 figure
fig, axes = plt.subplots(1, 3, figsize=(13, 4.8))

stacked_percent_plot(
    axes[0],
    df,
    "cf_noncf",
    ["CF", "Non-CF"],
    "CF vs Non-CF",
    ["#E76F51", "#A8DADC"]
)

stacked_percent_plot(
    axes[1],
    df,
    "human_nonhuman",
    ["Human", "Environmental"],
    "Human vs Environmental",
    ["#457B9D", "#90BE6D"],
    drop_other=True
)

stacked_percent_plot(
    axes[2],
    df,
    "chronic_other_human",
    ["Chronic", "Other human"],
    "Chronic vs Other human",
    ["#E9C46A", "#2A9D8F"],
    drop_other=True
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.88, wspace=0.35)

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {out_png}")