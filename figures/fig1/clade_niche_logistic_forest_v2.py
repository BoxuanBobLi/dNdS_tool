#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

# =========================
# Load data
# =========================
df = pd.read_csv(
    "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_4groups_clade_seperated/metadata/allGenomesAllNichesMetadata_BronchMapped.csv"
)

df["Clade"] = df["Clade"].astype(str).str.strip()
df["Niche"] = df["Niche"].astype(str).str.strip()

df = df[df["Clade"].isin(["CladeA", "CladeB"])].copy()

# OR > 1 means enriched in CladeA relative to CladeB
df["Clade_bin"] = (df["Clade"] == "CladeA").astype(int)

# =========================
# Define niche sets
# =========================
CF = {"CF", "early.CF"}

CHRONIC = {"CF", "early.CF", "Bronchiectasis"}

NON_HUMAN = {
    "Animal",
    "Aquatic",
    "Built.Environment",
    "Hospital.Environment",
    "Ocean",
    "Plant",
    "Terrestrial",
    "Waste.Water",
}

HUMAN = {
    "Blood",
    "Bronchiectasis",
    "Burn",
    "Ear",
    "Eye",
    "Lung",
    "Pneumonia",
    "Rectal/Feces",
    "Skin",
    "Sputum",
    "Throat",
    "Upper.Respiratory",
    "Urinary",
    "Wound",
    "CF",
    "early.CF",
}

# =========================
# Build outcomes
# =========================
df["is_CF"] = df["Niche"].isin(CF).astype(int)
df["is_human"] = df["Niche"].isin(HUMAN).astype(int)

# Chronic vs other human:
# keep only human isolates, then compare chronic vs non-chronic human
df_chronic_human = df[df["Niche"].isin(HUMAN)].copy()
df_chronic_human["is_chronic_vs_other_human"] = (
    df_chronic_human["Niche"].isin(CHRONIC)
).astype(int)

# =========================
# Logistic regression helper
# =========================
def run_logistic(y: pd.Series, X: pd.DataFrame):
    X = sm.add_constant(X, has_constant="add")
    model = sm.Logit(y, X).fit(disp=0)

    coef = model.params["Clade_bin"]
    se = model.bse["Clade_bin"]
    pval = model.pvalues["Clade_bin"]

    OR = float(np.exp(coef))
    CI_low = float(np.exp(coef - 1.96 * se))
    CI_high = float(np.exp(coef + 1.96 * se))

    return OR, CI_low, CI_high, float(pval)

# =========================
# Count table helper
# =========================
def print_2x2_table(
    title: str,
    frame: pd.DataFrame,
    outcome_col: str,
    positive_label: str,
    negative_label: str
):
    tab = pd.crosstab(frame["Clade"], frame[outcome_col])
    tab = tab.reindex(index=["CladeA", "CladeB"], fill_value=0)
    tab = tab.reindex(columns=[0, 1], fill_value=0)
    tab.columns = [negative_label, positive_label]

    print(f"\n=== {title} ===")
    print(tab)

# =========================
# Run models
# =========================
results = []

comparisons = [
    ("CF vs Non-CF", df, "is_CF"),
    ("Human vs Environmental", df, "is_human"),
    ("Chronic vs other human", df_chronic_human, "is_chronic_vs_other_human"),
]

for name, frame, ycol in comparisons:
    OR, low, high, pval = run_logistic(frame[ycol], frame[["Clade_bin"]])
    results.append((name, OR, low, high, pval, len(frame)))

res_df = pd.DataFrame(
    results,
    columns=["Comparison", "OR", "CI_low", "CI_high", "P_value", "N_used"]
)

print("\n=== Logistic regression results ===")
print(res_df.to_string(index=False))

print_2x2_table("CF vs Non-CF", df, "is_CF", "CF", "Non-CF")
print_2x2_table("Human vs Environmental", df, "is_human", "Human", "Environmental")
print_2x2_table(
    "Chronic vs other human",
    df_chronic_human,
    "is_chronic_vs_other_human",
    "Chronic",
    "Other human"
)

# =========================
# Plot forest plot
# =========================
plot_order = [
    "Human vs Environmental",
    "Chronic vs other human",
    "CF vs Non-CF",
]

plot_df = res_df.set_index("Comparison").loc[plot_order].reset_index()

fig, ax = plt.subplots(figsize=(7, 4.2))

y_pos = np.arange(len(plot_df))

ax.errorbar(
    plot_df["OR"],
    y_pos,
    xerr=[
        plot_df["OR"] - plot_df["CI_low"],
        plot_df["CI_high"] - plot_df["OR"]
    ],
    fmt="o",
    color="black",
    ecolor="black",
    capsize=4,
    markersize=8,
    linewidth=2
)

ax.axvline(1, linestyle="--", color="gray", linewidth=1.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df["Comparison"])

ax.set_xlabel("Odds ratio", fontsize=11, labelpad=14)
ax.set_title("Clade predicts niche association")

# Put Clade labels on the same baseline as the x-axis label
ax.text(
    0.0, -0.16, "Clade B",
    transform=ax.transAxes,
    ha="left",
    va="center",
    fontsize=10
)

ax.text(
    1.0, -0.16, "Clade A",
    transform=ax.transAxes,
    ha="right",
    va="center",
    fontsize=10
)

plt.tight_layout()
plt.savefig("clade_niche_logistic_forest_v2.png", dpi=300, bbox_inches="tight")
plt.close()

print("\nSaved: clade_niche_logistic_forest_v2.png")