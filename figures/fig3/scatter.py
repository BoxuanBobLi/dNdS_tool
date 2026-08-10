#!/usr/bin/env python3

import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from adjustText import adjust_text
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# -----------------------
# Paths
# -----------------------
DIRS = {
    "01": "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/0_vs_1",
    "00": "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/0_vs_0",
    "11": "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_1",
    "10": "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_0",
}

BLACKLIST_GENES = {"VFG014147"}
VF_ANYWHERE_RE = re.compile(r"(VFG\d+)", re.IGNORECASE)

# -----------------------
# Annotation
# -----------------------
annot_csv = "/storage/home/hcoda1/2/bli629/scratch/vfdb/db/vfdb_names.csv"
annot = pd.read_csv(annot_csv)
annot["VFG_ID"] = annot["ID"].astype(str).str.extract(r"(VFG\d+)", flags=re.I)

vfg_to_gene = (
    annot.dropna(subset=["VFG_ID"])
         .assign(VFG_ID=lambda x: x["VFG_ID"].str.upper())
         .set_index("VFG_ID")["Gene"]
         .to_dict()
)

def _paren_or_none(text):
    m = re.search(r"\(([^)]+)\)", str(text))
    return m.group(1).strip() if m else None

def label_from_vfg(val):
    m = re.search(r"(VFG\d+)", str(val), flags=re.I)
    if not m:
        return val
    vfg = m.group(1).upper()
    g = vfg_to_gene.get(vfg)
    if g:
        return _paren_or_none(g) or g.split(",")[0]
    return vfg

def full_gene_from_vfg(val):
    m = re.search(r"(VFG\d+)", str(val), flags=re.I)
    return vfg_to_gene.get(m.group(1).upper(), val) if m else val

# -----------------------
# Load data
# -----------------------
def read_one_vf_mean_dnds(path):
    df = pd.read_csv(path, sep=None, engine="python")
    norm = {str(c).strip().lower(): c for c in df.columns}

    col = norm.get("dn/ds") or norm.get("dnds") or norm.get("omega") or df.columns[-1]

    s = pd.to_numeric(df[col], errors="coerce")
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    s = s[s >= 0]

    return float(s.mean()) if not s.empty else np.nan

def folder_to_gene_means(folder):
    files = glob.glob(os.path.join(folder, "*.csv"))
    rows = []

    for f in files:
        if not f.endswith("_groupwise_dnds.csv"):
            continue
        m = VF_ANYWHERE_RE.search(f)
        if not m:
            continue
        gene = m.group(1).upper()
        rows.append((gene, read_one_vf_mean_dnds(f)))

    df = pd.DataFrame(rows, columns=["gene", "avg_dnds"])
    df = df.groupby("gene", as_index=False)["avg_dnds"].mean()
    return df[~df["gene"].isin(BLACKLIST_GENES)]

# -----------------------
# Pseudo-log
# -----------------------
def pseudolog_forward(x, linthresh=0.01):
    return np.log10(1 + x / linthresh)

def pseudolog_inverse(y, linthresh=0.01):
    return linthresh * (10**y - 1)

# -----------------------
# Plot
# -----------------------
def plot_pair(df_x, df_y, out_png, title, xlabel, ylabel, gamma, linthresh):

    m = df_x.merge(df_y, on="gene", suffixes=("_x", "_y"))
    m["Gene"] = m["gene"].apply(full_gene_from_vfg)

    x = m["avg_dnds_x"].values
    y = m["avg_dnds_y"].values

    above = y > x
    dist = np.abs(y - x)
    dmax = np.percentile(dist, 95)

    norm = Normalize(0, 1)
    cmap_above = cm.YlOrBr
    cmap_below = cm.GnBu

    colors = [
        (cmap_above if above[i] else cmap_below)(
            (dist[i] / dmax) ** gamma if dmax > 0 else 1
        )
        for i in range(len(x))
    ]

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(x, y, c=colors, s=28, edgecolors="white", linewidths=0.25)

    # ----------- COLORBARS -----------
    sm_above = cm.ScalarMappable(norm=norm, cmap=cmap_above)
    sm_below = cm.ScalarMappable(norm=norm, cmap=cmap_below)

    # RIGHT (auto)
    cbar_right = fig.colorbar(sm_below, ax=ax, fraction=0.046, pad=0.04)
    cbar_right.set_label("Strength of Relaxed Selection")

    # LEFT (manual, SAME HEIGHT)
    pos = ax.get_position()
    cax_left = fig.add_axes([
        pos.x0 - 0.08,   # move left/right here
        pos.y0,
        0.02,
        pos.height
    ])

    cbar_left = fig.colorbar(sm_above, cax=cax_left)
    cbar_left.set_label("Strength of Directional Selection")
    cbar_left.ax.yaxis.set_label_position("left")
    cbar_left.ax.yaxis.tick_left()

    # ----------- AXES -----------
    ax.set_xscale("function", functions=(
        lambda v: pseudolog_forward(v, linthresh),
        lambda v: pseudolog_inverse(v, linthresh)
    ))
    ax.set_yscale("function", functions=(
        lambda v: pseudolog_forward(v, linthresh),
        lambda v: pseudolog_inverse(v, linthresh)
    ))

    lim = max(np.nanmax(x), np.nanmax(y)) * 1.05
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)

    ax.axhline(1, color="black")
    ax.axvline(1, color="black")
    ax.plot([0, lim], [0, lim], "--", color="black")

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ----------- DIRECTION TEXT -----------
    ax.text(
        0.03, 0.96,
        "Directional",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
        color="#222222",
        path_effects=[pe.withStroke(linewidth=3, foreground="white", alpha=0.9)]
    )

    ax.text(
        0.97, 0.04,
        "Relaxed",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color="#222222",
        path_effects=[pe.withStroke(linewidth=3, foreground="white", alpha=0.9)]
    )
    # ----------- LABELS -----------
    diff_above = y - x
    diff_below = x - y

    idx_above = np.where(above)[0]
    idx_below = np.where(~above)[0]

    label_top_n_above = 15
    label_top_n_below = 10

    top_above = (
        idx_above[np.argsort(diff_above[idx_above])[::-1][:label_top_n_above]]
        if idx_above.size else np.array([], dtype=int)
    )
    top_below = (
        idx_below[np.argsort(diff_below[idx_below])[::-1][:label_top_n_below]]
        if idx_below.size else np.array([], dtype=int)
    )

    over1 = np.where((x > 1) | (y > 1))[0]
    label_idx = np.unique(np.concatenate([top_above, top_below, over1]))

    texts = []
    for li in label_idx:
        ann = ax.annotate(
            label_from_vfg(m.iloc[li]["gene"]),
            (x[li], y[li]),
            xytext=(x[li] + lim * 0.01, y[li] + lim * 0.01),
            textcoords="data",
            fontsize=8,
            ha="left",
            va="bottom",
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white", alpha=0.9)],
            arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5)
        )
        texts.append(ann)

    if texts:
        adjust_text(
            texts,
            x=x,
            y=y,
            ax=ax,
            expand_points=(1.2, 1.2),
            expand_text=(1.05, 1.2),
            force_points=0.3,
            force_text=0.5,
            only_move={"points": "y", "text": "xy"},
            autoalign=True,
            lim=300,
            arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5)
        )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------
# Main
# -----------------------
def main():
    df00 = folder_to_gene_means(DIRS["00"])
    df01 = folder_to_gene_means(DIRS["01"])
    df10 = folder_to_gene_means(DIRS["10"])
    df11 = folder_to_gene_means(DIRS["11"])

    plot_pair(
        df00, df01,
        "viz_00_vs_01.png",
        "Average dN/dS of Clade B Samples against Clade A/B References",
        "Clade B Sample against Clade B Reference",
        "Clade B Sample against Clade A Reference",
        gamma=0.75, linthresh=0.01
    )

    plot_pair(
        df11, df10,
        "viz_11_vs_10.png",
        "Average dN/dS of Clade A Samples against Clade B/A References",
        "Clade A Sample against Clade A Reference",
        "Clade A Sample against Clade B Reference",
        gamma=0.75, linthresh=0.01
    )

if __name__ == "__main__":
    main()