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
PARENT_DIR = "/storage/home/hcoda1/2/bli629/scratch/vfdb/workflow_4groups_clade_seperated/clade_A/cladeA_seperated_acute_vs_chronic/dnds_output"

# 0 = chronic, 1 = acute
DIRS = {
    "00": f"{PARENT_DIR}/0_vs_0",   # chronic query vs chronic ref
    "01": f"{PARENT_DIR}/0_vs_1",   # chronic query vs acute ref
    "10": f"{PARENT_DIR}/1_vs_0",   # acute query vs chronic ref
    "11": f"{PARENT_DIR}/1_vs_1",   # acute query vs acute ref
}

BLACKLIST_GENES = {
    "VFG014147",
}

VF_ANYWHERE_RE = re.compile(r"(VFG\d+)", re.IGNORECASE)

# -----------------------
# VFG -> gene name map
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

def _paren_or_none(text: str):
    m = re.search(r"\(([^)]+)\)", str(text))
    return m.group(1).strip() if m else None

# def label_from_vfg(val: str) -> str:
#     s = str(val)
#     m = re.search(r"(VFG\d+)", s, flags=re.I)
#     if not m:
#         return s
#     vfg = m.group(1).upper()
#     gene_name = vfg_to_gene.get(vfg)
#     if gene_name:
#         inside = _paren_or_none(gene_name)
#         if inside:
#             return inside
#         return str(gene_name).split(",")[0]
#     return vfg

def label_from_vfg(val: str) -> str:
    s = str(val)
    m = re.search(r"(VFG\d+)", s, flags=re.I)
    if not m:
        return s
    return m.group(1).upper()

def full_gene_from_vfg(val: str) -> str:
    s = str(val)
    m = re.search(r"(VFG\d+)", s, flags=re.I)
    if not m:
        return s
    return vfg_to_gene.get(m.group(1).upper(), s)

# -----------------------
# Parsing functions
# -----------------------
def read_one_vf_mean_dnds(path: str) -> float:
    df = pd.read_csv(path, sep=None, engine="python")
    norm = {str(c).strip().lower(): c for c in df.columns}

    if "dn/ds" in norm:
        col = norm["dn/ds"]
    elif "dnds" in norm:
        col = norm["dnds"]
    elif "omega" in norm:
        col = norm["omega"]
    else:
        col = df.columns[-1]

    s = pd.to_numeric(df[col], errors="coerce")
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    s = s[s >= 0]

    if s.empty:
        return np.nan
    return float(s.mean())

def folder_to_gene_means(folder: str) -> pd.DataFrame:
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Folder not found: {folder}")

    files = glob.glob(os.path.join(folder, "*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files in: {folder}")

    rows = []
    for f in sorted(files):
        base = os.path.basename(f)
        m = VF_ANYWHERE_RE.search(base)
        if not m:
            continue
        gene = m.group(1).upper()

        try:
            avg = read_one_vf_mean_dnds(f)
        except Exception:
            avg = np.nan

        rows.append((gene, avg))

    if not rows:
        raise RuntimeError(f"No VF gene files found in: {folder}")

    df = pd.DataFrame(rows, columns=["gene", "avg_dnds"])

    df = (
        df.groupby("gene", as_index=False)["avg_dnds"]
          .agg(lambda x: float(np.nanmean(x.values)) if np.isfinite(x).any() else np.nan)
    )

    df = df[~df["gene"].isin(BLACKLIST_GENES)].copy()
    return df

def assert_same_gene_set(tag: str, df1: pd.DataFrame, df2: pd.DataFrame):
    s1 = set(df1["gene"])
    s2 = set(df2["gene"])
    if s1 != s2:
        only1 = sorted(list(s1 - s2))[:20]
        only2 = sorted(list(s2 - s1))[:20]
        raise RuntimeError(
            f"[{tag}] Gene sets differ!\n"
            f"  in first only:  {len(s1 - s2)} (examples: {only1})\n"
            f"  in second only: {len(s2 - s1)} (examples: {only2})\n"
        )

# -----------------------
# Pseudo-log helpers
# -----------------------
def pseudolog_forward(x, linthresh=0.01):
    x = np.asarray(x, dtype=float)
    return np.log10(1.0 + x / linthresh)

def pseudolog_inverse(y, linthresh=0.01):
    y = np.asarray(y, dtype=float)
    return linthresh * (np.power(10.0, y) - 1.0)

def make_pseudolog_ticks(lim):
    base = np.array([
        0.0,
        0.001, 0.002, 0.005,
        0.01, 0.02, 0.05,
        0.1, 0.2, 0.5,
        1.0, 2.0, 5.0
    ], dtype=float)
    ticks = base[base <= lim]
    if lim not in ticks:
        ticks = np.r_[ticks, lim]
    ticks = np.unique(np.round(ticks, 6))
    return ticks

def format_tick(v):
    if abs(v) < 1e-12:
        return "0"
    if v < 0.01:
        return f"{v:.3f}".rstrip("0").rstrip(".")
    if v < 0.1:
        return f"{v:.2f}".rstrip("0").rstrip(".")
    if v < 1:
        return f"{v:.2f}".rstrip("0").rstrip(".")
    if abs(v - round(v)) < 1e-9:
        return str(int(round(v)))
    return f"{v:.1f}".rstrip("0").rstrip(".")

# -----------------------
# General plotting function
# -----------------------
def plot_pairwise_scatter(
    df_x,
    df_y,
    out_png: str,
    title: str,
    xlabel: str,
    ylabel: str,
    label_top_n: int,
    gamma: float,
    linthresh: float
):
    assert_same_gene_set(title, df_x, df_y)

    m = df_x.merge(
        df_y,
        on="gene",
        how="outer",
        suffixes=("_x", "_y")
    )

    m["Gene"] = m["gene"].apply(full_gene_from_vfg)

    plot_df = m[
        np.isfinite(m["avg_dnds_x"]) &
        np.isfinite(m["avg_dnds_y"])
    ].copy()

    if plot_df.empty:
        print(f"No valid genes to plot for {title}")
        return

    x = plot_df["avg_dnds_x"].to_numpy()
    y = plot_df["avg_dnds_y"].to_numpy()

    above = y > x
    dist = np.abs(y - x)
    dmax = np.nanpercentile(dist, 95) if dist.size else 1.0

    # Two colormaps:
    # above diagonal -> directional selection
    # below diagonal -> relaxed selection
    norm = Normalize(0, 1)
    cmap_above = cm.YlOrBr
    cmap_below = cm.GnBu

    colors = [
        (cmap_above if above[i] else cmap_below)(
            (dist[i] / dmax) ** gamma if dmax > 0 else 1
        )
        for i in range(len(x))
    ]

    lim = np.nanmax(np.r_[x, y])
    if not np.isfinite(lim) or lim == 0:
        lim = 1.0
    lim *= 1.05

    fig, ax = plt.subplots(figsize=(10, 8))
    fig.subplots_adjust(left=0.18, right=0.84, top=0.90, bottom=0.12)

    ax.scatter(
        x, y,
        c=colors,
        s=28,
        edgecolors="white",
        linewidths=0.25,
        alpha=0.90
    )

    # -----------------------
    # Color bars
    # -----------------------
    sm_above = cm.ScalarMappable(norm=norm, cmap=cmap_above)
    sm_below = cm.ScalarMappable(norm=norm, cmap=cmap_below)
    sm_above.set_array([])
    sm_below.set_array([])

    # Need current axis position after subplots_adjust
    pos = ax.get_position()

    # Left colorbar
    cax_left = fig.add_axes([
        pos.x0 - 0.08,
        pos.y0,
        0.025,
        pos.height
    ])
    cbar_left = fig.colorbar(sm_above, cax=cax_left)
    cbar_left.set_label("Strength of Directional Selection")
    cbar_left.ax.yaxis.set_label_position("left")
    cbar_left.ax.yaxis.tick_left()

    # Right colorbar
    cax_right = fig.add_axes([
        pos.x1 + 0.03,
        pos.y0,
        0.025,
        pos.height
    ])
    cbar_right = fig.colorbar(sm_below, cax=cax_right)
    cbar_right.set_label("Strength of Relaxed Selection")

    # -----------------------
    # Axes
    # -----------------------
    ax.set_xscale("function", functions=(
        lambda v: pseudolog_forward(v, linthresh=linthresh),
        lambda v: pseudolog_inverse(v, linthresh=linthresh)
    ))
    ax.set_yscale("function", functions=(
        lambda v: pseudolog_forward(v, linthresh=linthresh),
        lambda v: pseudolog_inverse(v, linthresh=linthresh)
    ))

    ax.axhline(1, color="black", linestyle="-", linewidth=1)
    ax.axvline(1, color="black", linestyle="-", linewidth=1)
    ax.plot([0, lim], [0, lim], color="black", linestyle="--", linewidth=1)

    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")

    ticks = make_pseudolog_ticks(lim)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([format_tick(t) for t in ticks])
    ax.set_yticklabels([format_tick(t) for t in ticks])

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # -----------------------
    # Direction labels inside figure box
    # -----------------------
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

    # -----------------------
    # Labels
    # -----------------------
    if label_top_n and label_top_n > 0:
        plot_df["diag_dist"] = np.abs(plot_df["avg_dnds_y"] - plot_df["avg_dnds_x"])
        top = plot_df.sort_values("diag_dist", ascending=False).head(label_top_n)

        texts = []
        for _, row in top.iterrows():
            ann = ax.annotate(
                label_from_vfg(row["gene"]),
                (row["avg_dnds_x"], row["avg_dnds_y"]),
                xytext=(row["avg_dnds_x"] + lim * 0.01, row["avg_dnds_y"] + lim * 0.01),
                textcoords="data",
                fontsize=8,
                ha="left", va="bottom",
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white", alpha=0.9)],
                arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5)
            )
            texts.append(ann)

        if texts:
            adjust_text(
                texts,
                x=x, y=y, ax=ax,
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
    plt.close(fig)

    out_csv = os.path.splitext(out_png)[0] + ".points.csv"
    m[["gene", "Gene", "avg_dnds_x", "avg_dnds_y"]].to_csv(out_csv, index=False)

    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_csv}")
    print(f"Genes plotted: {len(plot_df)}")

# -----------------------
# Main
# -----------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="cross_reference_plots_v2")
    ap.add_argument("--label_top_n", type=int, default=20)
    ap.add_argument("--gamma", type=float, default=0.75)
    ap.add_argument("--linthresh", type=float, default=0.01)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df00 = folder_to_gene_means(DIRS["00"])
    df01 = folder_to_gene_means(DIRS["01"])
    df10 = folder_to_gene_means(DIRS["10"])
    df11 = folder_to_gene_means(DIRS["11"])

    # Plot 1: chronic query, compare chronic ref vs acute ref
    plot_pairwise_scatter(
        df_x=df00,
        df_y=df01,
        out_png=os.path.join(args.outdir, "cladeA_chronic_00_vs_01.scatter.png"),
        title="Clade A chronic: chronic ref (00) vs acute ref (01)",
        xlabel="Clade A chronic vs chronic reference (00)",
        ylabel="Clade A chronic vs acute reference (01)",
        label_top_n=args.label_top_n,
        gamma=args.gamma,
        linthresh=args.linthresh
    )

    # Plot 2: acute query, compare acute ref vs chronic ref
    plot_pairwise_scatter(
        df_x=df11,
        df_y=df10,
        out_png=os.path.join(args.outdir, "cladeA_acute_11_vs_10.scatter.png"),
        title="Clade A acute: acute ref (11) vs chronic ref (10)",
        xlabel="Clade A acute vs acute reference (11)",
        ylabel="Clade A acute vs chronic reference (10)",
        label_top_n=args.label_top_n,
        gamma=args.gamma,
        linthresh=args.linthresh
    )

if __name__ == "__main__":
    main()