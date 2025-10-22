import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Normalize
import matplotlib.patheffects as pe
import numpy as np
import re
from adjustText import adjust_text

# -----------------------
# Load dN/dS summary
# -----------------------
csv_path = "/dir/to/multi_folder_dnds_summary_all_inf.csv"
df = pd.read_csv(csv_path)

# Ensure numeric columns
for col in ["00_avg_dnds", "01_avg_dnds", "10_avg_dnds", "11_avg_dnds"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# Filter all four metrics into [0, 2.5]
df = df[(df["00_avg_dnds"].between(0, 2.5)) &
        (df["01_avg_dnds"].between(0, 2.5)) &
        (df["10_avg_dnds"].between(0, 2.5)) &
        (df["11_avg_dnds"].between(0, 2.5))].copy()

# -----------------------
# Load VFG -> Gene-name map
# -----------------------
annot_csv = "/data1/B_Li/vfdb/db/vfdb_names.csv"
annot = pd.read_csv(annot_csv)
annot["VFG_ID"] = annot["ID"].astype(str).str.extract(r"(VFG\d+)", flags=re.I)

vfg_to_gene = (
    annot.dropna(subset=["VFG_ID"])
         .assign(VFG_ID=lambda x: x["VFG_ID"].str.upper())
         .set_index("VFG_ID")["Gene"]
         .to_dict()
)

def _paren_or_none(text: str) -> str | None:
    m = re.search(r"\(([^)]+)\)", str(text))
    return m.group(1).strip() if m else None

def _label_from_vfg(val: str) -> str:
    """Short label for the PLOT: content inside parentheses if available."""
    s = str(val)
    m = re.search(r"(VFG\d+)", s, flags=re.I)
    if not m:
        return s
    vfg = m.group(1).upper()
    gene_name = vfg_to_gene.get(vfg)
    if gene_name:
        inside = _paren_or_none(gene_name)
        if inside:
            return inside
        return str(gene_name).split(",")[0]
    return vfg

def _full_gene_from_vfg(val: str) -> str:
    """Full gene name for the CSV: use the entire 'Gene' field from annotation."""
    s = str(val)
    m = re.search(r"(VFG\d+)", s, flags=re.I)
    if not m:
        return s
    return vfg_to_gene.get(m.group(1).upper(), s)

# -----------------------
# Plot helper (writes labeled points to CSV)
# -----------------------
def plot_pair_like_AA_AB(
    df_all, xcol, ycol, out_png,
    gene_col="Gene",
    label_top_n_above=10, label_top_n_below=10,
    label_diff_min_above=None, label_diff_min_below=None,
    title=None,
    label_fontsize=7,
    fig_size=(10, 7),
):
    # Prepare data
    if xcol not in df_all.columns or ycol not in df_all.columns or gene_col not in df_all.columns:
        print(f"Skip: missing one of columns [{xcol}, {ycol}, {gene_col}]")
        return

    data = df_all[[gene_col, xcol, ycol]].copy()
    mask = np.isfinite(data[xcol]) & np.isfinite(data[ycol]) & (data[xcol] >= 0) & (data[ycol] >= 0)
    data = data[mask].reset_index(drop=True)
    if data.empty:
        print(f"Skip: no valid rows for {xcol} vs {ycol}")
        return

    x = data[xcol].to_numpy()
    y = data[ycol].to_numpy()

    # Above/below diagonal
    above = y > x
    below = ~above

    # Distance for color intensity
    dist = np.abs(y - x)
    dist95 = np.nanpercentile(dist, 95) if dist.size else 1.0
    distnorm = np.clip(dist / (dist95 if dist95 > 0 else 1.0), 0, 1)

    cmap_above = matplotlib.colormaps['Reds']
    cmap_below = matplotlib.colormaps['Blues']
    norm = Normalize(vmin=0, vmax=1)

    fig, ax = plt.subplots(figsize=fig_size)

    # Scatter points
    ax.scatter(x[above], y[above], c=cmap_above(distnorm[above]), s=20, edgecolors="none", alpha=0.8, label="Above diagonal")
    ax.scatter(x[below], y[below], c=cmap_below(distnorm[below]), s=20, edgecolors="none", alpha=0.8, label="Below diagonal")

    # Reference lines & diagonal
    ax.axhline(1, color='black', linestyle='-', linewidth=1)
    ax.axvline(1, color='black', linestyle='-', linewidth=1)
    ax.plot([0, 2.5], [0, 2.5], color="black", linestyle="--", linewidth=1)
    ax.set_xlim(0, 2.5); ax.set_ylim(0, 2.5)

    ax.set_xlabel(xcol); ax.set_ylabel(ycol)
    if title: ax.set_title(title)

    # Colorbars
    sm_above = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_above); sm_above.set_array([])
    sm_below = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_below); sm_below.set_array([])
    fig.colorbar(sm_below, ax=ax, location="right", fraction=0.046, pad=0.04).set_label("Relaxation")
    fig.colorbar(sm_above, ax=ax, location="left",  fraction=0.046, pad=0.10).set_label("Selection Towards Outgroup")

    # Candidate labeling
    diff_above = y - x
    diff_below = x - y

    def _pick_candidates(side_mask, diffs, topn, thr):
        loc_idx = np.where(side_mask)[0]
        if thr is not None:
            return loc_idx[diffs[loc_idx] >= thr]
        if topn is not None and loc_idx.size:
            order = np.argsort(diffs[loc_idx])[::-1]
            return loc_idx[order[:topn]]
        return loc_idx

    cand_above = _pick_candidates(above, diff_above, label_top_n_above, label_diff_min_above)
    cand_below = _pick_candidates(below, diff_below, label_top_n_below, label_diff_min_below)

    # Always label any point with x>1 or y>1
    over1 = np.where((x > 1) | (y > 1))[0]
    cand_above = np.unique(np.concatenate([cand_above, over1[y[over1] > x[over1]]]))
    cand_below = np.unique(np.concatenate([cand_below, np.setdiff1d(over1, over1[y[over1] > x[over1]])]))
    label_idx = np.unique(np.concatenate([cand_above, cand_below]))

    # Annotate & repel labels on the plot (short parentheses-only names)
    texts = []
    for li in label_idx:
        gx, gy = float(x[li]), float(y[li])
        name = _label_from_vfg(data.iloc[li][gene_col])
        ann = ax.annotate(
            name, (gx, gy),
            xytext=(gx + 0.01, gy + 0.01),
            textcoords="data",
            fontsize=label_fontsize,
            ha="left", va="bottom",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
            arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5)
        )
        texts.append(ann)
    if texts:
        adjust_text(
            texts, x=x, y=y, ax=ax,
            expand_points=(1.2, 1.2), expand_text=(1.05, 1.2),
            force_points=0.3, force_text=0.5,
            only_move={'points': 'y', 'text': 'xy'},
            autoalign=True, lim=300,
            arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5)
        )

    # ---- Save the figure ----
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Saved {out_png}")

    # ---- Also save labeled points to CSV (with FULL gene names) ----
    if label_idx.size:
        labeled = data.iloc[label_idx].copy()

        # Ensure we have a Series (not a DataFrame) even if there are duplicate column names
        gene_series = labeled[gene_col]
        if isinstance(gene_series, pd.DataFrame):
            gene_series = gene_series.iloc[:, 0]
        gene_series = gene_series.astype(str)

        # VFG ID and full name
        labeled["VFG_ID"]   = gene_series.str.extract(r"(VFG\d+)", flags=re.I, expand=False).str.upper()
        labeled["FullName"] = gene_series.apply(_full_gene_from_vfg)

        # Use the plotted numeric coordinates
        labeled[xcol] = x[label_idx]
        labeled[ycol] = y[label_idx]
        labeled["Side"]    = np.where(y[label_idx] > x[label_idx], "above", "below")
        labeled["AbsDiff"] = np.abs(y[label_idx] - x[label_idx])

        out_csv = out_png.replace(".png", "_labeled_points.csv")
        os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
        cols = ["VFG_ID", "FullName", xcol, ycol, "Side", "AbsDiff"]
        (
            labeled[cols]
            .sort_values("AbsDiff", ascending=False)
            .to_csv(out_csv, index=False)
        )
        print(f"Saved labels CSV: {out_csv}")


# -----------------------
# Scatterplots
# -----------------------
plot_pair_like_AA_AB(
    df,
    xcol="00_avg_dnds",
    ycol="01_avg_dnds",
    out_png="viz_scatter_00_vs_01_like_AA_AB.png",
    gene_col="Gene",
    title="Average dN/dS of Clade B Samples against Clade A/B References",
    label_top_n_above=15,
    label_top_n_below=10,
    label_fontsize=7,
    fig_size=(10, 7)
)

plot_pair_like_AA_AB(
    df,
    xcol="11_avg_dnds",
    ycol="10_avg_dnds",
    out_png="viz_scatter_11_vs_10_like_AA_AB.png",
    gene_col="Gene",
    title="Average dN/dS of Clade A Samples against Clade B/A References",
    label_top_n_above=15,
    label_top_n_below=10,
    label_fontsize=7,
    fig_size=(10, 7)
)
