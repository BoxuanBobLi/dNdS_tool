#!/usr/bin/env python3
import os, re, numpy as np, pandas as pd
import matplotlib.pyplot as plt

# ================== USER SETTINGS ==================
CSV_A = "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/visualization/scatter_fdr/effects_cladeA.csv"
CSV_B = "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/visualization/scatter_fdr/effects_cladeB.csv"
ANNOT = "/data1/B_Li/vfdb/db/vfdb_names.csv"

FC_THR      = 1.0      # label if |log2FC| >= FC_THR  (e.g., 1.0 = 2x)
Q_THR       = 0.05     # label if q_value < Q_THR
VLINE_AT    = 1.0      # vertical guide at average dN/dS = 1.0; set None to remove
LABEL_MAX   = 30       # cap #labels; set None to label all passing thresholds
SAVE_TABLES = False     # we’re reading from CSVs, but keep in case you tweak
# Coloring by q:
Q_CLIP_MIN  = 1e-300   # avoid inf in -log10(q)
Q_NEGLOG_VMAX = 50.0   # colorbar upper cap for -log10(q); tweak if you want (e.g., 30–100)
CMAP        = "viridis"  # any matplotlib cmap

ABS_A_LABEL_THR = 1.0   # label genes with average dN/dS (A_raw) > 1
LABEL_MAX_ABS   = 20    # cap for the extra A_raw>1 labels (set None for no cap)

# ===================================================


# ---------- helpers ----------
def load_annotations(annot_path):
    """
    Map VFG_ID -> short name (no parentheses).
    Fallbacks: full 'Gene' -> VFG_ID. Robust to pandas versions & missing cols.
    """
    try:
        annot = pd.read_csv(annot_path)
    except Exception as e:
        print(f"[WARN] failed to read annotation CSV: {e}")
        return {}

    cols = {c.lower(): c for c in annot.columns}
    col_id   = cols.get("id")
    col_gene = cols.get("gene")

    # ---- VFG_ID as a Series
    if col_id is None:
        if "VFG_ID" in annot.columns:
            vfg_df = annot["VFG_ID"].astype(str).str.extract(r"(VFG\d+)", flags=re.I, expand=True)
        else:
            print("[WARN] No 'ID' or 'VFG_ID' column found in annotations; no labels mapped.")
            return {}
    else:
        vfg_df = annot[col_id].astype(str).str.extract(r"(VFG\d+)", flags=re.I, expand=True)
    vfg = vfg_df.iloc[:, 0].str.upper()

    # ---- Label (prefer text inside parentheses, else whole Gene), as a Series
    if col_gene is not None:
        paren_df = annot[col_gene].astype(str).str.extract(r"\(([^)]+)\)", expand=True)
        paren = paren_df.iloc[:, 0]  # Series (may be all NaN)
        label = paren.fillna(annot[col_gene].astype(str))
    else:
        label = pd.Series([None] * len(annot))

    # Coerce to Series if any op above yielded a DataFrame on your pandas
    if isinstance(label, pd.DataFrame):
        label = label.iloc[:, 0]

    # Clean empties, then fallback to VFG_ID
    label = label.astype(str).str.strip()
    label = label.replace({"": np.nan})
    label = label.where(label.notna(), vfg)

    # Build mapping, drop rows with missing VFG_ID
    mapping = pd.Series(label.values, index=vfg)
    mapping = mapping[~mapping.index.isna()]
    mapping = mapping[~mapping.index.duplicated(keep="first")]

    return mapping.to_dict()



def load_effects(csv_path):
    """
    Load effects CSV with columns:
    VFG_ID, Mean_internal, Mean_external, log2FC, p_value, q_value
    """
    df = pd.read_csv(csv_path)
    # ensure numeric
    for col in ["Mean_internal", "Mean_external", "log2FC", "p_value", "q_value"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["VFG_ID"] = df["VFG_ID"].astype(str)
    # Compute A (raw) and M for plotting
    df["A_raw"] = 0.5 * (df["Mean_internal"] + df["Mean_external"])
    df["M_log2FC"] = df["log2FC"]
    # color metric: -log10(q), clipped
    df["_q_clipped"] = np.maximum(df["q_value"].values, Q_CLIP_MIN)
    df["_neglog10q"] = -np.log10(df["_q_clipped"])
    return df


def ma_plot_qcolor(
    df, out_png, title="",
    fc_thr=FC_THR, q_thr=Q_THR,
    vline_at=VLINE_AT,
    label_max=LABEL_MAX,
    id_to_label=None,
    cmap=CMAP, neglog_vmax=Q_NEGLOG_VMAX
):
    """
    MA plot without logging dN/dS, colored by -log10(q) from CSV.

      X (A_raw) = 0.5 * (Mean_external + Mean_internal)
      Y (M)     = log2FC (external / internal)
      Color     = -log10(q) (clipped), capped at neglog_vmax

    Labels: genes with (q < q_thr) & (|log2FC| >= fc_thr) using short names.
    """
    if df is None or df.empty:
        print(f"[WARN] empty dataframe for {title}")
        return

    # Plot data
    A = df["A_raw"].astype(float).values
    M = df["M_log2FC"].astype(float).values
    C_raw = df["_neglog10q"].astype(float).values
    C = np.minimum(C_raw, float(neglog_vmax))  # cap for display

    fig = plt.figure(figsize=(9, 5.8))
    ax = plt.gca()

    sc = ax.scatter(
        A, M, s=18, c=C, cmap=cmap, vmin=0.0, vmax=float(neglog_vmax),
        alpha=0.9, edgecolors="none", zorder=2
    )

    # Guides
    ax.axhline(0, ls="--", lw=1, c="k", alpha=0.6, zorder=1)  # no change
    if vline_at is not None:
        ax.axvline(vline_at, ls=":", lw=1, c="k", alpha=0.5, zorder=1)

    # Colorbar
    cbar = plt.colorbar(sc, pad=0.1)
    cbar.set_label(r"$-\log_{10}(q)$")

    ax.set_xlabel("A = average dN/dS (raw)")
    ax.set_ylabel("M = log2 fold-change (external / internal)")
    if title:
        ax.set_title(title)

    # -------- Labels pass #1: significant & meaningful (q < q_thr AND |M| >= fc_thr)
    labeled_ids = set()

    sig_mask = (df["q_value"].values < q_thr) & (np.abs(M) >= fc_thr)
    if np.any(sig_mask):
        lab_df = df.loc[sig_mask, ["VFG_ID", "A_raw", "M_log2FC", "q_value"]].copy()
        # prioritize by |M| desc, then smaller q; drop duplicate VFG_IDs
        order_idx = np.lexsort((lab_df["q_value"].values, -np.abs(lab_df["M_log2FC"].values)))
        lab_df = lab_df.iloc[order_idx].drop_duplicates("VFG_ID", keep="first")
        if label_max is not None and len(lab_df) > label_max:
            lab_df = lab_df.iloc[:label_max]

        base_dx_pts, base_dy_pts = 18, 10
        for i, r in lab_df.iterrows():
            x0 = float(r["A_raw"]); y0 = float(r["M_log2FC"])
            vid = str(r["VFG_ID"]); labeled_ids.add(vid)
            txt = id_to_label.get(vid, vid) if id_to_label else vid

            sx = 1 if (i % 2 == 0) else -1
            sy = 1 if (y0 >= 0) else -1
            dx_pts = int(base_dx_pts * sx * (1 + 0.3 * (i % 3)))
            dy_pts = int(base_dy_pts * sy * (1 + 0.3 * (i % 2)))

            ann = ax.annotate(
                txt,
                xy=(x0, y0), xycoords="data",
                xytext=(dx_pts, dy_pts), textcoords="offset points",
                fontsize=7,
                ha="left" if sx > 0 else "right",
                va="bottom" if sy > 0 else "top",
                arrowprops=dict(
                    arrowstyle="->", lw=0.7, color="black",
                    shrinkA=0, shrinkB=4, mutation_scale=8,
                    connectionstyle="arc3,rad=0.0"
                ),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65),
                zorder=3, annotation_clip=False
            )
            ann.set_in_layout(False)
            if getattr(ann, "arrow_patch", None) is not None:
                ann.arrow_patch.set_in_layout(False)
                ann.arrow_patch.set_clip_on(False)
            ann.set_clip_on(False)

    # -------- Labels pass #2: average dN/dS > 1, excluding already-labeled
    high_mask = (df["A_raw"].values > ABS_A_LABEL_THR)
    if np.any(high_mask):
        lab2 = df.loc[high_mask, ["VFG_ID", "A_raw", "M_log2FC"]].copy()
        lab2["VFG_ID"] = lab2["VFG_ID"].astype(str)
        lab2 = lab2[~lab2["VFG_ID"].isin(labeled_ids)]                # exclude already labeled
        lab2 = lab2.sort_values("A_raw", ascending=False)             # prioritize higher A
        lab2 = lab2.drop_duplicates("VFG_ID", keep="first")           # no dup IDs
        if LABEL_MAX_ABS is not None and len(lab2) > LABEL_MAX_ABS:
            lab2 = lab2.iloc[:LABEL_MAX_ABS]

        base_dx_pts2, base_dy_pts2 = 16, 10
        for j, r in lab2.iterrows():
            x0 = float(r["A_raw"]); y0 = float(r["M_log2FC"])
            vid = str(r["VFG_ID"]); labeled_ids.add(vid)
            txt = id_to_label.get(vid, vid) if id_to_label else vid

            sx = -1 if (j % 2 == 0) else 1
            sy = 1 if (y0 >= 0) else -1
            dx_pts = int(base_dx_pts2 * sx * (1 + 0.25 * (j % 3)))
            dy_pts = int(base_dy_pts2 * sy * (1 + 0.25 * (j % 2)))

            ann = ax.annotate(
                txt,
                xy=(x0, y0), xycoords="data",
                xytext=(dx_pts, dy_pts), textcoords="offset points",
                fontsize=7,
                ha="left" if sx > 0 else "right",
                va="bottom" if sy > 0 else "top",
                arrowprops=dict(
                    arrowstyle="->", lw=0.7, color="black",
                    shrinkA=0, shrinkB=4, mutation_scale=8,
                    connectionstyle="arc3,rad=0.0"
                ),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65),
                zorder=3, annotation_clip=False
            )
            ann.set_in_layout(False)
            if getattr(ann, "arrow_patch", None) is not None:
                ann.arrow_patch.set_in_layout(False)
                ann.arrow_patch.set_clip_on(False)
            ann.set_clip_on(False)

    # Layout & save (with fallback)
    try:
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
    except Exception as e:
        print(f"[WARN] savefig failed due to arrows ({e}). Retrying without tight_layout…")
        plt.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.12)
        plt.savefig(out_png, dpi=300)

    plt.close(fig)
    print(f"Saved {out_png}")


def load_full_gene_map(annot_path):
    """
    Map VFG_ID -> full Gene string (no shortening to parentheses).
    Case/column-name robust; falls back to VFG_ID when missing.
    """
    try:
        annot = pd.read_csv(annot_path)
    except Exception as e:
        print(f"[WARN] failed to read annotation CSV for full names: {e}")
        return {}

    cols = {c.lower(): c for c in annot.columns}
    col_id   = cols.get("id") or ("VFG_ID" if "VFG_ID" in annot.columns else None)
    col_gene = cols.get("gene")

    if col_id is None or col_gene is None:
        print("[WARN] Annotation CSV missing required columns ('ID'/'VFG_ID' and 'Gene'); full-name mapping will be empty.")
        return {}

    vfg = annot[col_id].astype(str).str.extract(r"(VFG\d+)", flags=re.I, expand=True).iloc[:, 0].str.upper()
    gene_full = annot[col_gene].astype(str)

    mapping = pd.Series(gene_full.values, index=vfg)
    mapping = mapping[~mapping.index.isna()]
    mapping = mapping[~mapping.index.duplicated(keep="first")]
    return mapping.to_dict()


def write_threshold_csv(df, vfg_to_fullname, out_csv, a_thr=1.0, fc_thr=2.0):
    """
    Save rows passing A_raw >= a_thr OR |log2FC| >= fc_thr, with full gene names.
    """
    if df is None or df.empty:
        print(f"[WARN] empty dataframe; skip writing {out_csv}")
        return

    # OR condition (your new requirement)
    mask = (df["A_raw"].values >= a_thr) | (np.abs(df["M_log2FC"].values) >= fc_thr)

    sub = df.loc[mask, ["VFG_ID", "A_raw", "M_log2FC", "p_value", "q_value"]].copy()
    if sub.empty:
        print(f"[INFO] no genes passed thresholds for {out_csv}")
        sub = pd.DataFrame(columns=["VFG_ID", "Gene_full", "A_raw", "M_log2FC", "p_value", "q_value"])

    sub["Gene_full"] = sub["VFG_ID"].map(vfg_to_fullname).fillna(sub["VFG_ID"])
    sub = sub[["VFG_ID", "Gene_full", "A_raw", "M_log2FC", "p_value", "q_value"]]
    sub.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")



# ---------- main ----------
def main():
    # Load short-name mapping (no parentheses)
    vfg_to_short = load_annotations(ANNOT)

    # Read precomputed effects (A and B)
    dfA = load_effects(CSV_A)
    dfB = load_effects(CSV_B)

    # Optional: save after enrichment
    if SAVE_TABLES:
        dfA.to_csv("effects_cladeA_enriched.csv", index=False)
        dfB.to_csv("effects_cladeB_enriched.csv", index=False)

    # Plots colored by -log10(q)
    ma_plot_qcolor(
        dfA, "MA_cladeA_qcolor.png",
        title="Clade A — MA plot",
        fc_thr=FC_THR, q_thr=Q_THR, vline_at=VLINE_AT,
        label_max=LABEL_MAX, id_to_label=vfg_to_short,
        cmap=CMAP, neglog_vmax=Q_NEGLOG_VMAX
    )
    ma_plot_qcolor(
        dfB, "MA_cladeB_qcolor.png",
        title="Clade B — MA plot",
        fc_thr=FC_THR, q_thr=Q_THR, vline_at=VLINE_AT,
        label_max=LABEL_MAX, id_to_label=vfg_to_short,
        cmap=CMAP, neglog_vmax=Q_NEGLOG_VMAX
    )

        # ---- Extra outputs: full gene names for threshold-passers (does not affect original outputs)
    vfg_to_fullname = load_full_gene_map(ANNOT)
    write_threshold_csv(dfA, vfg_to_fullname, "genes_pass_threshold_cladeA.csv", a_thr=1.0, fc_thr=2.0)
    write_threshold_csv(dfB, vfg_to_fullname, "genes_pass_threshold_cladeB.csv", a_thr=1.0, fc_thr=2.0)


if __name__ == "__main__":
    main()
