#!/usr/bin/env python3
import os, re, numpy as np, pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# ================== USER SETTINGS ==================
DIRS = {
    "00": "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/0_vs_0",  # internal B
    "01": "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/0_vs_1",  # external B (vs A)
    "11": "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_1",  # internal A
    "10": "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_0",  # external A (vs B)
}
ANNOT = "/data1/B_Li/vfdb/db/vfdb_names.csv"  # VFG_ID + Gene names

FC_THR = 1.0     # |log2FC| threshold for “meaningful” effect (e.g., 1 = 2x)
Q_THR  = 0.005   # q-value threshold for significance (dashed line)
EPS    = 1e-8    # numeric stability
LABEL_MAX = 40
# ===================================================

def gene_from_filename(fname):
    m = re.search(r"(VFG\d+)", fname, flags=re.I)
    return m.group(1).upper() if m else None

def load_dnds_dist(folder):
    """Return {VFG_ID: np.array of dN/dS}, dropping NaN/inf/negatives."""
    out = {}
    if not os.path.isdir(folder):
        return out
    for fname in os.listdir(folder):
        if not fname.endswith("_groupwise_dnds.csv"):
            continue
        gid = gene_from_filename(fname)
        if gid is None:
            continue
        path = os.path.join(folder, fname)
        try:
            df = pd.read_csv(path, sep=None, engine="python")
            if "dN/dS" not in df.columns:
                continue
            s = pd.to_numeric(df["dN/dS"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            s = s[s >= 0]
            if len(s):
                out[gid] = s.to_numpy(float)
        except Exception:
            pass
    return out

def _cohen_d(xe, xi):
    ne, ni = len(xe), len(xi)
    if ne < 2 or ni < 2: return np.nan
    me, mi = float(np.mean(xe)), float(np.mean(xi))
    ve, vi = float(np.var(xe, ddof=1)), float(np.var(xi, ddof=1))
    sp = np.sqrt(((ne-1)*ve + (ni-1)*vi) / (ne + ni - 2)) if (ne+ni-2) > 0 else np.nan
    return np.nan if (sp == 0 or np.isnan(sp)) else (me - mi) / sp

def _hedges_g(d, ne, ni):
    if np.isnan(d): return np.nan
    n = ne + ni
    J = 1 - 3.0/(4.0*n - 9.0) if n > 2 else 1.0
    return J * d

def _median_log2_ratio(xe, xi, eps=1e-8):
    me = np.median(xe) if len(xe) else np.nan
    mi = np.median(xi) if len(xi) else np.nan
    if np.isnan(me) or np.isnan(mi): return np.nan
    return np.log2((me + eps) / (mi + eps))

def _cliffs_delta(xe, xi):
    """O(n log n) implementation of Cliff's delta."""
    ne, ni = len(xe), len(xi)
    if ne == 0 or ni == 0: return np.nan
    xe = np.sort(np.asarray(xe))
    xi = np.sort(np.asarray(xi))
    i = j = 0
    wins = ties = 0
    while i < ne and j < ni:
        if xe[i] > xi[j]:
            wins += (ne - i)      # all remaining xe are > this xi[j]
            j += 1
        elif xe[i] < xi[j]:
            i += 1
        else:
            # handle tie block
            v = xe[i]
            i0, j0 = i, j
            while i < ne and xe[i] == v: i += 1
            while j < ni and xi[j] == v: j += 1
            ties += (i - i0) * (j - j0)
    losses = ne*ni - wins - ties
    return (wins - losses) / (ne*ni)


def per_gene_effects(dist_internal, dist_external, eps=EPS):
    """
    Return per-gene table with:
      n_internal, n_external, Mean_internal, Mean_external, SD_internal, SD_external,
      Delta_mean, log2FC, Cohen_d, Hedges_g, Median_log2_ratio, Cliffs_delta,
      p_value (Welch), q_value (BH)
    """
    rows = []
    genes = sorted(set(dist_internal) & set(dist_external))
    for gid in genes:
        xi, xe = dist_internal[gid], dist_external[gid]
        ni, ne = len(xi), len(xe)
        if ni < 3 or ne < 3:
            continue

        mi = float(np.mean(xi));  me = float(np.mean(xe))
        sdi = float(np.std(xi, ddof=1)); sde = float(np.std(xe, ddof=1))

        # primary effect (your plots use this)
        log2fc = np.log2((me + eps) / (mi + eps))
        delta_mean = me - mi

        # standardized effects
        d = _cohen_d(xe, xi)
        g = _hedges_g(d, ne, ni)

        # robust effects
        med_log2 = _median_log2_ratio(xe, xi)
        cd = _cliffs_delta(xe, xi)

        # test (Welch t)
        _, p = ttest_ind(xe, xi, equal_var=False, nan_policy="omit")

        rows.append({
            "VFG_ID": gid,
            "n_internal": ni, "n_external": ne,
            "Mean_internal": mi, "Mean_external": me,
            "SD_internal": sdi, "SD_external": sde,
            "Delta_mean": delta_mean,
            "log2FC": log2fc,
            "Cohen_d": d, "Hedges_g": g,
            "Median_log2_ratio": med_log2,
            "Cliffs_delta": cd,
            "p_value": p
        })

    df = pd.DataFrame(rows)
    if df.empty:
        df["q_value"] = []
        return df

    df["q_value"] = multipletests(df["p_value"].values, method="fdr_bh")[1]
    return df

def load_annotations(annot_path):
    """
    Map VFG_ID -> short name (no parentheses). Fallbacks: full 'Gene' -> VFG_ID.
    Robust across pandas versions.
    """
    try:
        annot = pd.read_csv(annot_path)
    except Exception as e:
        print(f"[WARN] failed to read annotation CSV: {e}")
        return {}

    cols = {c.lower(): c for c in annot.columns}
    col_id   = cols.get("id")
    col_gene = cols.get("gene")

    # VFG_ID as Series
    if col_id is None:
        if "VFG_ID" in annot.columns:
            vfg_df = annot["VFG_ID"].astype(str).str.extract(r"(VFG\d+)", flags=re.I, expand=True)
        else:
            print("[WARN] No 'ID' or 'VFG_ID' column in annotations; no labels mapped.")
            return {}
    else:
        vfg_df = annot[col_id].astype(str).str.extract(r"(VFG\d+)", flags=re.I, expand=True)
    vfg = vfg_df.iloc[:, 0].str.upper()

    # Label: prefer text inside (...) if present, else full Gene; no parentheses in output
    if col_gene is not None:
        paren_df = annot[col_gene].astype(str).str.extract(r"\(([^)]+)\)", expand=True)
        paren = paren_df.iloc[:, 0]  # Series (may be NaN)
        label = paren.fillna(annot[col_gene].astype(str))
    else:
        label = pd.Series([None] * len(annot))

    if isinstance(label, pd.DataFrame):
        label = label.iloc[:, 0]
    label = label.astype(str).str.strip().replace({"": np.nan})
    label = label.where(label.notna(), vfg)

    mapping = pd.Series(label.values, index=vfg)
    mapping = mapping[~mapping.index.isna()]
    mapping = mapping[~mapping.index.duplicated(keep="first")]
    return mapping.to_dict()

def volcano_color_by_fc(df, out_png, title="", fc_thr=FC_THR, q_thresh=Q_THR,
                        id_to_label=None, label_max=LABEL_MAX):
    """
    Volcano with:
      - GREEN: non-significant (q >= q_thresh), regardless of log2FC
      - RED:   significant AND log2FC >= fc_thr   (external↑)
      - BLUE:  significant AND log2FC <= -fc_thr  (internal↑)
      - GRAY:  significant BUT |log2FC| < fc_thr  (small effect)
    Labels: significant & meaningful genes with short names (no parentheses).
    """
    if df is None or df.empty:
        print(f"[WARN] empty dataframe for {title}")
        return

    x = df["log2FC"].astype(float)
    q = df["q_value"].astype(float)
    y = -np.log10(q.clip(lower=1e-300))

    sig = q < q_thresh
    ns  = ~sig
    pos = sig & (x >=  fc_thr)
    neg = sig & (x <= -fc_thr)
    neu = sig & ~(pos | neg)

    fig = plt.figure(figsize=(8.2, 5.8))
    ax = plt.gca()

    # draw points (ns first so sig draws on top)
    ax.scatter(x[ns],  y[ns],  s=16, c="#33a02c", alpha=0.85, label=f"q ≥ {q_thresh:g} (non-sig)", zorder=1)
    ax.scatter(x[neu], y[neu], s=16, c="#BFBFBF", alpha=0.85, label=f"q < {q_thresh:g}, |log2FC| < {fc_thr:g}", zorder=2)
    ax.scatter(x[neg], y[neg], s=18, c="#377eb8", alpha=0.95, label=f"q < {q_thresh:g}, ≤ -{fc_thr:g} (internal↑)", zorder=3)
    ax.scatter(x[pos], y[pos], s=18, c="#e41a1c", alpha=0.95, label=f"q < {q_thresh:g}, ≥ +{fc_thr:g} (external↑)", zorder=3)

    if q_thresh is not None:
        ax.axhline(-np.log10(q_thresh), ls="--", lw=1, c="k", alpha=0.6, zorder=0)

    ax.set_xlabel("log2 fold-change (external / internal)")
    ax.set_ylabel(r"$-\log_{10}(q)$")
    if title:
        ax.set_title(title)

    # --- Label significant & meaningful genes (no labels for green layer)
    lab_mask = sig & (np.abs(x) >= fc_thr)
    if lab_mask.any():
        lab_df = df.loc[lab_mask, ["VFG_ID", "log2FC", "q_value"]].copy()
        # prioritize by |log2FC| (desc), then by smaller q
        order_idx = np.lexsort((lab_df["q_value"].values, -np.abs(lab_df["log2FC"].values)))
        lab_df = lab_df.iloc[order_idx]

        # >>> CAP number of labels <<<
        if label_max is not None and len(lab_df) > label_max:
            lab_df = lab_df.iloc[:label_max]

        base_dx_pts, base_dy_pts = 16, 10
        for i, (vid, x0, q0) in enumerate(zip(lab_df["VFG_ID"], lab_df["log2FC"], lab_df["q_value"])):
            y0 = -np.log10(max(float(q0), 1e-300))
            name = id_to_label.get(str(vid), str(vid)) if id_to_label else str(vid)

            sx = 1 if (i % 2 == 0) else -1
            sy = 1 if (x0 >= 0) else -1
            dx_pts = int(base_dx_pts * sx * (1 + 0.3 * (i % 3)))
            dy_pts = int(base_dy_pts * sy * (1 + 0.3 * (i % 2)))

            ann = plt.annotate(
                name,
                xy=(x0, y0), xycoords="data",
                xytext=(dx_pts, dy_pts), textcoords="offset points",
                fontsize=7,
                ha="left" if sx > 0 else "right",
                va="bottom" if sy > 0 else "top",
                arrowprops=dict(
                    arrowstyle="->", lw=0.7, color="black",
                    shrinkA=0, shrinkB=4, mutation_scale=8,
                    connectionstyle="arc3"
                ),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.65),
                zorder=4, annotation_clip=True
            )
            # keep tight_layout stable
            ann.set_in_layout(False)
            if getattr(ann, "arrow_patch", None) is not None:
                ann.arrow_patch.set_in_layout(False)

    # legend
    ax.legend(
        loc="lower right",
        frameon=True, fancybox=True, framealpha=0.9,
        borderpad=0.4, labelspacing=0.4, handlelength=1.2
    )

    # robust layout: try tight_layout; if it fails, fall back to manual margins
    try:
        plt.tight_layout()
    except Exception as e:
        print(f"[WARN] tight_layout failed ({e}); using subplots_adjust fallback.")
        plt.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.12)

    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved {out_png}")


def main():
    # Load distributions
    dist_00 = load_dnds_dist(DIRS["00"])  # internal B
    dist_01 = load_dnds_dist(DIRS["01"])  # external B
    dist_11 = load_dnds_dist(DIRS["11"])  # internal A
    dist_10 = load_dnds_dist(DIRS["10"])  # external A

    # Stats
    res_B = per_gene_effects(dist_00, dist_01)
    res_A = per_gene_effects(dist_11, dist_10)

    # Load short-name mapping (no parentheses)
    vfg_to_short = load_annotations(ANNOT)

    # Optional: save
    if not res_B.empty: res_B.to_csv("effects_cladeB.csv", index=False)
    if not res_A.empty: res_A.to_csv("effects_cladeA.csv", index=False)

    # Plots
    volcano_color_by_fc(
        res_B, "volcano_cladeB_fcdir.png",
        title="Clade B (volcano: q-significance)",
        fc_thr=FC_THR, q_thresh=Q_THR, id_to_label=vfg_to_short
    )
    volcano_color_by_fc(
        res_A, "volcano_cladeA_fcdir.png",
        title="Clade A (volcano: q-significance)",
        fc_thr=FC_THR, q_thresh=Q_THR, id_to_label=vfg_to_short
    )

if __name__ == "__main__":
    main()
