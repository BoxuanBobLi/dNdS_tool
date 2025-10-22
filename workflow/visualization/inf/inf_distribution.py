import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# -------- CONFIG --------
dir1 = "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_1"  # CladeA vs CladeA
dir2 = "/data1/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/dnds_output/1_vs_0"  # CladeA vs CladeB
out_summary_csv = "inf_counts_summary.csv"
out_hist_png = "inf_counts_histogram.png"

required_col = "dN/dS"   # exact column name
# ------------------------

_vfg_re = re.compile(r"(VFG\d+)", re.I)

def vfg_id_from_name(name: str):
    m = _vfg_re.search(name)
    return m.group(1).upper() if m else None

def collect_by_gene(directory):
    """
    Return {VFG_ID: filepath} for all *.csv in directory.
    If a VFG appears multiple times, keep the longest filename (usually most specific).
    """
    mapping = {}
    if not os.path.isdir(directory):
        print(f"[WARN] Directory not found: {directory}")
        return mapping
    for fname in os.listdir(directory):
        if not fname.endswith(".csv"):
            continue
        gid = vfg_id_from_name(fname)
        if not gid:
            continue
        fpath = os.path.join(directory, fname)
        if gid not in mapping or len(fname) > len(os.path.basename(mapping[gid])):
            mapping[gid] = fpath
    return mapping

def read_csv_file(path):
    """Read CSV/TSV with auto delimiter detection."""
    return pd.read_csv(path, sep=None, engine="python")

def count_infs(series: pd.Series) -> int:
    """
    Count numeric and string forms of infinity in a pandas Series.
    """
    # 1) normalize common string forms
    s = series.replace(
        ["inf", "Inf", "INF", "+inf", "-inf", "+Inf", "-Inf", "+INF", "-INF", "infinity", "Infinity", "∞"],
        np.nan  # temporary; we'll coerce numeric next and catch infs further below
    )

    # 2) try to coerce to numeric (strings that are 'inf' will become NaN here, so handle separately)
    s_num = pd.to_numeric(s, errors="coerce")

    # 3) numeric inf count
    # if original series had actual np.inf values, they survived in s_num
    numeric_inf_count = int(np.isinf(s_num).sum())

    # 4) explicit parse on original strings to catch any odd cases (e.g., "  inf ", " +inf ")
    def is_inf_text(x):
        if isinstance(x, str):
            xx = x.strip().lower()
            if xx in {"inf", "+inf", "-inf", "infinity"} or "∞" in xx:
                return True
            try:
                return np.isinf(float(xx))
            except Exception:
                return False
        return False

    str_inf_count = int(series.apply(is_inf_text).sum())

    # max of both to avoid double-counting (usually one will capture it)
    return max(numeric_inf_count, str_inf_count)

# ---------- Pair by VFG ID ----------
map1 = collect_by_gene(dir1)
map2 = collect_by_gene(dir2)

print(f"[DEBUG] dir1 files: {len(os.listdir(dir1)) if os.path.isdir(dir1) else 'NA'} "
      f"-> mapped genes: {len(map1)}")
print(f"[DEBUG] dir2 files: {len(os.listdir(dir2)) if os.path.isdir(dir2) else 'NA'} "
      f"-> mapped genes: {len(map2)}")

common_genes = sorted(set(map1.keys()) & set(map2.keys()))
print(f"[DEBUG] common VFG IDs: {len(common_genes)}")
if len(common_genes) < 5:
    print("[DEBUG] A few examples of common VFG IDs:", common_genes[:5])

inf_counts = defaultdict(dict)
missing_col = []
errors = []
missing_side = []

for gid in common_genes:
    p1 = map1.get(gid)
    p2 = map2.get(gid)
    if not p1 or not p2:
        missing_side.append((gid, bool(p1), bool(p2)))
        continue

    try:
        df1 = read_csv_file(p1)
        df2 = read_csv_file(p2)

        # Verify required column
        if required_col not in df1.columns or required_col not in df2.columns:
            missing_col.append((gid, p1, list(df1.columns), p2, list(df2.columns)))
            continue

        c1 = count_infs(df1[required_col])
        c2 = count_infs(df2[required_col])

        inf_counts[gid]["CladeA_vs_CladeA"] = c1
        inf_counts[gid]["CladeA_vs_CladeB"] = c2

    except Exception as e:
        errors.append((gid, p1, p2, str(e)))

# ---------- Build summary ----------
counts_df = pd.DataFrame.from_dict(inf_counts, orient="index").fillna(0).astype(int)
counts_df.index.name = "Gene"
counts_df = counts_df.sort_index()
counts_df.to_csv(out_summary_csv)

print(f"[DEBUG] Summary rows written: {len(counts_df)} to {out_summary_csv}")

# ---------- Histogram ----------
plt.figure(figsize=(10, 6))
if "CladeA_vs_CladeA" in counts_df.columns:
    plt.hist(counts_df["CladeA_vs_CladeA"], bins=100, alpha=0.5, label="CladeA_vs_CladeA")
if "CladeA_vs_CladeB" in counts_df.columns:
    plt.hist(counts_df["CladeA_vs_CladeB"], bins=100, alpha=0.5, label="CladeA_vs_CladeB")
plt.xlabel("Count of 'inf' dN/dS values in one VF gene")
plt.ylabel("Frequency of VF genes with this count of 'inf'")
plt.title("Distribution of 'inf' counts across Clade A genomes")
plt.legend()
plt.tight_layout()
plt.savefig(out_hist_png, dpi=300)
plt.close()

# ---------- Logs ----------
if missing_col:
    with open("missing_dNdS_files.log", "w") as fh:
        for gid, p1, cols1, p2, cols2 in missing_col:
            fh.write(f"{gid}\n  dir1: {p1}\n  cols1: {cols1}\n  dir2: {p2}\n  cols2: {cols2}\n\n")
if missing_side:
    with open("missing_gene_side.log", "w") as fh:
        for gid, has1, has2 in missing_side:
            fh.write(f"{gid}: present_in_dir1={has1}, present_in_dir2={has2}\n")
if errors:
    with open("read_errors.log", "w") as fh:
        for gid, p1, p2, err in errors:
            fh.write(f"{gid}\n  dir1: {p1}\n  dir2: {p2}\n  error: {err}\n\n")

print(f"Done.\nSummary: {out_summary_csv}\nFigure: {out_hist_png}")
if missing_col:
    print(f"{len(missing_col)} genes skipped (missing '{required_col}'). See missing_dNdS_files.log")
if missing_side:
    print(f"{len(missing_side)} genes missing in one side. See missing_gene_side.log")
if errors:
    print(f"{len(errors)} genes had errors. See read_errors.log")
