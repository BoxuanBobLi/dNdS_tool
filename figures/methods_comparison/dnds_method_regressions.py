#!/usr/bin/env python3
"""
Pairwise regressions of per-gene mean dN/dS across three methods, for four conditions.

Fix for your case:
- Filenames contain VFGxxxxx or vfgxxxxx.
- We extract gene ID case-insensitively and normalize to lowercase,
  so merges work across methods.

Outputs:
- plots/*.png   (12 plots if enough shared genes)
- regression_summary.csv
"""

import argparse
import os
import re
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


CONDITIONS_DEFAULT = ["0_vs_0", "0_vs_1", "1_vs_0", "1_vs_1"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Pairwise regressions of per-gene mean dN/dS across methods.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--aggregated", required=True, help="Aggregated/concatenated method dnds_output directory")
    p.add_argument("--consensus", required=True, help="Consensus-average method dnds_output directory")
    p.add_argument("--ancestral", required=True, help="Ancestral-average method directory")
    p.add_argument("--outdir", required=True, help="Output directory for plots + summary CSV")
    p.add_argument("--PAO1PA14reference", required=True, help="PAO1+PA14 reference method directory")

    p.add_argument("--conditions", nargs="+", default=CONDITIONS_DEFAULT, help="Condition subfolders to analyze")

    # IMPORTANT FIX: case-insensitive VFG/vfg
    p.add_argument(
        "--gene-regex",
        default=r"(?i)(vfg\d+)",
        help="Regex to extract gene id from filename (first capture group used). Case-insensitive by default.",
    )

    p.add_argument("--drop-nonfinite", action="store_true", help="Drop rows where x or y is inf/-inf before plotting")
    p.add_argument("--drop-negative", action="store_true", help="Drop rows where x or y is negative before plotting")
    p.add_argument("--min-pairs", type=int, default=20, help="Minimum # shared genes required to generate a plot")

    p.add_argument(
        "--max-genes",
        type=int,
        default=None,
        help="If set, randomly subsample this many shared genes per plot (speed/legibility).",
    )
    p.add_argument("--seed", type=int, default=0, help="Random seed used if --max-genes is set")

    # Debugging helpers
    p.add_argument(
        "--debug-counts",
        action="store_true",
        help="Print per-method/per-condition gene counts and example gene IDs.",
    )

    p.add_argument(
        "--log10",
        action="store_true",
        help="Use log10(x + eps) and log10(y + eps)",
    )
    p.add_argument(
        "--log-eps",
        type=float,
        default=1e-6,
        help="Small pseudocount added before log10",
    )


    return p.parse_args()


def safe_mkdir(path: str):
    os.makedirs(path, exist_ok=True)


def extract_gene_id(filename: str, gene_re: re.Pattern) -> str | None:
    """
    Extract gene ID like VFG000115 or vfg000115, normalize to lowercase.
    """
    m = gene_re.search(filename)
    if not m:
        return None
    return m.group(1).lower()  # normalize so all methods share the same key (vfg000115)


def read_gene_mean_dnds(csv_path: str) -> float | None:
    """
    Compute per-gene mean dN/dS.
    Drops rows with inf/-inf/NaN, but does NOT drop the gene
    unless all rows are invalid.
    """
    try:
        df = pd.read_csv(csv_path, sep=None, engine="python")
    except Exception:
        return None

    # normalize column names
    df.columns = [c.strip() for c in df.columns]

    if "dN/dS" not in df.columns:
        return None

    vals = pd.to_numeric(df["dN/dS"], errors="coerce")

    # KEY LINE: drop bad rows, not the gene
    vals = vals.replace([np.inf, -np.inf], np.nan)

    if vals.dropna().empty:
        # all rows invalid → drop gene
        return None

    return float(vals.mean(skipna=True))




def load_method_condition(method_dir: str, condition: str, gene_re: re.Pattern) -> pd.DataFrame:
    """
    Reads all *.csv in method_dir/condition; returns DataFrame with columns:
      gene, dnds
    """
    cond_dir = os.path.join(method_dir, condition)
    if not os.path.isdir(cond_dir):
        return pd.DataFrame(columns=["gene", "dnds"])

    rows = []
    for fn in os.listdir(cond_dir):
        if not fn.endswith(".csv"):
            continue

        gene = extract_gene_id(fn, gene_re)
        if gene is None:
            continue

        path = os.path.join(cond_dir, fn)
        mean_dnds = read_gene_mean_dnds(path)
        if mean_dnds is None:
            continue

        rows.append((gene, mean_dnds))

    return pd.DataFrame(rows, columns=["gene", "dnds"])


def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    return pearson_r(rx, ry)


def fit_line(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """
    OLS fit: y = slope*x + intercept
    """
    if len(x) < 2:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def clean_pairs(df, drop_nonfinite, drop_negative, log10, log_eps):
    out = df.copy()

    # ensure numeric
    out["x"] = pd.to_numeric(out["x"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")

    # drop NaN rows (regression cannot use them)
    out = out.dropna(subset=["x", "y"])

    if drop_nonfinite:
        mask = np.isfinite(out["x"].to_numpy()) & np.isfinite(out["y"].to_numpy())
        out = out[mask]

    if log10:
        # allow zeros, just not negative
        out = out[(out["x"] >= 0) & (out["y"] >= 0)]

        # LOG WITH SMALL EPSILON
        out["x"] = np.log10(out["x"] + log_eps)
        out["y"] = np.log10(out["y"] + log_eps)

    elif drop_negative:
        out = out[(out["x"] >= 0) & (out["y"] >= 0)]

    return out




def make_plot(df_xy: pd.DataFrame, xlab: str, ylab: str, title: str, out_png: str):
    x = df_xy["x"].to_numpy()
    y = df_xy["y"].to_numpy()

    r = pearson_r(x, y)
    rho = spearman_rho(x, y)
    slope, intercept = fit_line(x, y)

    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    ymin = float(np.nanmin(y))
    ymax = float(np.nanmax(y))
    lo = min(xmin, ymin)
    hi = max(xmax, ymax)

    xs = np.linspace(xmin, xmax, 200)
    ys = slope * xs + intercept

    plt.figure()
    plt.scatter(x, y)              # no explicit colors
    stats_text = (
    f"N = {len(df_xy)}\n"
    f"Pearson r = {r:.3f}\n"
    f"Spearman ρ = {rho:.3f}\n"
    f"Slope = {slope:.3g}\n"
    f"Intercept = {intercept:.3g}"
    )

    plt.text(
        0.05, 0.95,
        stats_text,
        transform=plt.gca().transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", alpha=0.8)
    )

    plt.plot(xs, ys)               # OLS line
    plt.plot([lo, hi], [lo, hi], linestyle="--")  # y=x reference

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    return r, rho, slope, intercept


def main():
    args = parse_args()
    safe_mkdir(args.outdir)

    gene_re = re.compile(args.gene_regex)

    METHODS = {
        "concatenated": args.aggregated,
        "consensus_avg": args.consensus,
        "ancestral_avg": args.ancestral,
        "PAO1PA14_ref": args.PAO1PA14reference,
    }

    data = {m: {} for m in METHODS}
    for m, base in METHODS.items():
        for cond in args.conditions:
            data[m][cond] = load_method_condition(base, cond, gene_re)

    if args.debug_counts:
        print("\n[DEBUG] per-method/per-condition gene counts (after mean dN/dS extraction):")
        for cond in args.conditions:
            print(f"  Condition: {cond}")
            for m in METHODS:
                df = data[m][cond]
                ex = ", ".join(df["gene"].head(5).tolist()) if not df.empty else ""
                print(f"    {m:14s} n={len(df):6d}  ex=[{ex}]")
        print()

    plot_dir = os.path.join(args.outdir, "plots")
    safe_mkdir(plot_dir)

    summary_rows = []
    rng = np.random.default_rng(args.seed)
    method_pairs = list(combinations(METHODS.keys(), 2))

    for cond in args.conditions:
        for m1, m2 in method_pairs:
            df1 = data[m1][cond].rename(columns={"dnds": "x"})
            df2 = data[m2][cond].rename(columns={"dnds": "y"})

            merged = df1.merge(df2, on="gene", how="inner")
            merged = clean_pairs(
                merged,
                args.drop_nonfinite,
                args.drop_negative,
                args.log10,
                args.log_eps,
            )



            if args.max_genes is not None and len(merged) > args.max_genes:
                take = rng.choice(len(merged), size=args.max_genes, replace=False)
                merged = merged.iloc[take].copy()

            if len(merged) < args.min_pairs:
                summary_rows.append(
                    {
                        "condition": cond,
                        "x_method": m1,
                        "y_method": m2,
                        "n_shared_genes": len(merged),
                        "pearson_r": np.nan,
                        "spearman_rho": np.nan,
                        "slope": np.nan,
                        "intercept": np.nan,
                        "plot_png": "",
                        "note": f"SKIPPED: < min_pairs ({args.min_pairs})",
                    }
                )
                continue

            out_png = os.path.join(plot_dir, f"{cond}__{m1}__vs__{m2}.png")
            title = f"{cond}: {m1} vs {m2}" + (" (log10)" if args.log10 else "")
            if args.log10:
                xlab = f"log10({m1} mean dN/dS + {args.log_eps:g})"
                ylab = f"log10({m2} mean dN/dS + {args.log_eps:g})"
            else:
                xlab = f"{m1} mean dN/dS"
                ylab = f"{m2} mean dN/dS"



            r, rho, slope, intercept = make_plot(merged, xlab, ylab, title, out_png)

            summary_rows.append(
                {
                    "condition": cond,
                    "x_method": m1,
                    "y_method": m2,
                    "n_shared_genes": len(merged),
                    "pearson_r": r,
                    "spearman_rho": rho,
                    "slope": slope,
                    "intercept": intercept,
                    "plot_png": out_png,
                    "note": "",
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(args.outdir, "regression_summary.csv")
    summary_df.to_csv(summary_csv, index=False)

    print(f"[OK] Wrote plots to: {plot_dir}")
    print(f"[OK] Wrote summary to: {summary_csv}")


if __name__ == "__main__":
    main()
