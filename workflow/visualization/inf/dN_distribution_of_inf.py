import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, MaxNLocator, FormatStrFormatter
from matplotlib.ticker import ScalarFormatter, StrMethodFormatter

# -------- CONFIG --------
# Dir where Inf values infered from 
dir1 = "/dir/to/dnds_output/1_vs_1"
dir2 = "/dir/to/dnds_output/1_vs_0"

required_ratio_col = "dN/dS"
required_dn_col    = "dN"

fig_all_out = "dN_values_distribution_sidepanels_ALL.png"
fig_inf_out = "dN_values_distribution_sidepanels_INF.png"
# ------------------------

def read_csv_auto(path):
    return pd.read_csv(path, sep=None, engine="python")

def is_inf_value(x):
    if isinstance(x, (int, float, np.floating)): return np.isinf(x)
    if isinstance(x, str):
        s = x.strip().lower()
        if s in {"inf", "+inf", "-inf", "infinity"} or "∞" in s: return True
        try: return np.isinf(float(s))
        except Exception: return False
    return False

# -------- Collect data --------
dn1_all, dn2_all, dn1_inf, dn2_inf = [], [], [], []

for directory, out_all, out_inf in [(dir1, dn1_all, dn1_inf), (dir2, dn2_all, dn2_inf)]:
    for fname in sorted(f for f in os.listdir(directory) if f.endswith(".csv")):
        path = os.path.join(directory, fname)
        try:
            df = read_csv_auto(path)
            if required_dn_col not in df.columns: 
                continue
            # ALL rows
            dn = pd.to_numeric(df[required_dn_col], errors="coerce")
            dn = dn.replace([np.inf, -np.inf], np.nan).dropna()
            out_all.extend([v for v in dn if v >= 0])
            # INF rows only
            if required_ratio_col in df.columns:
                mask_inf = df[required_ratio_col].apply(is_inf_value)
                dn_inf = pd.to_numeric(df.loc[mask_inf, required_dn_col], errors="coerce")
                dn_inf = dn_inf.replace([np.inf, -np.inf], np.nan).dropna()
                out_inf.extend([v for v in dn_inf if v >= 0])
        except Exception:
            continue

# -------- Plot ALL rows --------
fig = plt.figure(figsize=(14, 5))
gs = GridSpec(1, 3, width_ratios=[2.2, 1.0, 1.0], wspace=0.25)

ax_main = fig.add_subplot(gs[0, 0])
ax_main.hist(dn1_all, bins=200, alpha=0.5, label="CladeA_vs_CladeA")
ax_main.hist(dn2_all, bins=200, alpha=0.5, label="CladeA_vs_CladeB")
ax_main.set_xlim(0, 2.5)
ax_main.set_yscale("log")
ax_main.set_xlabel("dN values (0–2.5)")
ax_main.set_ylabel("Count (log scale)")
ax_main.set_title("dN distribution (ALL rows)")
ax_main.legend()

ax_a = fig.add_subplot(gs[0, 1])
zoom_a = (0, 0.02)
ax_a.hist([v for v in dn1_all if zoom_a[0] <= v <= zoom_a[1]], bins=50, alpha=0.6, label="A_vs_A")
ax_a.hist([v for v in dn2_all if zoom_a[0] <= v <= zoom_a[1]], bins=50, alpha=0.6, label="A_vs_B")
ax_a.set_xlim(*zoom_a)
ax_a.set_ylim(bottom=0)
ax_a.xaxis.set_major_locator(MultipleLocator(0.001))
ax_a.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
ax_a.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
ax_a.set_xlabel("dN (0–0.02)")
ax_a.set_title("Zoom: 0–0.02")
ax_a.tick_params(axis="x", rotation=90)

ax_b = fig.add_subplot(gs[0, 2])
zoom_b = (0.90, 1.10)
ax_b.hist([v for v in dn1_all if zoom_b[0] <= v <= zoom_b[1]], bins=30, alpha=0.6, label="A_vs_A")
ax_b.hist([v for v in dn2_all if zoom_b[0] <= v <= zoom_b[1]], bins=30, alpha=0.6, label="A_vs_B")
ax_b.set_xlim(*zoom_b)
ax_b.set_ylim(bottom=0)
ax_b.xaxis.set_major_locator(MultipleLocator(0.01))
ax_b.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
ax_b.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
ax_b.set_xlabel("dN (0.90–1.10)")
ax_b.set_title("Zoom: ~1")
ax_b.tick_params(axis="x", rotation=90)

plt.tight_layout()
plt.savefig(fig_all_out, dpi=300)
plt.close()

# -------- Plot INF rows only --------
fig = plt.figure(figsize=(14, 5))
gs = GridSpec(1, 3, width_ratios=[2.2, 1.0, 1.0], wspace=0.45)

ax_main = fig.add_subplot(gs[0, 0])

# Use a common set of bins for both series (24 bins across 0–2.5)
bins_main = np.linspace(0, 2.5, 50)  # 50 edges -> 49 bins

# Single call with both datasets ensures identical bins for both
ax_main.hist([dn1_inf, dn2_inf],
             bins=bins_main, range=(0, 2.5),
             alpha=0.5, rwidth=0.9,
             label=["CladeA_vs_CladeA", "CladeA_vs_CladeB"])

ax_main.set_xlim(0, 2.5)
ax_main.set_yscale("log")
ax_main.set_xlabel("dN values (0–2.5)")
ax_main.set_ylabel("Count (log scale)")
ax_main.set_title("dN distribution (INF rows only)")
ax_main.legend()

# Keep your zooms unchanged
zoom_c = (0, 0.02)
zoom_b = (0.90, 1.10)

ax_a = fig.add_subplot(gs[0, 1])
ax_a.hist([v for v in dn1_inf if zoom_c[0] <= v <= zoom_c[1]], bins=50, alpha=0.6, label="A_vs_A")
ax_a.hist([v for v in dn2_inf if zoom_c[0] <= v <= zoom_c[1]], bins=50, alpha=0.6, label="A_vs_B")
ax_a.set_xlim(*zoom_c)
ax_a.set_ylim(bottom=0)
ax_a.xaxis.set_major_locator(MultipleLocator(0.001))
ax_a.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
ax_a.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
ax_a.set_xlabel("dN (0–0.02)")
ax_a.set_title("Zoom: 0–0.02")
ax_a.tick_params(axis="x", rotation=90)

ax_b = fig.add_subplot(gs[0, 2])
ax_b.hist([v for v in dn1_inf if zoom_b[0] <= v <= zoom_b[1]], bins=30, alpha=0.6, label="A_vs_A")
ax_b.hist([v for v in dn2_inf if zoom_b[0] <= v <= zoom_b[1]], bins=30, alpha=0.6, label="A_vs_B")
ax_b.set_xlim(*zoom_b)
ax_b.set_ylim(0, 1)
ax_b.xaxis.set_major_locator(MultipleLocator(0.01))
ax_b.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
ax_b.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
ax_b.set_xlabel("dN (0.90–1.10)")
ax_b.set_title("Zoom: ~1")
ax_b.tick_params(axis="x", rotation=90)

plt.tight_layout()
plt.savefig(fig_inf_out, dpi=300)
plt.close()


print("Saved:", fig_all_out, "and", fig_inf_out)
