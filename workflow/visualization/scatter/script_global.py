import os, re, numpy as np, pandas as pd

# this is the script used to generate overall dnds statistics used for following visualization script

folders = {
    "01": "/dir/to/dnds_output/0_vs_1",
    "00": "/dir/to/dnds_output/0_vs_0",
    "11": "/dir/to/dnds_output/1_vs_1",
    "10": "/dir/to/dnds_output/1_vs_0",
}
out_csv = "/dir/to/output/multi_folder_dnds_summary_all_inf.csv"
EXP = {"dN","dS","dN/dS"}

def read_df(p):
    for sep in (None, "\t", ",", r"\s+"):
        try:
            df = pd.read_csv(p, sep=sep, engine="python")
            df.columns = [str(c).strip() for c in df.columns]
            if EXP.issubset(df.columns): return df
        except: pass
    return None

def gene_from_fname(fn):
    m = re.search(r"(VFG\d{6})", fn, re.I)
    return m.group(1).upper() if m else fn.split(".")[0]

rows = {}
for label, d in folders.items():
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if not fn.endswith(".csv"): continue
        df = read_df(os.path.join(d, fn))
        if df is None: continue
        df["dN/dS"] = df["dN/dS"].replace(
            ["inf","Inf","INF","+inf","+Inf","+INF","-inf","-Inf","-INF"],
            [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,-np.inf,-np.inf,-np.inf]
        )
        for c in ["dN","dS","dN/dS"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.replace([np.inf,-np.inf], np.nan)
        df = df[df["dN/dS"] >= 0]
        avg = df["dN/dS"].dropna().mean()
        rows.setdefault(gene_from_fname(fn), {})[f"{label}_avg_dnds"] = avg

out = pd.DataFrame.from_dict(rows, orient="index").reset_index().rename(columns={"index":"Gene"})
cols = ["Gene","00_avg_dnds","01_avg_dnds","10_avg_dnds","11_avg_dnds"]
out = out[[c for c in cols if c in out.columns] + [c for c in out.columns if c not in cols]]
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
out.to_csv(out_csv, index=False)
print("Saved", out_csv)
