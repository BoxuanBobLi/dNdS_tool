#!/usr/bin/env python3

import os
import re
import time
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO

# ==========================
# EDIT THESE
# ==========================
INCSV = "/data2/B_Li/vfdb/metadata/core_VF_labeled.csv"
OUTCSV = "core_VF_labeled_with_locus_tag.csv"

# PAO1 protein FASTA and GFF downloaded from NCBI
PAO1_PROTEIN_FASTA = "GCF_000006765.1_ASM676v1_protein.faa"
PAO1_GFF = "GCF_000006765.1_ASM676v1_genomic.gff"

# Temporary files
TMP_QUERY_FASTA = "unmapped_vfdb_proteins_for_blast.faa"
BLAST_DB_PREFIX = "PAO1_protein_blastdb"
BLAST_OUT = "unmapped_vs_PAO1_blastp.tsv"

# NCBI requires an email for Entrez requests
Entrez.email = "boxuanl4@uci.edu"

# Optional but recommended if you have one
# Entrez.api_key = "YOUR_NCBI_API_KEY"

# BLAST cutoffs
MIN_PIDENT = 80.0
MIN_QCOV = 80.0
MAX_EVALUE = 1e-10

# Entrez settings
BATCH_SIZE = 50
SLEEP_SEC = 0.34


# ==========================
# Helpers
# ==========================
def extract_protein_accession(vfdb_accession):
    """
    Example:
    VFG000122(gb|NP_252230) -> NP_252230
    VFGxxxxx(gb|WP_012614635.1) -> WP_012614635.1
    """
    if pd.isna(vfdb_accession):
        return None

    s = str(vfdb_accession)

    m = re.search(r"gb\|([^)\s]+)", s)
    if m:
        return m.group(1)

    # fallback: catch NP_252230 / WP_012614635.1 / YP_XXXX etc.
    m = re.search(r"\b([A-Z]{2,3}_\d+(?:\.\d+)?)\b", s)
    if m:
        return m.group(1)

    return None


def accession_no_version(acc):
    if pd.isna(acc) or acc is None:
        return None
    return str(acc).split(".")[0]


def is_valid_text(x):
    return pd.notna(x) and str(x).strip() != ""


def extract_pao1_locus_tag_from_header(header):
    """
    Try to extract PAO1 locus tag from different FASTA header formats.
    This is only a backup. The GFF parser is more reliable.
    """
    if header is None:
        return None

    h = str(header)

    patterns = [
        r"\[locus_tag=(PA\d+)\]",
        r"locus_tag=(PA\d+)",
        r"\[old_locus_tag=(PA\d+)\]",
        r"old_locus_tag=(PA\d+)",
        r"\b(PA\d{4,5})\b",
    ]

    for pat in patterns:
        m = re.search(pat, h)
        if m:
            return m.group(1)

    return None


def parse_pao1_fasta_headers(pao1_fasta):
    """
    Build mapping:
    PAO1 protein FASTA record id -> header and possible PA locus tag.

    In many NCBI protein FASTA files, the header has NP_ accession but not PA####,
    so this alone may not recover locus tags.
    """
    mapping = {}

    for record in SeqIO.parse(pao1_fasta, "fasta"):
        subject_id_full = record.id
        subject_id_base = accession_no_version(record.id)
        header = record.description
        locus_tag = extract_pao1_locus_tag_from_header(header)

        value = {
            "blast_hit_locus_tag_from_header": locus_tag,
            "blast_hit_header": header,
        }

        mapping[subject_id_full] = value
        mapping[subject_id_base] = value

    return mapping


def parse_pao1_gff(gff_file):
    """
    Build mapping:
    PAO1 protein_id -> PAO1 locus_tag

    Example:
    NP_252230.1 -> PA3540
    NP_252230   -> PA3540

    This is the important part that fixes cases where BLAST finds NP_252230.1
    but the FASTA header does not contain PA3540.
    """
    mapping = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            feature_type = parts[2]
            attrs = parts[8]

            if feature_type != "CDS":
                continue

            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v

            locus_tag = attr_dict.get("locus_tag")
            protein_id = attr_dict.get("protein_id")

            if locus_tag and protein_id:
                mapping[protein_id] = locus_tag
                mapping[accession_no_version(protein_id)] = locus_tag

    print(f"Parsed PAO1 GFF protein_id -> locus_tag mappings: {len(mapping)}")
    return mapping


def fetch_ncbi_records(accessions, batch_size=50, sleep_sec=0.34):
    """
    Batch fetch protein GenBank records from NCBI.

    Returns:
    accession_no_version -> direct annotation + protein sequence
    """
    mapping = {}

    accessions = [a for a in accessions if pd.notna(a)]
    accessions = sorted(set(accessions))

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        print(f"Fetching {i+1}-{i+len(batch)} / {len(accessions)}")

        try:
            handle = Entrez.efetch(
                db="protein",
                id=",".join(batch),
                rettype="gb",
                retmode="text"
            )

            for record in SeqIO.parse(handle, "genbank"):
                acc_full = record.id
                acc_base = accession_no_version(record.id)

                locus_tag = None
                old_locus_tag = None
                gene = None
                product = record.description

                for feature in record.features:
                    q = feature.qualifiers

                    if "locus_tag" in q and locus_tag is None:
                        locus_tag = q["locus_tag"][0]

                    if "old_locus_tag" in q and old_locus_tag is None:
                        old_locus_tag = ";".join(q["old_locus_tag"])

                    if "gene" in q and gene is None:
                        gene = q["gene"][0]

                    if "product" in q and product is None:
                        product = q["product"][0]

                seq = str(record.seq)

                mapping[acc_base] = {
                    "protein_accession_ncbi": acc_full,
                    "locus_tag_direct": locus_tag,
                    "old_locus_tag_direct": old_locus_tag,
                    "ncbi_gene": gene,
                    "ncbi_product": product,
                    "protein_sequence": seq,
                }

        except Exception as e:
            print(f"Error fetching batch {batch}: {e}")

        try:
            handle.close()
        except Exception:
            pass

        time.sleep(sleep_sec)

    return mapping


def choose_direct_locus_tag(row):
    """
    Prefer direct PA-style locus tag if available.
    """
    locus = row.get("locus_tag_direct", None)
    old = row.get("old_locus_tag_direct", None)

    if is_valid_text(locus):
        return str(locus).strip()

    if is_valid_text(old):
        # Sometimes old_locus_tag can contain multiple values separated by ;
        vals = [x.strip() for x in str(old).split(";") if x.strip()]
        pa_vals = [x for x in vals if re.match(r"^PA\d+$", x)]

        if pa_vals:
            return pa_vals[0]
        if vals:
            return vals[0]

    return None


def write_unmapped_query_fasta(df, out_fasta):
    """
    Write protein sequences for rows that need BLAST fallback.
    Uses protein_accession_base as FASTA ID.
    """
    sub = df[
        df["locus_tag_direct_final"].isna()
        & df["protein_sequence"].notna()
        & (df["protein_sequence"].astype(str).str.len() > 0)
    ].copy()

    sub = sub.drop_duplicates("protein_accession_base")

    with open(out_fasta, "w") as out:
        for _, row in sub.iterrows():
            acc = row["protein_accession_base"]
            seq = str(row["protein_sequence"])

            if not is_valid_text(acc):
                continue

            out.write(f">{acc}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    print(f"Wrote BLAST query FASTA: {out_fasta}")
    print(f"Sequences for BLAST fallback: {len(sub)}")

    return len(sub)


def make_blast_db_if_needed(pao1_fasta, db_prefix):
    """
    Create BLAST database if needed.
    """
    pin = db_prefix + ".pin"
    phr = db_prefix + ".phr"
    psq = db_prefix + ".psq"

    if os.path.exists(pin) and os.path.exists(phr) and os.path.exists(psq):
        print("BLAST DB already exists.")
        return

    cmd = [
        "makeblastdb",
        "-in", pao1_fasta,
        "-dbtype", "prot",
        "-out", db_prefix
    ]

    print("Running makeblastdb:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def run_blastp(query_fasta, db_prefix, blast_out):
    """
    Run BLASTP against PAO1 protein database.
    """
    if not os.path.exists(query_fasta) or os.path.getsize(query_fasta) == 0:
        print("No query FASTA for BLAST. Skipping BLAST.")
        return False

    cmd = [
        "blastp",
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", blast_out,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-max_target_seqs", "5",
        "-evalue", str(MAX_EVALUE),
        "-num_threads", "4"
    ]

    print("Running BLASTP:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    return True


def get_pao1_locus_from_blast_hit(hit_accession, pao1_gff_map, pao1_header_map):
    """
    Given a BLAST hit protein accession, recover PAO1 locus tag.

    Priority:
    1. GFF protein_id -> locus_tag mapping
    2. FASTA header extraction backup
    """
    if not is_valid_text(hit_accession):
        return None

    hit_full = str(hit_accession)
    hit_base = accession_no_version(hit_full)

    # Best source: GFF
    if hit_full in pao1_gff_map:
        return pao1_gff_map[hit_full]

    if hit_base in pao1_gff_map:
        return pao1_gff_map[hit_base]

    # Backup: FASTA header
    if hit_full in pao1_header_map:
        tag = pao1_header_map[hit_full].get("blast_hit_locus_tag_from_header")
        if is_valid_text(tag):
            return tag

    if hit_base in pao1_header_map:
        tag = pao1_header_map[hit_base].get("blast_hit_locus_tag_from_header")
        if is_valid_text(tag):
            return tag

    return None


def get_pao1_header_from_blast_hit(hit_accession, pao1_header_map):
    if not is_valid_text(hit_accession):
        return None

    hit_full = str(hit_accession)
    hit_base = accession_no_version(hit_full)

    if hit_full in pao1_header_map:
        return pao1_header_map[hit_full].get("blast_hit_header")

    if hit_base in pao1_header_map:
        return pao1_header_map[hit_base].get("blast_hit_header")

    return None


def parse_blast_results(blast_out, pao1_header_map, pao1_gff_map):
    """
    Parse BLAST output and keep best acceptable hit per query.

    qcov = alignment length / query length * 100
    """
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen"
    ]

    empty_df = pd.DataFrame(columns=[
        "protein_accession_base",
        "blast_hit_accession",
        "blast_hit_locus_tag",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "blast_hit_bitscore",
        "blast_hit_header"
    ])

    if not os.path.exists(blast_out) or os.path.getsize(blast_out) == 0:
        return empty_df

    blast = pd.read_csv(blast_out, sep="\t", names=cols)

    if blast.empty:
        return empty_df

    blast["qcov"] = blast["length"] / blast["qlen"] * 100

    blast = blast[
        (blast["pident"] >= MIN_PIDENT)
        & (blast["qcov"] >= MIN_QCOV)
        & (blast["evalue"] <= MAX_EVALUE)
    ].copy()

    if blast.empty:
        return empty_df

    blast = blast.sort_values(
        ["qseqid", "evalue", "bitscore", "pident", "qcov"],
        ascending=[True, True, False, False, False]
    )

    best = blast.drop_duplicates("qseqid").copy()

    best["blast_hit_accession"] = best["sseqid"]
    best["blast_hit_pident"] = best["pident"]
    best["blast_hit_qcov"] = best["qcov"]
    best["blast_hit_evalue"] = best["evalue"]
    best["blast_hit_bitscore"] = best["bitscore"]

    best["blast_hit_locus_tag"] = best["blast_hit_accession"].apply(
        lambda x: get_pao1_locus_from_blast_hit(x, pao1_gff_map, pao1_header_map)
    )

    best["blast_hit_header"] = best["blast_hit_accession"].apply(
        lambda x: get_pao1_header_from_blast_hit(x, pao1_header_map)
    )

    best = best.rename(columns={"qseqid": "protein_accession_base"})

    return best[[
        "protein_accession_base",
        "blast_hit_accession",
        "blast_hit_locus_tag",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "blast_hit_bitscore",
        "blast_hit_header"
    ]]


def final_locus_tag(row):
    """
    Direct mapping first, then BLAST best-hit mapping.
    """
    direct = row.get("locus_tag_direct_final", None)
    blast = row.get("blast_hit_locus_tag", None)

    if is_valid_text(direct):
        return str(direct).strip()

    if is_valid_text(blast):
        return str(blast).strip()

    return None


def mapping_method(row):
    direct = row.get("locus_tag_direct_final", None)
    blast = row.get("blast_hit_locus_tag", None)
    hit = row.get("blast_hit_accession", None)

    if is_valid_text(direct):
        return "direct_ncbi_record"

    if is_valid_text(blast):
        return "blastp_best_hit_to_PAO1"

    if is_valid_text(hit):
        return "blast_hit_no_locus_tag"

    return "unmapped"


# ==========================
# Main
# ==========================
def main():
    if not os.path.exists(PAO1_PROTEIN_FASTA):
        raise FileNotFoundError(
            f"PAO1 protein FASTA not found: {PAO1_PROTEIN_FASTA}\n"
            "Download it first:\n"
            "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/"
            "GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_protein.faa.gz\n"
            "gunzip GCF_000006765.1_ASM676v1_protein.faa.gz"
        )

    if not os.path.exists(PAO1_GFF):
        raise FileNotFoundError(
            f"PAO1 GFF not found: {PAO1_GFF}\n"
            "Download it first:\n"
            "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/"
            "GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz\n"
            "gunzip GCF_000006765.1_ASM676v1_genomic.gff.gz"
        )

    df = pd.read_csv(INCSV)

    df["protein_accession"] = df["Accession"].apply(extract_protein_accession)
    df["protein_accession_base"] = df["protein_accession"].apply(accession_no_version)

    unique_accs = df["protein_accession"].dropna().unique().tolist()

    print(f"Input rows: {len(df)}")
    print(f"Unique protein accessions: {len(unique_accs)}")

    # --------------------------
    # 0. Parse PAO1 references
    # --------------------------
    print("Parsing PAO1 protein FASTA headers...")
    pao1_header_map = parse_pao1_fasta_headers(PAO1_PROTEIN_FASTA)

    print("Parsing PAO1 GFF protein_id -> locus_tag map...")
    pao1_gff_map = parse_pao1_gff(PAO1_GFF)

    # --------------------------
    # 1. Direct NCBI mapping
    # --------------------------
    ncbi_mapping = fetch_ncbi_records(
        unique_accs,
        batch_size=BATCH_SIZE,
        sleep_sec=SLEEP_SEC
    )

    map_df = (
        pd.DataFrame.from_dict(ncbi_mapping, orient="index")
        .reset_index()
        .rename(columns={"index": "protein_accession_base"})
    )

    out = df.merge(map_df, on="protein_accession_base", how="left")

    out["locus_tag_direct_final"] = out.apply(choose_direct_locus_tag, axis=1)

    print(
        f"Direct NCBI mapped: "
        f"{out['locus_tag_direct_final'].notna().sum()} / {len(out)}"
    )

    # --------------------------
    # 2. BLAST fallback for direct-unmapped
    # --------------------------
    n_query = write_unmapped_query_fasta(out, TMP_QUERY_FASTA)

    blast_best_df = pd.DataFrame(columns=[
        "protein_accession_base",
        "blast_hit_accession",
        "blast_hit_locus_tag",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "blast_hit_bitscore",
        "blast_hit_header"
    ])

    if n_query > 0:
        make_blast_db_if_needed(PAO1_PROTEIN_FASTA, BLAST_DB_PREFIX)
        ran_blast = run_blastp(TMP_QUERY_FASTA, BLAST_DB_PREFIX, BLAST_OUT)

        if ran_blast:
            blast_best_df = parse_blast_results(
                BLAST_OUT,
                pao1_header_map=pao1_header_map,
                pao1_gff_map=pao1_gff_map
            )

    out = out.merge(blast_best_df, on="protein_accession_base", how="left")

    # --------------------------
    # 3. Final assignment
    # --------------------------
    out["locus_tag"] = out.apply(final_locus_tag, axis=1)
    out["mapping_method"] = out.apply(mapping_method, axis=1)

    out.to_csv(OUTCSV, index=False)

    print(f"Saved: {OUTCSV}")
    print("")
    print("Mapping summary:")
    print(out["mapping_method"].value_counts(dropna=False))
    print("")
    print(f"Final mapped locus_tag: {out['locus_tag'].notna().sum()} / {len(out)}")

    preview_cols = [
        "Accession",
        "protein_accession",
        "protein_accession_base",
        "gene",
        "locus_tag",
        "mapping_method",
        "locus_tag_direct_final",
        "blast_hit_accession",
        "blast_hit_locus_tag",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "ncbi_gene",
        "ncbi_product",
    ]

    preview_cols = [c for c in preview_cols if c in out.columns]
    print("")
    print("Preview:")
    print(out[preview_cols].head(30))

    # Extra diagnostics
    print("")
    print("Rows with BLAST hit but no PA locus tag:")
    problem = out[
        out["blast_hit_accession"].notna()
        & out["blast_hit_locus_tag"].isna()
        & out["locus_tag_direct_final"].isna()
    ]
    print(len(problem))

    if len(problem) > 0:
        cols = [
            "Accession",
            "protein_accession",
            "gene",
            "blast_hit_accession",
            "blast_hit_pident",
            "blast_hit_qcov",
            "blast_hit_header",
        ]
        cols = [c for c in cols if c in problem.columns]
        print(problem[cols].head(20))


if __name__ == "__main__":
    main()