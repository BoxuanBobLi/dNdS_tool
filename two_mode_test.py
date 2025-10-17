import argparse
import os
import itertools
import csv
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
import multiprocessing as mp
import shlex
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.Align import analysis
from io import StringIO

import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonDeprecationWarning)


# NEW: extras for MK
import numpy as np
from collections import Counter
try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

# ---------------- main ----------------
def main():
    parser = argparse.ArgumentParser(description="Calculate dN/dS values (and MK test)")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Pairwise
    pairwise_parser = subparsers.add_parser("pairwise", help="Pairwise dN/dS across all sequences")
    pairwise_parser.add_argument("input_fasta", help="Input codon alignment FASTA file")
    pairwise_parser.add_argument("--format", choices=["long", "matrix"], default="long", help="Output format")
    pairwise_parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    pairwise_parser.add_argument("-t", "--threads", type=int, default=1, help="Threads")
    pairwise_parser.add_argument("--log", default=None, help="Directory to save log file")

    # Groupwise
    groupwise_parser = subparsers.add_parser("groupwise", help="Groupwise dN/dS to a reference sequence")
    groupwise_parser.add_argument("alignment", help="Input codon alignment FASTA file")
    groupwise_parser.add_argument("reference", help="Reference: either single FASTA or alignment")
    groupwise_parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    groupwise_parser.add_argument("-t", "--threads", type=int, default=1, help="Threads")
    groupwise_parser.add_argument("--fast", action="store_true", help="Use most frequent sequence as reference from MSA")
    groupwise_parser.add_argument("--consensus", action="store_true", help="Use consensus sequence as reference from MSA")
    groupwise_parser.add_argument("--log", default=None, help="Directory to save log file")

    # MK test (aligned with groupwise reference handling)
    mk_parser = subparsers.add_parser("mk", help="McDonald–Kreitman test (ingroup alignment vs reference)")
    mk_parser.add_argument("ingroup_alignment", help="Ingroup codon alignment FASTA (>=2 sequences)")
    mk_parser.add_argument("reference", help="Reference: either single FASTA or alignment (same behavior as groupwise)")
    mk_parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    mk_parser.add_argument("--outname", default=None, help="Output filename (TSV); default: <ingroup>_MK.tsv")
    mk_parser.add_argument("-t", "--threads", type=int, default=1, help="Threads for reference derivation")
    mk_parser.add_argument("--fast", action="store_true", help="Use most frequent sequence as reference from MSA")
    mk_parser.add_argument("--consensus", action="store_true", help="Use consensus sequence as reference from MSA")
    mk_parser.add_argument("--log", default=None, help="Directory to save log file")

    args = parser.parse_args()

    if args.mode == "pairwise":
        records = sanitize_records(args.input_fasta)
        output_base = os.path.splitext(os.path.basename(args.input_fasta))[0].split('.')[0]
        output_path = os.path.join(args.output_dir, f"{output_base}_pairwise_dnds.csv")
        log_path = _prepare_log(args.log, args.output_dir, output_base)
        write_to_log("Starting dN/dS calculation with cmd: " + ' '.join(os.sys.argv), log_path)
        dnds_matrix = compute_dnds_matrix(records, threads=args.threads, log_path=log_path)
        if args.format == "long":
            write_long_format(output_path, dnds_matrix)
        else:
            write_matrix_format(output_path, dnds_matrix, [r.id for r in records])
        print(f"Pairwise dN/dS written to {output_path}")

    elif args.mode == "groupwise":
        records = sanitize_records(args.alignment)
        # ref could be single or alignment
        if len(list(SeqIO.parse(args.reference, "fasta"))) == 1:
            ref_record = next(SeqIO.parse(args.reference, "fasta"))
        else:
            if args.fast and args.consensus:
                raise ValueError("Cannot use both --fast and --consensus.")
            elif args.consensus:
                ref_record = find_reference_consensus(args.reference, args.output_dir)
            elif args.fast:
                ref_record = find_reference_fast(args.reference, args.output_dir, threads=args.threads)
            else:
                ref_record = find_reference_from_msa(args.reference, args.output_dir, threads=args.threads)

        output_base = os.path.splitext(os.path.basename(args.alignment))[0].split('.')[0]
        output_path = os.path.join(args.output_dir, f"{output_base}_groupwise_dnds.csv")
        log_path = _prepare_log(args.log, args.output_dir, output_base)
        write_to_log("Starting dN/dS calculation with cmd: " + ' '.join(os.sys.argv), log_path)
        dnds_matrix = compute_groupwise_dnds(records, ref_record, threads=args.threads, log_path=log_path)
        write_long_format_ref(output_path, dnds_matrix)
        print(f"Groupwise dN/dS written to {output_path}")


    elif args.mode == "mk":
        # Load ingroup alignment (sanitized)
        ingroup = sanitize_records(args.ingroup_alignment)
        if len(ingroup) < 2:
            raise ValueError("MK test requires ≥2 ingroup sequences.")

        # Determine reference (same logic as groupwise)
        ref_records = list(SeqIO.parse(args.reference, "fasta"))
        if len(ref_records) == 1:
            ref_record = ref_records[0]
        else:
            if args.fast and args.consensus:
                raise ValueError("Cannot use both --fast and --consensus for MK.")
            elif args.consensus:
                ref_record = find_reference_consensus(args.reference, args.output_dir)
            elif args.fast:
                ref_record = find_reference_fast(args.reference, args.output_dir, threads=args.threads)
            else:
                ref_record = find_reference_from_msa(args.reference, args.output_dir, threads=args.threads)

        # Compute MK
        mk = compute_mk(ingroup, ref_record)

        # Output
        base = os.path.splitext(os.path.basename(args.ingroup_alignment))[0]
        outname = args.outname or f"{base}_MK.tsv"
        os.makedirs(args.output_dir, exist_ok=True)
        out_path = os.path.join(args.output_dir, outname)
        write_mk_tsv(out_path, mk, ingroup_name=base, out_name=ref_record.id)
        print(f"MK test written to {out_path}")

# -------------- helpers (shared) --------------
def _prepare_log(log_dir, out_dir, base):
    if log_dir:
        if os.path.exists(os.path.join(log_dir, f"{base}_errors.log")):
            os.remove(os.path.join(log_dir, f"{base}_errors.log"))
        os.makedirs(log_dir, exist_ok=True)
        return os.path.join(log_dir, f"{base}_errors.log")
    else:
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, f"{base}_errors.log")

def sanitize_records(fasta_path):
    sanitized = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id.startswith("_R_"):
            record.id = record.id.replace("_R_", "")
            record.description = record.description.replace("_R_", "")
            record.seq = record.seq.reverse_complement()
        sanitized.append(record)
    return sanitized

def sanitize_list(records):
    out = []
    for r in records:
        rr = SeqRecord(Seq(str(r.seq)), id=r.id.replace("_R_", ""), description=r.description.replace("_R_", ""))
        out.append(rr)
    return out

# (gap/quality funcs unchanged; omitted for brevity)
# ... your remove_invalid_codons / compute_dnds_pair / compute_dnds_matrix / etc. remain unchanged ...

# ---------------- MK test utilities (NEW) ----------------
_valid_nt = set("ACGT")

def _is_valid_codon(c):
    return len(c) == 3 and all(n in _valid_nt for n in c)

def _one_nt_diff(c1, c2):
    return sum(a != b for a, b in zip(c1, c2)) == 1

def _aa(codon):
    try:
        return str(Seq(codon).translate(table=11))  # bacterial table
    except Exception:
        return None

def _is_synonymous(c1, c2):
    if not (_is_valid_codon(c1) and _is_valid_codon(c2)):
        return None
    aa1, aa2 = _aa(c1), _aa(c2)
    if aa1 is None or aa2 is None or "*" in (aa1, aa2):
        return None
    return aa1 == aa2

def _consensus_codon(codons):
    # majority codon among valid triplets; fall back to most common including invalid
    valid = [c for c in codons if _is_valid_codon(c)]
    pool = valid if valid else codons
    if not pool:
        return None
    return Counter(pool).most_common(1)[0][0]

def _slice_codons(seq):
    s = str(seq).upper().replace("U", "T")
    return [s[i:i+3] for i in range(0, len(s) - (len(s)%3), 3)]

def compute_mk_with_correction(Pn, Ps, Dn, Ds, pseudo=0.5, alternative='two-sided'):
    """
    Compute MK test statistics with continuity correction to avoid division by zero.
    """
    # Fisher Exact Test (on raw counts)
    table = np.array([[Pn, Ps], [Dn, Ds]])
    _, fisher_p = fisher_exact(table, alternative=alternative)

    # Continuity correction on counts for NI/alpha
    Pn_c = Pn + pseudo
    Ps_c = Ps + pseudo
    Dn_c = Dn + pseudo
    Ds_c = Ds + pseudo

    NI = (Pn_c / Ps_c) / (Dn_c / Ds_c)
    alpha = 1.0 - NI

    return NI, alpha, fisher_p

def compute_mk(ingroup_records, outgroup_record):
    # Build per-codon arrays
    ing_codons = [ _slice_codons(r.seq) for r in ingroup_records ]
    out_codons = _slice_codons(outgroup_record.seq)
    L = min(min(len(c) for c in ing_codons), len(out_codons))

    Pn = Ps = Dn = Ds = 0

    for i in range(L):
        ing_site = [c[i] for c in ing_codons]
        out_c = out_codons[i]

        if not _is_valid_codon(out_c):
            continue

        cons = _consensus_codon(ing_site)
        if cons is None or not _is_valid_codon(cons):
            continue

        # --- Polymorphism classification (site-based, conservative) ---
        alleles = set(ing_site)
        if len(alleles) > 1:
            has_nonsyn = False
            only_syn = True
            for alt in alleles:
                if alt == cons:
                    continue
                if not _is_valid_codon(alt):
                    continue
                if not _one_nt_diff(cons, alt):
                    # skip multi-step changes conservatively
                    continue
                syn = _is_synonymous(cons, alt)
                if syn is None:
                    continue
                if syn:
                    pass
                else:
                    has_nonsyn = True
                    only_syn = False
            if has_nonsyn:
                Pn += 1
            elif only_syn:
                Ps += 1

        # --- Divergence classification ---
        if cons != out_c:
            if not _one_nt_diff(cons, out_c):
                continue
            syn_div = _is_synonymous(cons, out_c)
            if syn_div is None:
                continue
            if syn_div:
                Ds += 1
            else:
                Dn += 1

    # Compute Fisher's exact test on raw counts
    table = np.array([[Pn, Ps], [Dn, Ds]])
    try:
        _, fisher_p = fisher_exact(table, alternative='two-sided')
    except Exception:
        fisher_p = np.nan

    # Compute NI and alpha safely
    if Ps == 0 or Ds == 0 or Dn == 0:
        # Use continuity correction
        NI, alpha, _ = compute_mk_with_correction(Pn, Ps, Dn, Ds, pseudo=0.5, alternative='two-sided')
    else:
        NI = (Pn / Ps) / (Dn / Ds)
        alpha = 1.0 - NI

    return {
        "Pn": int(Pn), "Ps": int(Ps), "Dn": int(Dn), "Ds": int(Ds),
        "NI": float(NI), "alpha": float(alpha), "p_value": float(fisher_p)
    }


def write_mk_tsv(path, mk_dict, ingroup_name, out_name):
    with open(path, "w") as f:
        f.write("Ingroup\tOutgroup\tPn\tPs\tDn\tDs\tNI\talpha\tFisher_p\n")
        f.write(f"{ingroup_name}\t{out_name}\t{mk_dict['Pn']}\t{mk_dict['Ps']}\t"
                f"{mk_dict['Dn']}\t{mk_dict['Ds']}\t"
                f"{_fmt(mk_dict['NI'])}\t{_fmt(mk_dict['alpha'])}\t{_fmt(mk_dict['p_value'])}\n")

def _fmt(x):
    try:
        return f"{x:.4g}"
    except Exception:
        return "NA"

# ---------------- your existing functions (unchanged) ----------------
# remove_invalid_codons, compute_dnds_pair, compute_dnds_matrix, etc.
# (keep your implementations exactly as you had them)
# --------------------------------------------------------------------

# Reference derivation functions (unchanged)
def find_reference_from_msa(msa_path, work_dir=".", node_name=None, save_path=None, threads=1):
    cmd = f"cd {shlex.quote(work_dir)} && iqtree2 -s {shlex.quote(msa_path)} --ancestral -redo -T {threads}"
    subprocess.run(cmd, shell=True, check=True)
    base = os.path.basename(msa_path)
    anc_file = os.path.join(work_dir, base + ".state")
    tree_file = os.path.join(work_dir, base + ".treefile")
    if node_name is None:
        with open(tree_file) as f:
            tree_str = f.read().strip()
        last_paren = tree_str.rfind(")")
        root_str = tree_str[last_paren+1:].strip(" :;,")
        node_name = root_str if root_str.startswith("Node") else None
        if node_name is None:
            raise RuntimeError("Could not determine root node from tree file.")
    sequence, header = [], None
    with open(anc_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if header is None:
                header = parts
                idx_node = header.index("Node")
                idx_state = header.index("State")
                continue
            if parts[idx_node] == node_name:
                sequence.append(parts[idx_state])
    if not sequence:
        raise RuntimeError(f"No sequence reconstructed for node '{node_name}' in file {anc_file}")
    record = SeqRecord(Seq("".join(sequence)), id=node_name, description="Ancestral sequence from IQ-TREE")
    if save_path is None:
        save_path = os.path.join(work_dir, base + ".ancestral.fasta")
    SeqIO.write(record, save_path, "fasta")
    return record

def find_reference_fast(msa_path, work_dir=".", node_name=None, save_path=None, threads=1):
    records = list(SeqIO.parse(msa_path, "fasta"))
    seq_counts = Counter(str(r.seq).upper() for r in records)
    most_common_seq, _ = seq_counts.most_common(1)[0]
    ref_record = next(r for r in records if str(r.seq).upper() == most_common_seq)
    ref_record = SeqRecord(Seq(str(ref_record.seq)), id="ref", description="Most common sequence in MSA")
    if save_path is None:
        save_path = os.path.join(work_dir, os.path.basename(msa_path).split(".")[0] + "_fast_ref.fasta")
    SeqIO.write(ref_record, save_path, "fasta")
    return ref_record

def find_reference_consensus(msa_path, work_dir=".", save_path=None):
    records = list(SeqIO.parse(msa_path, "fasta"))
    alignment = MultipleSeqAlignment(records)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(threshold=0.7, ambiguous='N')
    ref_record = SeqRecord(consensus, id="ref", name="ref", description="Consensus sequence from MSA")
    if save_path is None:
        save_path = os.path.join(work_dir, os.path.basename(msa_path).split(".")[0] + "_consensus_ref.fasta")
    SeqIO.write(ref_record, save_path, "fasta")
    return ref_record

# dN/dS writers unchanged
def write_long_format(output_path, dnds_matrix):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Seq1", "Seq2", "dN", "dS", "dN/dS"])
        for (id1, id2), (dN, dS) in dnds_matrix.items():
            ratio = None
            if dN is not None and dS is not None:
                if dS != 0:
                    ratio = dN / dS
                elif dS == 0 and dN != 0:
                    ratio = float('inf')
                elif dS == 0 and dN == 0:
                    ratio = 0.0
            writer.writerow([
                id1, id2,
                f"{dN:.4f}" if dN is not None else "NA",
                f"{dS:.4f}" if dS is not None else "NA",
                f"{ratio:.4f}" if ratio is not None else "NA"
            ])

def write_long_format_ref(output_path, dnds_matrix):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Ref", "Seq", "dN", "dS", "dN/dS"])
        for (id1, id2), (dN, dS) in dnds_matrix.items():
            ratio = None
            if dN is not None and dS is not None:
                if dS != 0:
                    ratio = dN / dS
                elif dS == 0 and dN != 0:
                    ratio = float('inf')
                elif dS == 0 and dN == 0:
                    ratio = 0.0
            writer.writerow([
                id1, id2,
                f"{dN:.4f}" if dN is not None else "NA",
                f"{dS:.4f}" if dS is not None else "NA",
                f"{ratio:.4f}" if ratio is not None else "NA"
            ])

def write_matrix_format(output_path, dnds_matrix, ids):
    ids_sorted = sorted(ids)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = [""] + ids_sorted
        writer.writerow(header)
        for id1 in ids_sorted:
            row = [id1]
            for id2 in ids_sorted:
                if id1 == id2:
                    row.append("0.0000")
                else:
                    key = (id1, id2) if (id1, id2) in dnds_matrix else (id2, id1)
                    dN, dS = dnds_matrix.get(key, (None, None))
                    ratio = None
                    if dN is not None and dS is not None:
                        if dS != 0:
                            ratio = dN / dS
                        elif dS == 0 and dN != 0:
                            ratio = float("inf")
                        elif dS == 0 and dN == 0:
                            ratio = 0.0
                    row.append(f"{ratio:.4f}" if ratio is not None else "NA")
            writer.writerow(row)

def write_to_log(message, log_path):
    with open(log_path, "a") as log_file:
        log_file.write(message + "\n")

if __name__ == "__main__":
    main()
