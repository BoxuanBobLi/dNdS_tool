import argparse
import os
import itertools
import csv
import subprocess
from Bio import SeqIO
# from Bio.codonalign.codonseq import CodonSeq
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

def main():
    parser = argparse.ArgumentParser(description="Calculate dN/dS values")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Pairwise mode
    pairwise_parser = subparsers.add_parser("pairwise", help="Pairwise dN/dS across all sequences")
    pairwise_parser.add_argument("input_fasta", help="Input codon alignment FASTA file")
    pairwise_parser.add_argument("--format", choices=["long", "matrix"], default="long", help="Output format")
    pairwise_parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    pairwise_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for parallel processing")
    pairwise_parser.add_argument("--log", default=None, help="Directory to save log file")

    # Groupwise mode
    groupwise_parser = subparsers.add_parser("groupwise", help="Groupwise dN/dS to a reference sequence")
    groupwise_parser.add_argument("alignment", help="Input codon alignment FASTA file")
    groupwise_parser.add_argument("reference", help="Reference: either single FASTA or alignment")
    groupwise_parser.add_argument("-o", "--output_dir", default=".", help="Output directory")
    groupwise_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for parallel processing")
    groupwise_parser.add_argument("--fast", action="store_true", help="Use most frequent sequence as reference from MSA")
    groupwise_parser.add_argument("--consensus", action="store_true", help="Use consensus sequence as reference from MSA")
    groupwise_parser.add_argument("--log", default=None, help="Directory to save log file")

    args = parser.parse_args()


    if args.mode == "pairwise":
        #records = list(SeqIO.parse(args.input_fasta, "fasta"))
        records = sanitize_records(args.input_fasta)

        # check_alignment_quality(records, 
        #     sequence_gap_threshold=0.3,
        #     bad_sequence_fraction_threshold=0.3,
        #     column_gap_threshold=0.5)

        output_base = os.path.splitext(os.path.basename(args.input_fasta))[0].split('.')[0]
        output_path = os.path.join(args.output_dir, f"{output_base}_pairwise_dnds.csv")

        if args.log:
            if os.path.exists(os.path.join(args.log, f"{output_base}_errors.log")):
                os.remove(os.path.join(args.log, f"{output_base}_errors.log"))
            os.makedirs(args.log, exist_ok=True)
            log_path = os.path.join(args.log, f"{output_base}_errors.log")
        else:
            os.makedirs(args.output_dir, exist_ok=True)
            log_path = os.path.join(args.output_dir, f"{output_base}_errors.log")
        write_to_log("Starting dN/dS calculation with cmd: " + ' '.join(os.sys.argv), log_path)

        dnds_matrix = compute_dnds_matrix(records, threads=args.threads, log_path=log_path)
        if args.format == "long":
            write_long_format(output_path, dnds_matrix)
        else:
            write_matrix_format(output_path, dnds_matrix, [r.id for r in records])
        print(f"Pairwise dN/dS written to {output_path}")

    elif args.mode == "groupwise":
        records = sanitize_records(args.alignment)
        #records = list(SeqIO.parse(args.alignment, "fasta"))
        # count_of_sequences_with_gaps = sum(1 for record in records if record.seq.count("-") / len(record.seq) > 0.30)
        # if count_of_sequences_with_gaps / len(records) > 0.30:
        #     raise ValueError("Input sequences contain more than 30% gaps. Please check your input alignment.")

        if len(list(SeqIO.parse(args.reference, "fasta"))) == 1:
            ref_record = next(SeqIO.parse(args.reference, "fasta"))
        else:
            if args.fast and args.consensus:
                raise ValueError("Cannot use both --fast and --consensus. Choose one method to derive reference from MSA.")
            elif args.consensus:
                ref_record = find_reference_consensus(args.reference, args.output_dir)
            elif args.fast:
                ref_record = find_reference_fast(args.reference, args.output_dir, threads=args.threads)
            else:
                ref_record = find_reference_from_msa(args.reference, args.output_dir, threads=args.threads)

        output_base = os.path.splitext(os.path.basename(args.alignment))[0].split('.')[0]
        output_path = os.path.join(args.output_dir, f"{output_base}_groupwise_dnds.csv")
        
        if args.log:
            if os.path.exists(os.path.join(args.log, f"{output_base}_errors.log")):
                os.remove(os.path.join(args.log, f"{output_base}_errors.log"))
            os.makedirs(args.log, exist_ok=True)
            log_path = os.path.join(args.log, f"{output_base}_errors.log")
        else:
            os.makedirs(args.output_dir, exist_ok=True)
            log_path = os.path.join(args.output_dir, f"{output_base}_errors.log")
        write_to_log("Starting dN/dS calculation with cmd: " + ' '.join(os.sys.argv), log_path)

        dnds_matrix = compute_groupwise_dnds(records, ref_record, threads=args.threads, log_path=log_path)
        write_long_format_ref(output_path, dnds_matrix)
        print(f"Groupwise dN/dS written to {output_path}")


def sanitize_records(fasta_path):
    """
    Parses and sanitizes records from a FASTA file.
    - Removes '_R_' from IDs/descriptions.
    - Reverse complements sequences if they start with '_R_'.

    Returns a list of sanitized SeqRecord objects.
    """
    sanitized = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id.startswith("_R_"):
            record.id = record.id.replace("_R_", "")
            record.description = record.description.replace("_R_", "")
            record.seq = record.seq.reverse_complement()
        sanitized.append(record)
    return sanitized

def check_alignment_quality(records, sequence_gap_threshold=0.3, bad_sequence_fraction_threshold=0.3, column_gap_threshold=0.7):
    """
    Check alignment quality by identifying sequences with high gaps,
    ignoring alignment columns that are mostly gaps.
    """
    seqs = [str(r.seq) for r in records]
    num_seqs = len(seqs)
    seq_length = len(seqs[0])

    # Step 1: identify bad columns (too gappy)
    column_gaps = [sum(1 for s in seqs if s[i] in ['-', '.']) for i in range(seq_length)]
    bad_columns = {i for i, gap_count in enumerate(column_gaps) if gap_count / num_seqs > column_gap_threshold}

    # Step 2: per-sequence gap fraction (excluding bad columns)
    def effective_gap_fraction(seq):
        total, gaps = 0, 0
        for i, base in enumerate(seq):
            if i in bad_columns:
                continue
            total += 1
            if base in ['-', '.']:
                gaps += 1
        return gaps / total if total > 0 else 0

    num_bad_seqs = sum(1 for s in seqs if effective_gap_fraction(s) > sequence_gap_threshold)

    if num_bad_seqs / num_seqs > bad_sequence_fraction_threshold:
        raise ValueError(
            f"Alignment rejected: {num_bad_seqs}/{num_seqs} ({100*num_bad_seqs/num_seqs:.1f}%) "
            f"sequences contain >{sequence_gap_threshold*100:.0f}% gaps "
            f"(excluding columns where >{column_gap_threshold*100:.0f}% of sequences are gaps)."
        )


def clean_stop_codons(seq1, seq2):
    # Uppercase
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    # Convert to lists to allow mutation
    seq1_list = list(seq1)
    seq2_list = list(seq2)
    
    stop_codons = {"TAA", "TAG", "TGA"}
    
    # Loop through codons in-frame
    for i in range(0, len(seq1) - 2, 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]
        
        if codon1 in stop_codons or codon2 in stop_codons:
            # Replace codon with gaps in both sequences
            seq1_list[i:i+3] = ['-'] * 3
            seq2_list[i:i+3] = ['-'] * 3

    return ''.join(seq1_list), ''.join(seq2_list)


def resolve_ambiguity(seq, ref_seq):
    resolved = []
    for base, ref_base in zip(seq, ref_seq):
        if base.upper() in "ACGT-":
            resolved.append(base)
        elif ref_base.upper() in "ACGT-":
            resolved.append(ref_base)
        else:
            resolved.append("A")
    return ''.join(resolved)


def remove_shared_gaps_frame_preserving(seq1, seq2):
    if len(seq1) != len(seq2):
        print(len(seq1), len(seq2))
        raise ValueError("Sequences must be of the same length")

    kept1 = []
    kept2 = []

    i = 0
    while i < len(seq1):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]

        if len(codon1) < 3 or len(codon2) < 3:
            # Skip partial codons at the end (or raise error if needed)
            break

        # Check if all 3 positions are gaps in both sequences
        if codon1 == "---" and codon2 == "---":
            i += 3
            continue  # skip shared gap codon

        kept1.append(codon1)
        kept2.append(codon2)
        i += 3

    cleaned_seq1 = ''.join(kept1)
    cleaned_seq2 = ''.join(kept2)


    # Final check for codon frame
    if len(cleaned_seq1) % 3 != 0 or len(cleaned_seq2) % 3 != 0:
        raise ValueError("Resulting sequences are not in codon frame after gap removal")

    return cleaned_seq1, cleaned_seq2


def remove_invalid_codons(seq1, seq2):
    """
    Remove any codons that contain gaps or unequal length across seq1 and seq2.
    Ensures resulting seqs are equal-length and codon-aligned.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # remove shared gaps
    seq1, seq2 = remove_shared_gaps_frame_preserving(seq1, seq2)

    # Resolve ambiguity
    seq1 = resolve_ambiguity(str(seq1), str(seq2))
    seq2 = resolve_ambiguity(str(seq2), str(seq1))

    # Remove stop codons
    seq1, seq2 = clean_stop_codons(seq1, seq2)


    new_seq1 = []
    new_seq2 = []

    for i in range(0, min(len(seq1), len(seq2)) - 2, 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]
        if "-" in codon1 or "-" in codon2:
            continue  # skip codons with gaps in either sequence
        if len(codon1) == 3 and len(codon2) == 3:
            new_seq1.append(codon1)
            new_seq2.append(codon2)

    return ''.join(new_seq1), ''.join(new_seq2)

def compute_dnds_pair(args):
    rec1, rec2, log_path = args
    if rec1.id == rec2.id:
        return None
    # error handling for invalid sequences
    print(f"Processing pair: {rec1.id} vs {rec2.id}, lengths: {len(rec1.seq)}, {len(rec2.seq)}")
    seq1, seq2 = remove_invalid_codons(str(rec1.seq), str(rec2.seq))

    data = f">seq1\n{seq1}\n>seq2\n{seq2}\n"
    stream = StringIO(data)
    
    try:
        alignment = Align.read(stream, "fasta")
        alignment = alignment[:, :-3]  # Optionally trim codon tail
        dN, dS = analysis.calculate_dn_ds(alignment, method="NG86")
    except Exception as e:
        message = (
            f"Failed alignment for {rec1.id} vs {rec2.id}: {e}\n"
            f"seq1 length: {len(seq1)}, seq2 length: {len(seq2)}\n"
            f"seq1: {seq1}\n"
            f"seq2: {seq2}\n"
            f"seq1 codons: {[seq1[i:i+3] for i in range(0, len(seq1)-2, 3)]}\n"
            f"seq2 codons: {[seq2[i:i+3] for i in range(0, len(seq2)-2, 3)]}\n"
            + "-" * 80
        )
        write_to_log(message, log_path)
        dN, dS = None, None

    finally:
        stream.close()

    return ((rec1.id, rec2.id), (dN, dS))

def compute_dnds_matrix(records, threads=1, log_path=None):
    pairs = [(rec1, rec2, log_path) for rec1, rec2 in itertools.combinations(records, 2)]
    with mp.Pool(threads) as pool:
        results = list(tqdm(pool.imap(compute_dnds_pair, pairs), total=len(pairs), desc="Computing pairwise dN/dS"))
    return dict(results)



def compute_group_dnds_worker(args):
    rec, ref_record, log_path = args 

        
    if rec.id == ref_record.id:
        return None

    seq1, seq2 = remove_invalid_codons(str(rec.seq), str(ref_record.seq))

    data = f">seq1\n{seq1}\n>seq2\n{seq2}\n"
    stream = StringIO(data)

    try:
        alignment = Align.read(stream, "fasta")
        alignment = alignment[:, :-3]  # Optionally trim codon tail
        dN, dS = analysis.calculate_dn_ds(alignment, method="NG86")
    except Exception as e:
        message = (
            f"Failed alignment for {rec.id} vs {ref_record.id}: {e}\n"
            f"seq1 length: {len(seq1)}, seq2 length: {len(seq2)}\n"
            f"seq1: {seq1}\n"
            f"seq2: {seq2}\n"
            f"seq1 codons: {[seq1[i:i+3] for i in range(0, len(seq1)-2, 3)]}\n"
            f"seq2 codons: {[seq2[i:i+3] for i in range(0, len(seq2)-2, 3)]}\n"
            + "-" * 80
        )
        write_to_log(message, log_path)
        dN, dS = None, None
    finally:
        stream.close()

    return ((ref_record.id, rec.id), (dN, dS))


def compute_groupwise_dnds(records, ref_record, threads=1, log_path=None):
    args = [(rec, ref_record, log_path) for rec in records if rec.id != ref_record.id]
    with mp.Pool(threads) as pool:
        results = list(tqdm(pool.imap(compute_group_dnds_worker, args), total=len(args), desc="Computing groupwise dN/dS"))
    return {k: v for k, v in results if k is not None}


def find_reference_from_msa(msa_path, work_dir=".", node_name=None, save_path=None, threads=1):
    # Run IQ-TREE ancestral state reconstruction
    cmd = f"cd {shlex.quote(work_dir)} && iqtree2 -s {shlex.quote(msa_path)} --ancestral -redo -T {threads}"
    subprocess.run(cmd, shell=True, check=True)


    base = os.path.basename(msa_path)
    anc_file = os.path.join(work_dir, base + ".state")
    tree_file = os.path.join(work_dir, base + ".treefile")

    # Automatically find root node if not provided
    if node_name is None:
        with open(tree_file) as f:
            tree_str = f.read().strip()
        last_paren = tree_str.rfind(")")
        root_str = tree_str[last_paren+1:].strip(" :;,")
        node_name = root_str if root_str.startswith("Node") else None
        if node_name is None:
            raise RuntimeError("Could not determine root node from tree file.")
    print(f"Node name: {node_name}")
    # Reconstruct sequence
    sequence = []
    header = None
    with open(anc_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
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

    # Save automatically if no path specified
    if save_path is None:
        save_path = os.path.join(work_dir, base + ".ancestral.fasta")
    SeqIO.write(record, save_path, "fasta")
    print(f"Saved to: {save_path}")

    return record

# most frequent sequence
def find_reference_fast(msa_path, work_dir=".", node_name=None, save_path=None, threads=1):
    # find the most common sequence in the MSA
    records = list(SeqIO.parse(msa_path, "fasta"))
    seq_counts = {}
    for record in records:
        seq_str = str(record.seq).upper()
        if seq_str not in seq_counts:
            seq_counts[seq_str] = 0
        seq_counts[seq_str] += 1
    most_common_seq = max(seq_counts, key=seq_counts.get)
    most_common_count = seq_counts[most_common_seq]
    print(f"Most common sequence count: {most_common_count}")
    most_common_records = [record for record in records if str(record.seq).upper() == most_common_seq]
    if len(most_common_records) > 1:
        # If there are multiple records with the same sequence, choose the first one
        ref_record = most_common_records[0]
    else:
        ref_record = most_common_records[0]
    ref_record.id = "ref"
    ref_record.description = "Most common sequence in MSA"
    if save_path is None:
        save_path = os.path.join(work_dir, os.path.basename(msa_path).split(".")[0] + "_fast_ref.fasta")
    SeqIO.write(ref_record, save_path, "fasta")
    print(f"Saved reference sequence to: {save_path}")
    #print how long the sequence is
    print(f"Reference sequence length: {len(ref_record.seq)}")
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
    print(f"Saved consensus reference sequence to: {save_path}")
    #pring how long the sequence is
    print(f"Consensus sequence length: {len(consensus)}")
    return ref_record


def write_long_format(output_path, dnds_matrix):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Seq1", "Seq2", "dN", "dS", "dN/dS"])
        for (id1, id2), (dN, dS) in dnds_matrix.items():
            if dN is not None and dS is not None:
                if dS != 0 and dN != 0 or dS != 0 and dN == 0:
                    ratio = dN / dS
                elif dS == 0 and dN != 0:
                    ratio = float('inf')
                elif dS == 0 and dN == 0:
                    ratio = 0.0
            else:
                ratio = None
            writer.writerow([
                id1,
                id2,
                f"{dN:.4f}" if dN is not None else "NA",
                f"{dS:.4f}" if dS is not None else "NA",
                f"{ratio:.4f}" if ratio is not None else "NA"
            ])

# write long with sef generated reference (only difference in first column names)
def write_long_format_ref(output_path, dnds_matrix):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Ref", "Seq", "dN", "dS", "dN/dS"])
        for (id1, id2), (dN, dS) in dnds_matrix.items():
            if dN is not None and dS is not None:
                if dS != 0 and dN != 0 or dS != 0 and dN == 0:
                    ratio = dN / dS
                elif dS == 0 and dN != 0:
                    ratio = float('inf')
                elif dS == 0 and dN == 0:
                    ratio = 0.0
            else:
                ratio = None
            writer.writerow([
                id1,
                id2,
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
                    if dN is not None and dS is not None:
                        if dS != 0 and dN != 0 or dS != 0 and dN == 0:
                            ratio = dN / dS
                        elif dS == 0 and dN != 0:
                            ratio = float("inf")
                        elif dS == 0 and dN == 0:
                            ratio = 0.0
                        else:
                            ratio = None
                    else:
                        ratio = None
                    row.append(f"{ratio:.4f}" if ratio is not None else "NA")
            writer.writerow(row)

def write_to_log(message, log_path):
    with open(log_path, "a") as log_file:
        log_file.write(message + "\n")


if __name__ == "__main__":
    main()

### Example Usage:
# pairwise: python script.py pairwise alignment.fasta --format matrix -o results/
# groupwise: python script.py groupwise alignment.fasta ref.fasta -o results/
# groupwise & iqtree refered reference: python script.py groupwise alignment.fasta ref_alignment.fasta -o results/

# Note: remove intermediate files generated by iqtree
# Note: add a progress bar
# Add a error log output

# dealing with dN and dS is all 0
# dealing with dN is nonzero and dS is zero


## count all syn nonsyn sites and calculate overall dN/dS
## filter out 20% gap sequences