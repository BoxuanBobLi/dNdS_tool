#!/usr/bin/env python3
"""
Deduplicate a FASTA by header ID (first token) and keep only the longest sequence per ID.

Usage:
  python dedup_longest.py -i input.fasta -o output.dedup.fasta [--map kept_map.tsv]

Notes:
  - The header ID is taken as the first whitespace-delimited token from the FASTA header.
  - If multiple sequences share the same ID, only the longest is kept.
  - Optionally writes a 2-column mapping of kept IDs to selected lengths (and counts).

Requires:
  Biopython  (pip install biopython)
"""

import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def dedup_longest(in_fa: str, out_fa: str, map_out: str | None = None) -> None:
    # Track the best (longest) sequence for each first-token ID
    best = {}            # id -> (length, sequence, id)
    counts = defaultdict(int)

    for rec in SeqIO.parse(in_fa, "fasta"):
        hid = rec.id.split()[0]  # first token in header
        seq = str(rec.seq)
        counts[hid] += 1
        if hid not in best or len(seq) > best[hid][0]:
            best[hid] = (len(seq), seq, hid)

    # Write deduplicated FASTA
    with open(out_fa, "w") as fh:
        for hid, (_, seq, _) in best.items():
            SeqIO.write(SeqRecord(Seq(seq), id=hid, description=""), fh, "fasta")

    # Optional mapping report
    if map_out:
        with open(map_out, "w") as mh:
            mh.write("ID\tKeptLength\tNumDuplicates\n")
            for hid, (L, _, _) in sorted(best.items()):
                mh.write(f"{hid}\t{L}\t{counts[hid]}\n")

def main():
    ap = argparse.ArgumentParser(description="Deduplicate FASTA by header (first token), keep longest sequence.")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA (e.g., *.nogap.fasta)")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA (deduplicated)")
    ap.add_argument("--map", default=None, help="Optional TSV mapping (ID, kept length, duplicate count)")
    args = ap.parse_args()

    dedup_longest(args.input, args.output, args.map)

if __name__ == "__main__":
    main()
