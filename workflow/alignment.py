import pandas as pd
import os
import argparse
from collections import defaultdict
from Bio.Seq import Seq

def process_csv_to_fasta(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    try:
        df = pd.read_csv(input_file, header=None)
    except Exception as e:
        print(f"Failed to read {input_file}: {e}")
        return

    groups = defaultdict(list)
    for _, row in df.iterrows():
        accession = row[0].split('(')[0]
        genome_id = row[8]
        sequence = row[6]
        direction = row[7]
        if pd.notna(sequence) and pd.notna(genome_id):
            groups[accession].append((genome_id, sequence, direction))

    for accession, entries in groups.items():
        fasta_file = os.path.join(output_dir, f"{accession}.fasta")
        with open(fasta_file, "w") as f:
            for genome_id, seq, direction in entries:
                if direction == "-":
                    seq = str(Seq(seq).reverse_complement())
                f.write(f">{genome_id}\n{seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate FASTA from a single VFDB CSV with consistent direction")
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to a single VF CSV file")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Directory to write FASTA files")
    args = parser.parse_args()

    process_csv_to_fasta(args.input_file, args.output_dir)




# nice -19 mafft --auto /data1/B_Li/vfdb/all_samples/reformat_output_qseq_direction/All_in_one/VFG000122.fasta > /data1/B_Li/vfdb/all_samples/reformat_output_qseq_direction/All_in_one/VFG000122_aligned.fasta
# nice -19 mafft --auto --thread 100 /data1/B_Li/vfdb/all_samples/reformat_output_qseq_direction/All_in_one/VFG000122.fasta > /data1/B_Li/vfdb/all_samples/reformat_output_qseq_direction/All_in_one/VFG000122_aligned_threads.fasta