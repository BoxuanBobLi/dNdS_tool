import os
import pandas as pd
from Bio import SeqIO
import argparse

def split_fasta_by_host(metadata_csv, input_dir, output_dir):
    # Load metadata mapping
    trait_df = pd.read_csv(metadata_csv)
    trait_dict = dict(zip(trait_df['Genome_ID'], trait_df['Trait'].astype(str)))  # Ensure it's str: '1' or '2'

    host_dirs = {
        '0': os.path.join(output_dir, '0_splitted_aln'),      # e.g. animal/human
        '1': os.path.join(output_dir, '1_splitted_aln')   # e.g. environment
    }
    for d in host_dirs.values():
        os.makedirs(d, exist_ok=True)

    for root, _, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".nt_ali.fasta"):
                input_path = os.path.join(root, filename)
                host_records = {'0': [], '1': []}

                for record in SeqIO.parse(input_path, "fasta"):
                    genome_id = record.id
                    host_label = trait_dict.get(genome_id)
                    if host_label in host_records:
                        host_records[host_label].append(record)

                for host_label, records in host_records.items():
                    if records:
                        out_path = os.path.join(host_dirs[host_label], filename)
                        SeqIO.write(records, out_path, "fasta")
                        print(f"Wrote {len(records)} sequences to {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTA alignment files by host-relatedness")
    parser.add_argument("-m", "--metadata_csv", required=True, help="Path to metadata_with_host.csv with Genome_ID and Host_Related columns")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing alignment FASTA files")
    parser.add_argument("-o", "--output_dir", required=True, help="Base output directory for 0 and 1")
    args = parser.parse_args()

    split_fasta_by_host(args.metadata_csv, args.input_dir, args.output_dir)
