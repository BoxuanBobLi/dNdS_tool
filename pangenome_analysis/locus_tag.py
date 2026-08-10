#!/usr/bin/env python3

import ssl

# Work around the server's self-signed certificate chain
ssl._create_default_https_context = ssl._create_unverified_context

import os
import re
import time
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO


# ============================================================
# SETTINGS
# ============================================================

INPUT_CSV = "raw_vfdb_frequency.csv"
OUTPUT_CSV = "raw_vfdb_frequency_with_locus_tag.csv"

# Local PAO1 reference files
PAO1_PROTEIN_FASTA = "GCF_000006765.1_ASM676v1_protein.faa"
PAO1_GFF = "GCF_000006765.1_ASM676v1_genomic.gff"

# Temporary BLAST files
TMP_QUERY_FASTA = "raw_vfdb_unmapped_proteins.faa"
BLAST_DB_PREFIX = "PAO1_protein_blastdb"
BLAST_OUT = "raw_vfdb_unmapped_vs_PAO1_blastp.tsv"

Entrez.email = "boxuanl4@uci.edu"

# Optional:
# Entrez.api_key = "YOUR_NCBI_API_KEY"

MIN_PIDENT = 80.0
MIN_QCOV = 80.0
MAX_EVALUE = 1e-10

BATCH_SIZE = 50
SLEEP_SEC = 0.34


# ============================================================
# BASIC HELPERS
# ============================================================

def extract_protein_accession(value):
    """
    Examples:
        VFG000122(gb|NP_252230) -> NP_252230
        VFG000123(gb|WP_012614635.1) -> WP_012614635.1
    """
    if pd.isna(value):
        return None

    text = str(value)

    match = re.search(r"gb\|([^)\s]+)", text)

    if match:
        return match.group(1)

    match = re.search(
        r"\b([A-Z]{2,3}_\d+(?:\.\d+)?)\b",
        text,
    )

    if match:
        return match.group(1)

    return None


def accession_no_version(accession):
    if accession is None or pd.isna(accession):
        return None

    return str(accession).split(".")[0]


def is_valid_text(value):
    return (
        pd.notna(value)
        and str(value).strip() != ""
        and str(value).lower() != "nan"
    )


# ============================================================
# PAO1 REFERENCE PARSING
# ============================================================

def extract_pao1_locus_tag_from_header(header):
    if header is None:
        return None

    patterns = [
        r"\[locus_tag=(PA\d+)\]",
        r"locus_tag=(PA\d+)",
        r"\[old_locus_tag=(PA\d+)\]",
        r"old_locus_tag=(PA\d+)",
        r"\b(PA\d{4,5})\b",
    ]

    for pattern in patterns:
        match = re.search(pattern, str(header))

        if match:
            return match.group(1)

    return None


def parse_pao1_fasta_headers(pao1_fasta):
    """
    Build mappings from PAO1 protein accession to FASTA header.
    """
    mapping = {}

    for record in SeqIO.parse(pao1_fasta, "fasta"):
        accession_full = record.id
        accession_base = accession_no_version(record.id)

        value = {
            "blast_hit_locus_tag_from_header":
                extract_pao1_locus_tag_from_header(record.description),
            "blast_hit_header": record.description,
        }

        mapping[accession_full] = value
        mapping[accession_base] = value

    print(
        f"Parsed PAO1 FASTA protein records: "
        f"{len(mapping)}"
    )

    return mapping


def parse_pao1_gff(gff_file):
    """
    Build:
        PAO1 protein_id -> PAO1 locus_tag
    """
    mapping = {}

    with open(gff_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9:
                continue

            if parts[2] != "CDS":
                continue

            attributes = {}

            for item in parts[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = value

            locus_tag = attributes.get("locus_tag")
            protein_id = attributes.get("protein_id")

            if locus_tag and protein_id:
                mapping[protein_id] = locus_tag
                mapping[
                    accession_no_version(protein_id)
                ] = locus_tag

    print(
        f"Parsed PAO1 GFF protein-to-locus mappings: "
        f"{len(mapping)}"
    )

    return mapping


# ============================================================
# NCBI DIRECT RETRIEVAL
# ============================================================

def fetch_ncbi_records(
    accessions,
    batch_size=50,
    sleep_sec=0.34,
):
    """
    Return:
        accession without version ->
        NCBI annotation and protein sequence
    """
    mapping = {}

    accessions = sorted(
        {
            accession
            for accession in accessions
            if is_valid_text(accession)
        }
    )

    failed_batches = 0

    for start in range(
        0,
        len(accessions),
        batch_size,
    ):
        batch = accessions[
            start:start + batch_size
        ]

        print(
            f"Fetching {start + 1}-"
            f"{start + len(batch)} / "
            f"{len(accessions)}"
        )

        handle = None

        try:
            handle = Entrez.efetch(
                db="protein",
                id=",".join(batch),
                rettype="gb",
                retmode="text",
            )

            records_found = 0

            for record in SeqIO.parse(
                handle,
                "genbank",
            ):
                records_found += 1

                accession_full = record.id
                accession_base = accession_no_version(
                    record.id
                )

                locus_tag = None
                old_locus_tag = None
                gene = None
                product = record.description

                for feature in record.features:
                    qualifiers = feature.qualifiers

                    if (
                        "locus_tag" in qualifiers
                        and locus_tag is None
                    ):
                        locus_tag = qualifiers[
                            "locus_tag"
                        ][0]

                    if (
                        "old_locus_tag" in qualifiers
                        and old_locus_tag is None
                    ):
                        old_locus_tag = ";".join(
                            qualifiers["old_locus_tag"]
                        )

                    if (
                        "gene" in qualifiers
                        and gene is None
                    ):
                        gene = qualifiers["gene"][0]

                    if (
                        "product" in qualifiers
                        and not is_valid_text(product)
                    ):
                        product = qualifiers[
                            "product"
                        ][0]

                mapping[accession_base] = {
                    "protein_accession_ncbi":
                        accession_full,
                    "locus_tag_direct":
                        locus_tag,
                    "old_locus_tag_direct":
                        old_locus_tag,
                    "ncbi_gene":
                        gene,
                    "ncbi_product":
                        product,
                    "protein_sequence":
                        str(record.seq),
                }

            print(
                f"  Returned records: "
                f"{records_found}"
            )

        except Exception as error:
            failed_batches += 1

            print(
                f"Error fetching batch: {error}"
            )

        finally:
            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass

        time.sleep(sleep_sec)

    print(
        f"NCBI retrieval complete. "
        f"Mapped accessions: {len(mapping)}; "
        f"failed batches: {failed_batches}"
    )

    return mapping


def choose_direct_locus_tag(row):
    locus_tag = row.get(
        "locus_tag_direct"
    )

    old_locus_tag = row.get(
        "old_locus_tag_direct"
    )

    if is_valid_text(locus_tag):
        return str(locus_tag).strip()

    if is_valid_text(old_locus_tag):
        values = [
            value.strip()
            for value in str(
                old_locus_tag
            ).split(";")
            if value.strip()
        ]

        pa_values = [
            value
            for value in values
            if re.fullmatch(
                r"PA\d+",
                value,
            )
        ]

        if pa_values:
            return pa_values[0]

        if values:
            return values[0]

    return None


# ============================================================
# BLAST FALLBACK
# ============================================================

def write_unmapped_query_fasta(
    dataframe,
    output_fasta,
):
    subset = dataframe[
        dataframe[
            "locus_tag_direct_final"
        ].isna()
        & dataframe[
            "protein_sequence"
        ].notna()
        & (
            dataframe[
                "protein_sequence"
            ]
            .astype(str)
            .str.len()
            > 0
        )
    ].copy()

    subset = subset.drop_duplicates(
        "protein_accession_base"
    )

    with open(output_fasta, "w") as output:
        for _, row in subset.iterrows():
            accession = row[
                "protein_accession_base"
            ]

            sequence = str(
                row["protein_sequence"]
            )

            if not is_valid_text(accession):
                continue

            output.write(
                f">{accession}\n"
            )

            for start in range(
                0,
                len(sequence),
                60,
            ):
                output.write(
                    sequence[
                        start:start + 60
                    ]
                    + "\n"
                )

    print(
        f"Sequences written for BLAST: "
        f"{len(subset)}"
    )

    return len(subset)


def make_blast_db_if_needed(
    pao1_fasta,
    database_prefix,
):
    required_files = [
        database_prefix + ".pin",
        database_prefix + ".phr",
        database_prefix + ".psq",
    ]

    if all(
        os.path.exists(path)
        for path in required_files
    ):
        print("PAO1 BLAST database already exists.")
        return

    command = [
        "makeblastdb",
        "-in",
        pao1_fasta,
        "-dbtype",
        "prot",
        "-out",
        database_prefix,
    ]

    print("Running:")
    print(" ".join(command))

    subprocess.run(
        command,
        check=True,
    )


def run_blastp(
    query_fasta,
    database_prefix,
    blast_output,
):
    if (
        not os.path.exists(query_fasta)
        or os.path.getsize(query_fasta) == 0
    ):
        print(
            "No sequences available for BLAST."
        )
        return False

    command = [
        "blastp",
        "-query",
        query_fasta,
        "-db",
        database_prefix,
        "-out",
        blast_output,
        "-outfmt",
        (
            "6 qseqid sseqid pident length "
            "mismatch gapopen qstart qend "
            "sstart send evalue bitscore "
            "qlen slen"
        ),
        "-max_target_seqs",
        "5",
        "-evalue",
        str(MAX_EVALUE),
        "-num_threads",
        "4",
    ]

    print("Running:")
    print(" ".join(command))

    subprocess.run(
        command,
        check=True,
    )

    return True


def get_pao1_locus_from_blast_hit(
    hit_accession,
    pao1_gff_map,
    pao1_header_map,
):
    if not is_valid_text(hit_accession):
        return None

    accession_full = str(hit_accession)
    accession_base = accession_no_version(
        hit_accession
    )

    if accession_full in pao1_gff_map:
        return pao1_gff_map[
            accession_full
        ]

    if accession_base in pao1_gff_map:
        return pao1_gff_map[
            accession_base
        ]

    if accession_full in pao1_header_map:
        tag = pao1_header_map[
            accession_full
        ].get(
            "blast_hit_locus_tag_from_header"
        )

        if is_valid_text(tag):
            return tag

    if accession_base in pao1_header_map:
        tag = pao1_header_map[
            accession_base
        ].get(
            "blast_hit_locus_tag_from_header"
        )

        if is_valid_text(tag):
            return tag

    return None


def get_pao1_header_from_blast_hit(
    hit_accession,
    pao1_header_map,
):
    if not is_valid_text(hit_accession):
        return None

    accession_full = str(hit_accession)
    accession_base = accession_no_version(
        hit_accession
    )

    if accession_full in pao1_header_map:
        return pao1_header_map[
            accession_full
        ].get("blast_hit_header")

    if accession_base in pao1_header_map:
        return pao1_header_map[
            accession_base
        ].get("blast_hit_header")

    return None


def parse_blast_results(
    blast_output,
    pao1_header_map,
    pao1_gff_map,
):
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    ]

    output_columns = [
        "protein_accession_base",
        "blast_hit_accession",
        "blast_hit_locus_tag",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "blast_hit_bitscore",
        "blast_hit_header",
    ]

    if (
        not os.path.exists(blast_output)
        or os.path.getsize(blast_output) == 0
    ):
        return pd.DataFrame(
            columns=output_columns
        )

    blast = pd.read_csv(
        blast_output,
        sep="\t",
        names=columns,
    )

    blast["qcov"] = (
        blast["length"]
        / blast["qlen"]
        * 100
    )

    blast = blast[
        (blast["pident"] >= MIN_PIDENT)
        & (blast["qcov"] >= MIN_QCOV)
        & (blast["evalue"] <= MAX_EVALUE)
    ].copy()

    if blast.empty:
        return pd.DataFrame(
            columns=output_columns
        )

    blast = blast.sort_values(
        [
            "qseqid",
            "evalue",
            "bitscore",
            "pident",
            "qcov",
        ],
        ascending=[
            True,
            True,
            False,
            False,
            False,
        ],
    )

    best = blast.drop_duplicates(
        "qseqid"
    ).copy()

    best["blast_hit_accession"] = (
        best["sseqid"]
    )

    best["blast_hit_locus_tag"] = (
        best["blast_hit_accession"]
        .apply(
            lambda accession:
                get_pao1_locus_from_blast_hit(
                    accession,
                    pao1_gff_map,
                    pao1_header_map,
                )
        )
    )

    best["blast_hit_header"] = (
        best["blast_hit_accession"]
        .apply(
            lambda accession:
                get_pao1_header_from_blast_hit(
                    accession,
                    pao1_header_map,
                )
        )
    )

    best = best.rename(
        columns={
            "qseqid":
                "protein_accession_base",
            "pident":
                "blast_hit_pident",
            "qcov":
                "blast_hit_qcov",
            "evalue":
                "blast_hit_evalue",
            "bitscore":
                "blast_hit_bitscore",
        }
    )

    return best[output_columns]


# ============================================================
# FINAL ASSIGNMENT
# ============================================================

def determine_final_locus_tag(row):
    direct = row.get(
        "locus_tag_direct_final"
    )

    blast = row.get(
        "blast_hit_locus_tag"
    )

    if is_valid_text(direct):
        return str(direct).strip()

    if is_valid_text(blast):
        return str(blast).strip()

    return None


def determine_mapping_method(row):
    direct = row.get(
        "locus_tag_direct_final"
    )

    blast_tag = row.get(
        "blast_hit_locus_tag"
    )

    blast_hit = row.get(
        "blast_hit_accession"
    )

    protein_accession = row.get(
        "protein_accession"
    )

    if is_valid_text(direct):
        return "direct_ncbi_record"

    if is_valid_text(blast_tag):
        return "blastp_best_hit_to_PAO1"

    if is_valid_text(blast_hit):
        return "blast_hit_no_locus_tag"

    if not is_valid_text(protein_accession):
        return "no_protein_accession"

    return "unmapped"


# ============================================================
# MAIN
# ============================================================

def main():
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(
            f"Input CSV not found: {INPUT_CSV}"
        )

    if not os.path.exists(
        PAO1_PROTEIN_FASTA
    ):
        raise FileNotFoundError(
            f"PAO1 protein FASTA not found: "
            f"{PAO1_PROTEIN_FASTA}"
        )

    if not os.path.exists(PAO1_GFF):
        raise FileNotFoundError(
            f"PAO1 GFF not found: "
            f"{PAO1_GFF}"
        )

    dataframe = pd.read_csv(INPUT_CSV)

    if dataframe.empty:
        raise ValueError(
            "Input CSV contains no rows."
        )

    # The first input column contains the VFDB accessions
    accession_column = dataframe.columns[0]

    print(
        f"Using first column: "
        f"{accession_column}"
    )
    print(
        f"Input rows: "
        f"{len(dataframe)}"
    )

    dataframe["protein_accession"] = (
        dataframe[accession_column]
        .apply(
            extract_protein_accession
        )
    )

    dataframe[
        "protein_accession_base"
    ] = (
        dataframe["protein_accession"]
        .apply(accession_no_version)
    )

    unique_accessions = (
        dataframe["protein_accession"]
        .dropna()
        .unique()
        .tolist()
    )

    print(
        f"Unique protein accessions: "
        f"{len(unique_accessions)}"
    )

    # Parse local PAO1 references
    print(
        "Parsing PAO1 FASTA headers..."
    )

    pao1_header_map = (
        parse_pao1_fasta_headers(
            PAO1_PROTEIN_FASTA
        )
    )

    print(
        "Parsing PAO1 GFF..."
    )

    pao1_gff_map = parse_pao1_gff(
        PAO1_GFF
    )

    # Direct NCBI mapping
    ncbi_mapping = fetch_ncbi_records(
        unique_accessions,
        batch_size=BATCH_SIZE,
        sleep_sec=SLEEP_SEC,
    )

    mapping_dataframe = (
        pd.DataFrame.from_dict(
            ncbi_mapping,
            orient="index",
        )
        .reset_index()
        .rename(
            columns={
                "index":
                    "protein_accession_base"
            }
        )
    )

    # Ensure expected columns exist even if NCBI returned nothing
    expected_mapping_columns = [
        "protein_accession_base",
        "protein_accession_ncbi",
        "locus_tag_direct",
        "old_locus_tag_direct",
        "ncbi_gene",
        "ncbi_product",
        "protein_sequence",
    ]

    for column in expected_mapping_columns:
        if column not in mapping_dataframe.columns:
            mapping_dataframe[column] = None

    mapping_dataframe = mapping_dataframe[
        expected_mapping_columns
    ]

    output = dataframe.merge(
        mapping_dataframe,
        on="protein_accession_base",
        how="left",
    )

    output[
        "locus_tag_direct_final"
    ] = output.apply(
        choose_direct_locus_tag,
        axis=1,
    )

    print(
        f"Direct NCBI locus tags: "
        f"{output['locus_tag_direct_final'].notna().sum()} "
        f"/ {len(output)}"
    )

    # BLAST fallback
    query_count = write_unmapped_query_fasta(
        output,
        TMP_QUERY_FASTA,
    )

    blast_best = pd.DataFrame(
        columns=[
            "protein_accession_base",
            "blast_hit_accession",
            "blast_hit_locus_tag",
            "blast_hit_pident",
            "blast_hit_qcov",
            "blast_hit_evalue",
            "blast_hit_bitscore",
            "blast_hit_header",
        ]
    )

    if query_count > 0:
        make_blast_db_if_needed(
            PAO1_PROTEIN_FASTA,
            BLAST_DB_PREFIX,
        )

        blast_ran = run_blastp(
            TMP_QUERY_FASTA,
            BLAST_DB_PREFIX,
            BLAST_OUT,
        )

        if blast_ran:
            blast_best = parse_blast_results(
                BLAST_OUT,
                pao1_header_map,
                pao1_gff_map,
            )

    output = output.merge(
        blast_best,
        on="protein_accession_base",
        how="left",
    )

    output["locus_tag"] = output.apply(
        determine_final_locus_tag,
        axis=1,
    )

    output["mapping_method"] = output.apply(
        determine_mapping_method,
        axis=1,
    )

    # Move locus_tag into the third column
    locus_tag_values = output.pop(
        "locus_tag"
    )

    insertion_position = min(
        2,
        len(output.columns),
    )

    output.insert(
        insertion_position,
        "locus_tag",
        locus_tag_values,
    )

    output.to_csv(
        OUTPUT_CSV,
        index=False,
    )

    print("")
    print(
        f"Saved: {OUTPUT_CSV}"
    )

    print("")
    print("Mapping summary:")
    print(
        output["mapping_method"]
        .value_counts(
            dropna=False
        )
    )

    print("")
    print(
        f"Final mapped locus tags: "
        f"{output['locus_tag'].notna().sum()} "
        f"/ {len(output)}"
    )

    print("")
    print(
        "First output columns:"
    )
    print(
        output.columns[:8].tolist()
    )


if __name__ == "__main__":
    main()