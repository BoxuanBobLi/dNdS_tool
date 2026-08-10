#!/usr/bin/env python3

import os
import re
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm


# ============================================================
# PATHS
# ============================================================

SUMMARY_CSV = (
    "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/"
    "locus_check/core_VF_labeled_with_locus_tag.csv"
)

ALIGNMENT_DIR = (
    "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/"
    "final_aln"
)

OUTPUT_DIR = (
    "/data2/B_Li/vfdb/workflow_clade_translatorX_solved_dup_solved/"
    "locus_check/pick_VFs/output"
)


# ============================================================
# FILTERING SETTINGS
# ============================================================

# A sequence is usable if no more than 20% of aligned positions
# contain gaps or non-ATGC characters.
MAX_BAD_SEQUENCE_FRACTION = 0.20

# A codon column is usable if at least 80% of sequences contain
# an unambiguous ATGC codon at that position.
MIN_CODON_USABLE_PROPORTION = 0.80


# ============================================================
# HELPERS
# ============================================================

def extract_vfg_id(value):
    """
    Extract VFG ID such as VFG000115.
    """

    match = re.search(
        r"(VFG\d+)",
        str(value),
        re.IGNORECASE,
    )

    return match.group(1).upper() if match else None


def find_alignment(vfg_id):
    """
    Find the final TranslatorX nucleotide alignment for one VF.
    """

    expected_path = os.path.join(
        ALIGNMENT_DIR,
        f"{vfg_id}.nogap.dedup",
        f"{vfg_id}.nogap.dedup.nt_ali.fasta",
    )

    if os.path.isfile(expected_path):
        return expected_path

    # Fallback in case the directory structure varies slightly.
    matches = glob.glob(
        os.path.join(
            ALIGNMENT_DIR,
            f"{vfg_id}.nogap.dedup",
            f"{vfg_id}*.nt_ali.fasta",
        )
    )

    if len(matches) == 1:
        return matches[0]

    return None


def empty_quality_result(
    vfg_id,
    alignment_path=None,
    error_message="",
):
    """
    Return an empty result for missing or failed alignments.
    """

    return {
        "VFG_ID": vfg_id,
        "alignment_path": alignment_path,
        "alignment_found": False,
        "analysis_success": False,
        "error_message": error_message,
        "n_sequences_total": 0,
        "n_sequences_usable": 0,
        "proportion_sequences_usable": np.nan,
        "alignment_length_nt": np.nan,
        "alignment_length_codons": np.nan,
        "trailing_non_codon_nt": np.nan,
        "gap_fraction": np.nan,
        "ambiguous_fraction": np.nan,
        "bad_character_fraction": np.nan,
        "usable_codon_count": np.nan,
        "usable_codon_fraction": np.nan,
    }


def read_fasta_sequences(path):
    """
    Read FASTA sequences using SimpleFastaParser, which is faster
    than constructing full SeqRecord objects with SeqIO.parse.
    """

    sequences = []

    with open(path, "r") as handle:

        for _, sequence in SimpleFastaParser(handle):

            cleaned = (
                sequence
                .replace(" ", "")
                .replace("\n", "")
                .replace("\r", "")
                .upper()
            )

            sequences.append(cleaned)

    return sequences


def alignment_quality_worker(task):
    """
    Analyze one codon alignment.

    This function runs inside a worker process, so it must remain
    defined at module level.
    """

    vfg_id, alignment_path = task

    if alignment_path is None:

        return empty_quality_result(
            vfg_id=vfg_id,
            alignment_path=None,
            error_message="Alignment file not found",
        )

    try:

        sequences = read_fasta_sequences(
            alignment_path
        )

        if not sequences:

            return empty_quality_result(
                vfg_id=vfg_id,
                alignment_path=alignment_path,
                error_message="Alignment contains no sequences",
            )

        lengths = np.fromiter(
            (len(sequence) for sequence in sequences),
            dtype=np.int64,
        )

        unique_lengths = np.unique(lengths)

        if len(unique_lengths) != 1:

            return empty_quality_result(
                vfg_id=vfg_id,
                alignment_path=alignment_path,
                error_message=(
                    "Sequences have different alignment lengths: "
                    + ",".join(map(str, unique_lengths[:20]))
                ),
            )

        n_sequences = len(sequences)
        alignment_length = int(unique_lengths[0])

        if alignment_length == 0:

            return empty_quality_result(
                vfg_id=vfg_id,
                alignment_path=alignment_path,
                error_message="Alignment length is zero",
            )

        # Convert all sequences into one n_sequences × alignment_length
        # byte matrix. ASCII A, T, G, C, N and '-' are one byte each.
        joined = "".join(sequences).encode("ascii", errors="replace")

        matrix = np.frombuffer(
            joined,
            dtype=np.uint8,
        ).reshape(
            n_sequences,
            alignment_length,
        )

        # ASCII values:
        # A=65, C=67, G=71, T=84, -=45
        is_a = matrix == 65
        is_c = matrix == 67
        is_g = matrix == 71
        is_t = matrix == 84

        is_valid_base = (
            is_a
            | is_c
            | is_g
            | is_t
        )

        is_gap = matrix == 45

        # Ambiguous means anything other than ATGC or a gap.
        is_ambiguous = (
            ~is_valid_base
            & ~is_gap
        )

        is_bad = ~is_valid_base

        total_characters = matrix.size

        gap_count = int(
            np.count_nonzero(is_gap)
        )

        ambiguous_count = int(
            np.count_nonzero(is_ambiguous)
        )

        gap_fraction = (
            gap_count / total_characters
        )

        ambiguous_fraction = (
            ambiguous_count / total_characters
        )

        bad_character_fraction = (
            gap_fraction
            + ambiguous_fraction
        )

        # ----------------------------------------------------
        # Usable sequences
        # ----------------------------------------------------

        bad_count_per_sequence = (
            np.count_nonzero(
                is_bad,
                axis=1,
            )
        )

        bad_fraction_per_sequence = (
            bad_count_per_sequence
            / alignment_length
        )

        usable_sequence_mask = (
            bad_fraction_per_sequence
            <= MAX_BAD_SEQUENCE_FRACTION
        )

        n_sequences_usable = int(
            np.count_nonzero(
                usable_sequence_mask
            )
        )

        proportion_sequences_usable = (
            n_sequences_usable
            / n_sequences
        )

        # ----------------------------------------------------
        # Usable codon columns
        # ----------------------------------------------------

        n_codons = alignment_length // 3
        trailing_non_codon_nt = alignment_length % 3

        if n_codons > 0:

            codon_length = n_codons * 3

            valid_codon_matrix = (
                is_valid_base[:, :codon_length]
                .reshape(
                    n_sequences,
                    n_codons,
                    3,
                )
                .all(axis=2)
            )

            valid_sequences_per_codon = (
                valid_codon_matrix.sum(axis=0)
            )

            usable_proportion_per_codon = (
                valid_sequences_per_codon
                / n_sequences
            )

            usable_codon_mask = (
                usable_proportion_per_codon
                >= MIN_CODON_USABLE_PROPORTION
            )

            usable_codon_count = int(
                np.count_nonzero(
                    usable_codon_mask
                )
            )

            usable_codon_fraction = (
                usable_codon_count
                / n_codons
            )

        else:

            usable_codon_count = 0
            usable_codon_fraction = np.nan

        return {
            "VFG_ID": vfg_id,
            "alignment_path": alignment_path,
            "alignment_found": True,
            "analysis_success": True,
            "error_message": "",
            "n_sequences_total": n_sequences,
            "n_sequences_usable": n_sequences_usable,
            "proportion_sequences_usable": (
                proportion_sequences_usable
            ),
            "alignment_length_nt": alignment_length,
            "alignment_length_codons": n_codons,
            "trailing_non_codon_nt": trailing_non_codon_nt,
            "gap_fraction": gap_fraction,
            "ambiguous_fraction": ambiguous_fraction,
            "bad_character_fraction": (
                bad_character_fraction
            ),
            "usable_codon_count": usable_codon_count,
            "usable_codon_fraction": (
                usable_codon_fraction
            ),
        }

    except Exception as error:

        return empty_quality_result(
            vfg_id=vfg_id,
            alignment_path=alignment_path,
            error_message=repr(error),
        )


def safe_metric(value, digits=4):
    """
    Format metrics safely for the selection-reason text.
    """

    if pd.isna(value):
        return "NA"

    if isinstance(value, (int, np.integer)):
        return str(int(value))

    return f"{float(value):.{digits}f}"


def build_selection_reason(group):
    """
    Explain the selected representative for one locus group.
    """

    selected = group.iloc[0]

    return (
        f"Selected {selected['VFG_ID']}: "
        f"genome_count="
        f"{safe_metric(selected['genome_count'], 0)}, "
        f"usable_sequences="
        f"{safe_metric(selected['n_sequences_usable'], 0)}, "
        f"usable_sequence_fraction="
        f"{safe_metric(selected['proportion_sequences_usable'])}, "
        f"usable_codon_fraction="
        f"{safe_metric(selected['usable_codon_fraction'])}, "
        f"gap_fraction="
        f"{safe_metric(selected['gap_fraction'])}, "
        f"ambiguous_fraction="
        f"{safe_metric(selected['ambiguous_fraction'])}, "
        f"protein_length="
        f"{safe_metric(selected['protein_length'], 0)}"
    )


def parse_arguments():
    """
    Command-line options.
    """

    default_workers = min(
        32,
        os.cpu_count() or 1,
    )

    parser = argparse.ArgumentParser(
        description=(
            "Select one representative VF per repeated locus tag "
            "using genome coverage and codon-alignment quality."
        )
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=default_workers,
        help=(
            "Number of parallel worker processes. "
            f"Default: {default_workers}"
        ),
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=1,
        help=(
            "Task chunk size for multiprocessing. "
            "Default: 1"
        ),
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_arguments()

    workers = max(
        1,
        args.workers,
    )

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    print("=" * 68)
    print("REPRESENTATIVE VF SELECTION")
    print("=" * 68)

    print(f"Summary table:  {SUMMARY_CSV}")
    print(f"Alignment dir:  {ALIGNMENT_DIR}")
    print(f"Output dir:     {OUTPUT_DIR}")
    print(f"Worker count:   {workers}")

    # ========================================================
    # READ SUMMARY TABLE
    # ========================================================

    df = pd.read_csv(
        SUMMARY_CSV
    )

    required_columns = {
        "Accession",
        "locus_tag",
        "genome_count",
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    df["VFG_ID"] = (
        df["Accession"]
        .apply(extract_vfg_id)
    )

    df["locus_tag"] = (
        df["locus_tag"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    df["genome_count"] = pd.to_numeric(
        df["genome_count"],
        errors="coerce",
    )

    if "freq" in df.columns:

        df["freq"] = pd.to_numeric(
            df["freq"],
            errors="coerce",
        )

    if "protein_sequence" in df.columns:

        df["protein_length"] = (
            df["protein_sequence"]
            .fillna("")
            .astype(str)
            .str.replace(
                r"\s+",
                "",
                regex=True,
            )
            .str.len()
        )

    else:

        df["protein_length"] = np.nan

    df = df[
        df["VFG_ID"].notna()
        & df["locus_tag"].ne("")
    ].copy()

    # Avoid accidental repeated rows for one VFG-locus pair.
    df = df.drop_duplicates(
        subset=[
            "VFG_ID",
            "locus_tag",
        ],
        keep="first",
    )

    # ========================================================
    # IDENTIFY REPEATED LOCUS TAGS
    # ========================================================

    locus_counts = (
        df.groupby(
            "locus_tag"
        )["VFG_ID"]
        .transform("nunique")
    )

    repeated = df[
        locus_counts > 1
    ].copy()

    repeated_locus_count = (
        repeated[
            "locus_tag"
        ].nunique()
    )

    repeated_vfg_count = (
        repeated[
            "VFG_ID"
        ].nunique()
    )

    print()
    print(
        f"Repeated locus tags: "
        f"{repeated_locus_count}"
    )

    print(
        f"VF entries in repeated groups: "
        f"{repeated_vfg_count}"
    )

    if repeated.empty:

        raise RuntimeError(
            "No repeated locus-tag groups were found."
        )

    # ========================================================
    # LOCATE ALIGNMENTS
    # ========================================================

    repeated_vfg_ids = sorted(
        repeated[
            "VFG_ID"
        ].unique()
    )

    tasks = []

    for vfg_id in repeated_vfg_ids:

        alignment_path = find_alignment(
            vfg_id
        )

        tasks.append(
            (
                vfg_id,
                alignment_path,
            )
        )

    missing_alignment_count = sum(
        alignment_path is None
        for _, alignment_path in tasks
    )

    print(
        f"Missing alignment files: "
        f"{missing_alignment_count}"
    )

    # ========================================================
    # PARALLEL ALIGNMENT ANALYSIS
    # ========================================================

    quality_rows = []

    print()
    print(
        "Analyzing codon alignments..."
    )

    if workers == 1:

        for task in tqdm(
            tasks,
            total=len(tasks),
            desc="Alignment quality",
            unit="VF",
            dynamic_ncols=True,
        ):

            result = alignment_quality_worker(
                task
            )

            quality_rows.append(
                result
            )

    else:

        with ProcessPoolExecutor(
            max_workers=workers
        ) as executor:

            future_to_vfg = {
                executor.submit(
                    alignment_quality_worker,
                    task,
                ): task[0]
                for task in tasks
            }

            progress = tqdm(
                total=len(future_to_vfg),
                desc="Alignment quality",
                unit="VF",
                dynamic_ncols=True,
            )

            for future in as_completed(
                future_to_vfg
            ):

                vfg_id = future_to_vfg[
                    future
                ]

                try:

                    result = future.result()

                except Exception as error:

                    result = empty_quality_result(
                        vfg_id=vfg_id,
                        alignment_path=None,
                        error_message=(
                            f"Worker failure: {repr(error)}"
                        ),
                    )

                quality_rows.append(
                    result
                )

                progress.set_postfix_str(
                    vfg_id
                )

                progress.update(1)

            progress.close()

    quality_df = pd.DataFrame(
        quality_rows
    )

    quality_df = quality_df.sort_values(
        "VFG_ID"
    )

    quality_df.to_csv(
        os.path.join(
            OUTPUT_DIR,
            "alignment_quality_metrics.csv",
        ),
        index=False,
    )

    failed_quality = quality_df[
        ~quality_df[
            "analysis_success"
        ]
    ].copy()

    failed_quality.to_csv(
        os.path.join(
            OUTPUT_DIR,
            "alignment_quality_failures.csv",
        ),
        index=False,
    )

    print()
    print(
        f"Successful alignment analyses: "
        f"{quality_df['analysis_success'].sum()}"
    )

    print(
        f"Failed or missing alignments: "
        f"{len(failed_quality)}"
    )

    # ========================================================
    # MERGE QUALITY METRICS
    # ========================================================

    scored = repeated.merge(
        quality_df,
        on="VFG_ID",
        how="left",
    )

    # ========================================================
    # RANK WITHIN EACH LOCUS TAG
    # ========================================================
    #
    # Lexicographic priority:
    #
    # 1. Successful alignment analysis
    # 2. Highest genome count
    # 3. Most usable aligned sequences
    # 4. Highest usable-sequence fraction
    # 5. Highest usable-codon fraction
    # 6. Lowest gap fraction
    # 7. Lowest ambiguous fraction
    # 8. Longest protein as final tie-breaker
    #
    # This does not use PAO1 similarity.
    # ========================================================

    scored["analysis_success_rank"] = (
        scored["analysis_success"]
        .fillna(False)
        .astype(int)
    )

    scored["genome_count_rank"] = (
        scored["genome_count"]
        .fillna(-1)
    )

    scored["usable_sequence_count_rank"] = (
        scored["n_sequences_usable"]
        .fillna(-1)
    )

    scored["usable_sequence_fraction_rank"] = (
        scored[
            "proportion_sequences_usable"
        ]
        .fillna(-1)
    )

    scored["usable_codon_fraction_rank"] = (
        scored[
            "usable_codon_fraction"
        ]
        .fillna(-1)
    )

    scored["gap_fraction_rank"] = (
        scored["gap_fraction"]
        .fillna(np.inf)
    )

    scored["ambiguous_fraction_rank"] = (
        scored[
            "ambiguous_fraction"
        ]
        .fillna(np.inf)
    )

    scored["protein_length_rank"] = (
        scored["protein_length"]
        .fillna(-1)
    )

    scored = scored.sort_values(
        [
            "locus_tag",
            "analysis_success_rank",
            "genome_count_rank",
            "usable_sequence_count_rank",
            "usable_sequence_fraction_rank",
            "usable_codon_fraction_rank",
            "gap_fraction_rank",
            "ambiguous_fraction_rank",
            "protein_length_rank",
            "VFG_ID",
        ],
        ascending=[
            True,
            False,
            False,
            False,
            False,
            False,
            True,
            True,
            False,
            True,
        ],
    )

    scored["representative_rank"] = (
        scored.groupby(
            "locus_tag"
        )
        .cumcount()
        + 1
    )

    scored["selected_representative"] = (
        scored[
            "representative_rank"
        ] == 1
    )

    # ========================================================
    # SELECTION REASONS
    # ========================================================

    reason_rows = []

    for locus_tag, group in scored.groupby(
        "locus_tag",
        sort=False,
    ):

        reason_rows.append({
            "locus_tag": locus_tag,
            "selection_reason": (
                build_selection_reason(
                    group
                )
            ),
        })

    reasons = pd.DataFrame(
        reason_rows
    )

    scored = scored.merge(
        reasons,
        on="locus_tag",
        how="left",
    )

    # ========================================================
    # OUTPUT TABLES
    # ========================================================

    output_columns = [
        "locus_tag",
        "representative_rank",
        "selected_representative",
        "VFG_ID",
        "Accession",
        "gene",
        "genome_count",
        "freq",
        "n_sequences_total",
        "n_sequences_usable",
        "proportion_sequences_usable",
        "alignment_length_nt",
        "alignment_length_codons",
        "trailing_non_codon_nt",
        "gap_fraction",
        "ambiguous_fraction",
        "bad_character_fraction",
        "usable_codon_count",
        "usable_codon_fraction",
        "protein_length",
        "mapping_method",
        "blast_hit_pident",
        "blast_hit_qcov",
        "blast_hit_evalue",
        "blast_hit_bitscore",
        "alignment_found",
        "analysis_success",
        "error_message",
        "alignment_path",
        "selection_reason",
    ]

    output_columns = [
        column
        for column in output_columns
        if column in scored.columns
    ]

    scored[
        output_columns
    ].to_csv(
        os.path.join(
            OUTPUT_DIR,
            "repeated_locus_all_VF_quality_rankings.csv",
        ),
        index=False,
    )

    representatives = scored[
        scored[
            "selected_representative"
        ]
    ].copy()

    representatives[
        output_columns
    ].to_csv(
        os.path.join(
            OUTPUT_DIR,
            "selected_representative_VF_per_locus.csv",
        ),
        index=False,
    )

    nonrepresentatives = scored[
        ~scored[
            "selected_representative"
        ]
    ].copy()

    nonrepresentatives[
        output_columns
    ].to_csv(
        os.path.join(
            OUTPUT_DIR,
            "excluded_redundant_VF_entries.csv",
        ),
        index=False,
    )

    # ========================================================
    # CREATE FINAL NONREDUNDANT VF TABLE
    # ========================================================

    group_sizes = (
        df.groupby(
            "locus_tag"
        )["VFG_ID"]
        .transform("nunique")
    )

    unique_locus_rows = df[
        group_sizes == 1
    ].copy()

    selected_ids = set(
        representatives[
            "VFG_ID"
        ]
    )

    selected_repeated_rows = df[
        df[
            "VFG_ID"
        ].isin(selected_ids)
    ].copy()

    final_nonredundant = pd.concat(
        [
            unique_locus_rows,
            selected_repeated_rows,
        ],
        ignore_index=True,
    )

    final_nonredundant = (
        final_nonredundant
        .drop_duplicates(
            subset=[
                "locus_tag",
                "VFG_ID",
            ]
        )
        .sort_values(
            [
                "locus_tag",
                "VFG_ID",
            ]
        )
    )

    final_nonredundant.to_csv(
        os.path.join(
            OUTPUT_DIR,
            "core_VF_nonredundant_representatives.csv",
        ),
        index=False,
    )

    final_nonredundant[
        ["VFG_ID"]
    ].drop_duplicates().to_csv(
        os.path.join(
            OUTPUT_DIR,
            "representative_VFG_ID_list.txt",
        ),
        index=False,
        header=False,
    )

    # ========================================================
    # PRINT SUMMARY
    # ========================================================

    print()
    print("=" * 68)
    print("SELECTION COMPLETE")
    print("=" * 68)

    print(
        f"Repeated locus groups evaluated: "
        f"{scored['locus_tag'].nunique()}"
    )

    print(
        f"Representatives selected: "
        f"{len(representatives)}"
    )

    print(
        f"Redundant VF entries excluded: "
        f"{len(nonrepresentatives)}"
    )

    print(
        f"Unique-locus VF entries retained: "
        f"{unique_locus_rows['VFG_ID'].nunique()}"
    )

    print(
        f"Final nonredundant VF count: "
        f"{final_nonredundant['VFG_ID'].nunique()}"
    )

    print(
        f"Alignment analysis failures: "
        f"{len(failed_quality)}"
    )

    print()
    print("Outputs written to:")
    print(OUTPUT_DIR)


# Multiprocessing requires this guard.
if __name__ == "__main__":
    main()