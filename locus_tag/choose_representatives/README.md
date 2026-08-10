# Representative VF Selection

This script selects one representative VF entry for each repeated locus tag using genome coverage and codon-alignment quality.

It is intended to remove redundant VFDB entries that map to the same gene before downstream dN/dS analysis.

## Input

The script requires:

- `core_VF_labeled_with_locus_tag.csv`
- the corresponding codon-aligned FASTA files from the workflow `final_aln` directory

Expected alignment structure:

`final_aln/<VFG_ID>.nogap.dedup/<VFG_ID>.nogap.dedup.nt_ali.fasta`

The summary CSV must contain at least:

- `Accession`
- `locus_tag`
- `genome_count`

:contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}

## Selection Criteria

For locus tags represented by multiple VF entries, the script ranks candidates using the following priority:

1. Successful alignment analysis
2. Highest genome count
3. Highest number of usable aligned sequences
4. Highest usable-sequence fraction
5. Highest usable-codon fraction
6. Lowest gap fraction
7. Lowest ambiguous-base fraction
8. Longest protein as the final tie-breaker

PAO1 similarity is not used for representative selection. :contentReference[oaicite:2]{index=2}

A sequence is considered usable when no more than 20% of aligned positions contain gaps or non-ATGC characters. A codon position is considered usable when at least 80% of sequences contain an unambiguous ATGC codon at that position. :contentReference[oaicite:3]{index=3}

## Main Outputs

`alignment_quality_metrics.csv`  
Alignment-quality statistics for all VF entries in repeated locus-tag groups.

`alignment_quality_failures.csv`  
VF entries with missing or failed alignment analysis.

`repeated_locus_all_VF_quality_rankings.csv`  
All redundant VF candidates together with their quality metrics and representative ranking.

`selected_representative_VF_per_locus.csv`  
The selected representative for each repeated locus tag.

`excluded_redundant_VF_entries.csv`  
Redundant VF entries that were not selected.

`core_VF_nonredundant_representatives.csv`  
Final nonredundant core VF table containing unique-locus VFs plus the selected representatives from repeated locus groups.

`representative_VFG_ID_list.txt`  
List of VFG IDs retained in the final nonredundant set.

:contentReference[oaicite:4]{index=4}

## Usage

Run with the default number of worker processes:

`python representative_VF_selection.py`

Or specify the number of workers:

`python representative_VF_selection.py --workers 16`

> **Important:** Run the main workflow first so that the required `final_aln` codon alignments have been generated before running this script.
