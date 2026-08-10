# VFDB Locus Tag Mapping

This script maps VFDB protein accessions to *Pseudomonas aeruginosa* PAO1 locus tags.

It first retrieves locus tag annotations directly from NCBI protein records. For entries that cannot be mapped directly, it uses BLASTP against the PAO1 proteome and converts the best matching PAO1 protein to its locus tag using the PAO1 GFF annotation. :contentReference[oaicite:0]{index=0}

## Requirements

- Python 3
- pandas
- Biopython
- NCBI BLAST+

## Input

- VFDB annotation CSV
- PAO1 protein FASTA
- PAO1 GFF

Main output:

`raw_vfdb_frequency_with_locus_tag.csv`

This file contains the original VFDB information together with the mapped `locus_tag` and the mapping method (`direct_ncbi_record` or `blastp_best_hit_to_PAO1`).
