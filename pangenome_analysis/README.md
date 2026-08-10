# VF Presence/Absence and Pangenome Analysis

This folder contains four scripts used to calculate VF frequencies, map VF entries to PAO1 locus tags, remove redundant VF entries, and test for differences in VF presence/absence between Clade A and Clade B.

## Run Order

### 1. raw_vfdb_frequency.R

Run this script first.

It uses the VFDB annotation and genome metadata to calculate the frequency of each VF in Clade A and Clade B. 

Main output:

raw_vfdb_frequency.csv

---

### 2. locus_tag.py

Run this script second.

Input:

raw_vfdb_frequency.csv

The script maps VFDB protein accessions to PAO1 locus tags. It first attempts direct mapping using NCBI protein annotations. Entries that cannot be mapped directly are mapped by BLASTP against the PAO1 proteome, using the PAO1 GFF file to recover the corresponding locus tag.  

Required reference files:

GCF_000006765.1_ASM676v1_protein.faa
GCF_000006765.1_ASM676v1_genomic.gff

Main output:

raw_vfdb_frequency_with_locus_tag.csv 

---

### 3. raw_vfdb_frequency_nonredundant.R

Run this script third.

Input:

raw_vfdb_frequency_with_locus_tag.csv

Some different VFDB accessions map to the same locus tag. This script removes these redundant entries by keeping one representative VF for each mapped locus.

The representative is primarily selected based on the highest combined genome coverage across the two clades. VFs without a mapped locus tag are retained individually. 

Main outputs:

VF_representative_list_one_per_locus.csv
raw_vfdb_frequency_nonredundant.csv
VF_discarded_repeated_locus_entries.csv 

---

### 4. pangenome_analysis.R

Run this script last.

Input:

raw_vfdb_frequency_nonredundant.csv

This script compares VF presence/absence between Clade A and Clade B using Fisher's exact test.

It calculates:

* VF frequency in each clade
* frequency difference between clades
* odds ratio
* Fisher's exact test p-value
* Benjamini-Hochberg FDR
* which clade a VF is enriched in

A VF is considered meaningfully enriched when FDR ≤ 0.05 and the absolute frequency difference between clades is at least 0.05. 

Main outputs:

VF_presence_absence_enrichment_Clade_A_vs_B_all.csv
VF_presence_enriched_in_Clade_A.csv
VF_presence_enriched_in_Clade_B.csv
VF_statistically_significant_small_difference.csv
VF_enrichment_without_locus_tag.csv 

## Workflow

raw_vfdb_frequency.R
↓
raw_vfdb_frequency.csv
↓
locus_tag.py
↓
raw_vfdb_frequency_with_locus_tag.csv
↓
raw_vfdb_frequency_nonredundant.R
↓
raw_vfdb_frequency_nonredundant.csv
↓
pangenome_analysis.R
↓
VF presence/absence enrichment results

