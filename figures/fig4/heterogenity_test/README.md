# vgrG Acute vs Chronic Heterogeneity Analysis

This script compares nucleotide sequence heterogeneity of `vgrG` genes between **acute** and **chronic** samples in Clade A.

For each matching FASTA alignment, it calculates:

- Number of sequences
- Alignment length
- Variable sites
- Unique haplotypes
- Mean pairwise nucleotide distance
- Mean distance to acute and chronic consensus sequences

It also reports simple flags indicating whether the acute group is more heterogeneous and whether the acute consensus is less representative.

## Output

The main output is:

`vgrG_acute_vs_chronic_heterogeneity_summary.csv`

containing all heterogeneity and consensus-distance statistics for each `vgrG` alignment.
