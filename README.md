# dN/dS Calculation Script

This script calculates **dN/dS (nonsynonymous / synonymous substitution rate)** from codon-aligned FASTA files using either **pairwise** or **groupwise** comparisons.

> **Important:** First, run the `locus_tag` mapping/filtering step on `core_VF.csv` to identify and remove redundant VF entries that map to the same locus tag. Then run the main workflow to generate the `dnds_output` directory and all required dN/dS result files before performing any downstream analyses in the `figures` folder.

---

## Modes

### pairwise

Computes dN/dS between all pairs of sequences in a codon-aligned FASTA file.

### groupwise

Computes dN/dS between each sequence in an alignment and a reference sequence.

The reference can be:

* a single reference FASTA
* the most frequent sequence from a reference alignment using `--fast`
* a consensus sequence from a reference alignment using `--consensus`
* an ancestral sequence reconstructed with IQ-TREE if neither `--fast` nor `--consensus` is specified

---

## Requirements

* Python 3
* Biopython
* tqdm
* IQ-TREE 2.x, only when ancestral sequence reconstruction is used

Install Python dependencies with:

`pip install biopython tqdm`

---

## Usage

`python two_mode_test.py <mode> [arguments]`

### 1. Pairwise Mode

`python two_mode_test.py pairwise input.fasta -o OUTPUT_DIR -t THREADS --format long`

Arguments:

* `input.fasta`: codon-aligned FASTA file
* `--format`: output format, either `long` or `matrix`
* `-o`, `--output_dir`: output directory
* `-t`, `--threads`: number of parallel worker processes
* `--log`: optional directory for the error log

Example:

`python two_mode_test.py pairwise example_alignment.fasta -o results/ --format long -t 8`

Main output:

`*_pairwise_dnds.csv`

Long-format columns:

`Seq1    Seq2    dN    dS    dN/dS`

---

### 2. Groupwise Mode

`python two_mode_test.py groupwise alignment.fasta reference.fasta -o OUTPUT_DIR -t THREADS`

Arguments:

* `alignment.fasta`: codon-aligned sequences to evaluate
* `reference.fasta`: either a single reference sequence or a reference alignment
* `--fast`: use the most frequent sequence from the reference alignment
* `--consensus`: use the consensus sequence from the reference alignment
* `-o`, `--output_dir`: output directory
* `-t`, `--threads`: number of worker processes / IQ-TREE threads
* `--log`: optional directory for the error log

Example using a single reference:

`python two_mode_test.py groupwise alignment.fasta reference.fasta -o results/ -t 8`

Example using a consensus reference:

`python two_mode_test.py groupwise alignment.fasta reference_alignment.fasta -o results/ --consensus -t 8`

Example using the most frequent sequence:

`python two_mode_test.py groupwise alignment.fasta reference_alignment.fasta -o results/ --fast -t 8`

If neither `--fast` nor `--consensus` is provided and the reference contains multiple sequences, IQ-TREE is used for ancestral sequence reconstruction.

Main output:

`*_groupwise_dnds.csv`

Columns:

`Ref    Seq    dN    dS    dN/dS`

---

## Reference Methods

### Most Frequent Sequence

Using `--fast` selects the sequence that occurs most frequently in the reference alignment.

### Consensus Sequence

Using `--consensus` generates a consensus sequence from the reference alignment.

### Ancestral Sequence

If neither `--fast` nor `--consensus` is used with a multi-sequence reference alignment, the script runs IQ-TREE ancestral state reconstruction and uses the reconstructed root sequence as the reference.

---

## Sequence Processing

Before dN/dS calculation, the script:

* reverse-complements sequences marked with `_R_`
* removes shared gap codons
* resolves ambiguous bases using the paired reference sequence
* removes stop codons
* removes codons containing gaps
* maintains the codon reading frame

---

## Output Summary

Pairwise mode:

`*_pairwise_dnds.csv`

Groupwise mode:

`*_groupwise_dnds.csv`

Reference FASTA files may also be generated when using consensus, most-frequent-sequence, or ancestral reference inference.

---

## Notes

* Input sequences should already be codon-aligned.
* dN and dS are calculated using the NG86 method.
* If `dS = 0` and `dN > 0`, the reported dN/dS ratio is `inf`.
* If both `dN = 0` and `dS = 0`, the reported ratio is `0`.
* Failed sequence comparisons are written to the error log rather than stopping the entire analysis.
