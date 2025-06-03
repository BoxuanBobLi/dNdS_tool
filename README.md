
---

# dN/dS Calculation Script

This script calculates **dN/dS (nonsynonymous to synonymous substitution rate) ratios** using codon-aligned FASTA files. It supports pairwise and groupwise analyses for evolutionary comparisons.

## Modes

* **pairwise**: Computes dN/dS ratios between all pairs of sequences in a codon-aligned FASTA file.
* **groupwise**: Computes dN/dS ratios between all sequences and a reference sequence (either from a separate reference file or extracted from a multiple sequence alignment).

---

## Requirements

* Python 3.6+
* Biopython (`pip install biopython`)
* iqtree 2.4.0 (if neither `--fast` nor `--consensus` mode is used in `groupwise` mode)

---

## Usage

```bash
python script.py <mode> [arguments]
```

### Modes and Arguments

#### **pairwise**

```bash
python script.py pairwise input.fasta [-o OUTPUT_DIR] [-t THREADS] [--format FORMAT]
```

* `input.fasta`: Input codon alignment FASTA file.
* `--format`: Output format: `long` (default) or `matrix`.
* `-o`, `--output_dir`: Output directory (default: current directory).
* `-t`, `--threads`: Number of threads for parallel processing (default: 1).
* `--log`: Directory to save log file (default: output_dir).

#### **groupwise**

```bash
python script.py groupwise alignment.fasta reference.fasta [-o OUTPUT_DIR] [-t THREADS]
```

* `alignment.fasta`: Input codon alignment FASTA file.
* `reference.fasta`: Either a single-sequence reference or an alignment used to infer the reference.
* `-o`, `--output_dir`: Output directory (default: current directory).
* `-t`, `--threads`: Number of threads for parallel processing (default: 1).
* `--fast`: Use most frequent sequence as reference from MSA
* `--consensus`: Use consensus sequence as reference from MSA
* `--log`: Directory to save log file (default: output_dir).

---

## Output

* Pairwise: CSV file containing all pairwise dN/dS values and error log file
* Groupwise: CSV file with dN/dS values comparing each sequence to the reference, reference sequence generated, and error log file.

---

## Example Commands

**Pairwise comparison:**

```bash
python script.py pairwise example_alignment.fasta -o results/ --format long
```

**Groupwise comparison:**

```bash
python script.py groupwise example_alignment.fasta ref_seq.fasta -o results/ --consensus
```

---


