
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
* iqtree 2.4.0

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

#### **groupwise**

```bash
python script.py groupwise alignment.fasta reference.fasta [-o OUTPUT_DIR] [-t THREADS]
```

* `alignment.fasta`: Input codon alignment FASTA file.
* `reference.fasta`: Either a single-sequence reference or an alignment used to infer the reference.
* `-o`, `--output_dir`: Output directory (default: current directory).
* `-t`, `--threads`: Number of threads for parallel processing (default: 1).

---

## Output

* Pairwise: CSV file containing all pairwise dN/dS values.
* Groupwise: CSV file with dN/dS values comparing each sequence to the reference.

---

## Example Commands

**Pairwise comparison:**

```bash
python script.py pairwise example_alignment.fasta -o results/ --format long
```

**Groupwise comparison:**

```bash
python script.py groupwise example_alignment.fasta ref_seq.fasta -o results/
```

---


