

# dN/dS & MK Test Calculation Script

This script calculates:

* **dN/dS (nonsynonymous / synonymous substitution rate)** using codon-aligned FASTA files in **pairwise** and **groupwise** modes; and
* **McDonald–Kreitman (MK) test** statistics using an **ingroup codon alignment** and a **reference** (either a single sequence or a consensus from an alignment).

---

## Modes

* **pairwise**
  Computes dN/dS ratios between all pairs of sequences in a codon-aligned FASTA file.

* **groupwise**
  Computes dN/dS ratios between each sequence and a reference sequence (either provided as a single FASTA or inferred from a multiple sequence alignment via `--fast` or `--consensus`).

* **mk**
  Performs the **McDonald–Kreitman (MK) test** using an **ingroup codon alignment** (≥2 sequences) and a **reference** (either a single FASTA or derived from a multiple sequence alignment). Counts **Pn, Ps, Dn, Ds**, computes **Neutrality Index (NI)**, **alpha (α = 1 − NI)**, and **Fisher’s exact test p-value**. A continuity correction is applied when necessary to avoid division by zero.

---

## Requirements

* Python 3.6+
* Biopython (`pip install biopython`)
* SciPy (optional; for Fisher’s exact test in MK mode) (`pip install scipy`)
* IQ-TREE 2.x (only if you use neither `--fast` nor `--consensus` for reference inference in **groupwise** or **mk** modes)

---

## Usage

```bash
python two_mode_test.py <mode> [arguments]
```

### Modes and Arguments

#### 1) **pairwise**

```bash
python two_mode_test.py pairwise input.fasta [-o OUTPUT_DIR] [-t THREADS] [--format FORMAT] [--log LOG_DIR]
```

* `input.fasta`: Codon-aligned FASTA with ≥2 sequences
* `--format`: Output format: `long` (default) or `matrix`
* `-o`, `--output_dir`: Output directory (default: current directory)
* `-t`, `--threads`: Number of worker threads (default: 1)
* `--log`: Directory to save error log (default: output_dir)

**Output**:

* `*_pairwise_dnds.csv` containing pairwise dN, dS, and dN/dS

---

#### 2) **groupwise**

```bash
python two_mode_test.py groupwise alignment.fasta reference.fasta [-o OUTPUT_DIR] [-t THREADS] [--fast | --consensus] [--log LOG_DIR]
```

* `alignment.fasta`: Codon alignment of sequences to evaluate
* `reference.fasta`: Either a **single reference sequence** OR an **alignment** used to derive a reference
* `--fast`: Use the **most frequent sequence** in the reference MSA as the reference
* `--consensus`: Use the **consensus sequence** of the reference MSA as the reference
* (If neither `--fast` nor `--consensus` is provided and the reference is an alignment, the script will perform **ancestral state reconstruction** using **IQ-TREE**)
* `-o`, `--output_dir`: Output directory
* `-t`, `--threads`: Number of threads for reference inference
* `--log`: Directory to save log file (default: output_dir)

**Output**:

* `*_groupwise_dnds.csv` (tab-separated) with columns: `Ref`, `Seq`, `dN`, `dS`, `dN/dS`
* If using `--consensus` or `--fast`, a reference FASTA will be saved alongside

---

#### 3) **mk** (McDonald–Kreitman Test)

```bash
python two_mode_test.py mk ingroup_alignment.fasta reference.fasta [-o OUTPUT_DIR] [-t THREADS] [--fast | --consensus] [--outname OUT.tsv] [--log LOG_DIR]
```

* `ingroup_alignment.fasta`: Ingroup **codon-aligned** FASTA with **≥2 sequences**
* `reference.fasta`: Either a **single reference** FASTA or an **alignment** used to derive a reference (**same behavior as `groupwise`**)
* `--fast`: Use the **most frequent sequence** from the reference MSA as the outgroup reference
* `--consensus`: Use the **consensus** sequence from the reference MSA as the outgroup reference
  *(If neither flag is used and the reference is an alignment, an **ancestral sequence** is reconstructed via IQ-TREE.)*
* `-o`, `--output_dir`: Output directory (default: current directory)
* `--outname`: Custom output filename (default: `<ingroup_base>_MK.tsv`)
* `-t`, `--threads`: Number of threads for reference inference
* `--log`: Directory to save log file (default: output_dir)

**Output**:
A **TSV file** with header:

```
Ingroup   Outgroup   Pn   Ps   Dn   Ds   NI   alpha   Fisher_p
```

* **Pn**/**Ps**: polymorphic nonsynonymous/synonymous sites in the ingroup
* **Dn**/**Ds**: fixed nonsynonymous/synonymous differences between the **ingroup consensus** and the **reference**
* **NI = (Pn/Ps) / (Dn/Ds)**; **alpha = 1 − NI**
* **Fisher_p**: p-value from Fisher’s exact test on the 2×2 table `[[Pn, Ps], [Dn, Ds]]`
* A **continuity correction** (pseudo-count 0.5) is applied if necessary to avoid division by zero

> **Note**: Using a **consensus or most frequent sequence** as reference is typical in clade-level MK tests when a true outgroup is unavailable. For highly identical clades, Dn/Ds may be small, making MK less stable — consider applying minor allele frequency (MAF) filtering and verifying counts.

---

## Examples

**Pairwise comparison:**

```bash
python two_mode_test.py pairwise example_alignment.fasta -o results/ --format long
```

**Groupwise comparison using consensus reference:**

```bash
python two_mode_test.py groupwise example_alignment.fasta ref_alignment.fasta -o results/ --consensus
```

**MK test using ingroup alignment and consensus outgroup:**

```bash
python two_mode_test.py mk ingroup_aln.fasta outgroup_alignment.fasta -o results/ --consensus
```

**MK test using most frequent sequence as outgroup reference:**

```bash
python two_mode_test.py mk ingroup_aln.fasta outgroup_alignment.fasta -o results/ --fast
```

---

## Output Summary

* **Pairwise**:
  `*_pairwise_dnds.csv` (long/matrix format depending on `--format`)

* **Groupwise**:
  `*_groupwise_dnds.csv` with dN/dS vs reference, and (if applicable) a derived reference FASTA

* **MK test**:
  `*_MK.tsv` with **Pn, Ps, Dn, Ds, NI, alpha, Fisher_p**

---

## Notes

* All modes expect **codon-aligned** FASTA inputs (e.g., produced by TranslatorX or MACSE).
* In **mk** and **groupwise**, when using a reference **alignment**, `--consensus` is typically the most stable; `--fast` selects the most frequent sequence present in the MSA.
* If neither `--consensus` nor `--fast` is used with an alignment reference, the script will run **IQ-TREE** to reconstruct an **ancestral** sequence and use that as the reference.
* For large ingroup sets in MK tests, consider downstream filtering of rare variants to reduce polymorphism inflation.

---
