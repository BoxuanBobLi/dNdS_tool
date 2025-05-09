
---

# Sequence Processing Script

This Python script is a versatile tool for handling genomic FASTA files. It supports several modes for processing individual sequences, large collections, or groups of sequences in a readable format. The output is typically stored in a specified directory, optionally compressed into a `.tar.gz` archive.

## Features

* **Large mode**: Processes all FASTA files in a directory and saves each sequence as a separate file.
* **Single mode**: Processes a single FASTA file and saves each record as a separate file.
* **Groupwise mode**: Groups FASTA records by a specified attribute (e.g., prefix) and writes them accordingly.
* **Readable mode**: Converts a FASTA file into a more human-readable format (e.g., TSV with ID and sequence).

---

## Requirements

* Python 3.6+
* Biopython (`pip install biopython`)

---

## Usage

```bash
python script.py --input <input_path> --output <output_directory> --mode <mode> [--name <name>] [--compress]
```

### Arguments

| Argument     | Description                                                      |
| ------------ | ---------------------------------------------------------------- |
| `--input`    | Path to the input FASTA file or directory                        |
| `--output`   | Path to the output directory                                     |
| `--mode`     | Mode of operation: `large`, `single`, `groupwise`, or `readable` |
| `--name`     | Optional base name for the output files (used in single mode)    |
| `--compress` | Optional flag to compress output directory into `.tar.gz` format |

---

## Examples

### Large Mode

Process all FASTA files in a directory:

```bash
python script.py --input data/ --output output/ --mode large
```

### Single Mode

Process one FASTA file and save records individually:

```bash
python script.py --input genome.fasta --output output/ --mode single --name genome_seq
```

### Groupwise Mode

Process a FASTA file and group records by shared prefix or other logic:

```bash
python script.py --input grouped_sequences.fasta --output output/ --mode groupwise
```

### Readable Mode

Convert a FASTA file into a tab-separated readable format:

```bash
python script.py --input genome.fasta --output output/ --mode readable
```

---

## Output

* FASTA records are written as individual `.fasta` files named by their record ID.
* If `--compress` is used, output is archived as a `.tar.gz` file.

---
