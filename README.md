
---

# Sequence Processing Script

This Python script processes FASTA files in different ways depending on the selected mode. It supports batch processing of large directories, single-file extraction, groupwise organization, and converting FASTA to a more readable format.

## Features

* **large**: Processes all FASTA files in a given directory and writes each sequence to a separate file.
* **single**: Extracts all records from a single FASTA file, writing each to a separate file using the record ID as the filename.
* **groupwise**: Groups sequences in a single FASTA file by their prefix and writes grouped files to the output directory.
* **readable**: Converts a FASTA file into a tab-separated text format (record ID and sequence).

---

## Requirements

* Python 3.6+
* Biopython (`pip install biopython`)

---

## Usage

```bash
python script.py --input <input_path> --output <output_directory> --mode <mode>
```

### Arguments

| Argument   | Description                                                    |
| ---------- | -------------------------------------------------------------- |
| `--input`  | Path to a FASTA file or directory (depending on mode)          |
| `--output` | Path to the output directory                                   |
| `--mode`   | Processing mode: `large`, `single`, `groupwise`, or `readable` |

---

## Examples

### Process all FASTA files in a directory (`large` mode)

```bash
python script.py --input /path/to/folder --output /path/to/output --mode large
```

### Process a single FASTA file (`single` mode)

```bash
python script.py --input file.fasta --output output_dir --mode single
```

### Group sequences by prefix (`groupwise` mode)

```bash
python script.py --input file.fasta --output output_dir --mode groupwise
```

### Convert FASTA to readable TSV format (`readable` mode)

```bash
python script.py --input file.fasta --output output_dir --mode readable
```

---

## Output

* Output directory will contain `.fasta` or `.txt` files depending on the mode.
* Filenames are based on record IDs or group prefixes.

---


