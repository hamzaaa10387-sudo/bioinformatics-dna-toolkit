# Bioinformatics DNA Toolkit

A small, reproducible Python project for DNA sequence quality checks and first-pass sequence analysis. It reads FASTA files, calculates nucleotide statistics, identifies simple forward-strand open reading frame (ORF) candidates, and exports JSON/CSV reports. An optional flag submits sequences to NCBI BLAST.

> **Scope:** This is an educational portfolio project. ORF predictions and BLAST similarity are exploratory evidence, not proof of gene function and not a clinical diagnostic.

## Features

- Parse one or more records from a FASTA file using Biopython.
- Validate DNA symbols (`A`, `C`, `G`, `T`, `N`) and report clear errors.
- Calculate sequence length, A/C/G/T/N counts, GC%, AT%, and ambiguous-base count.
- Find complete ATG-to-stop ORF candidates in the three forward reading frames.
- Export detailed JSON and summary CSV reports.
- Optional remote NCBI BLAST search, saved as XML.
- Unit tests for calculations, validation, ORFs, and report writing.

## Repository structure

```text
bioinformatics-dna-toolkit/
├── bioinformatics_dna_toolkit.py
├── data/
│   └── sample.fasta
├── results/
│   └── .gitkeep
├── tests/
│   └── test_toolkit.py
├── .gitignore
├── README.md
└── requirements.txt
```

## Requirements

- Python 3.10 or newer recommended
- pip
- Internet connection only for optional NCBI BLAST

## Setup

### Windows (PowerShell)

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

If PowerShell blocks activation, use Command Prompt:

```bat
py -m venv .venv
.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### macOS/Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run the analysis

From the repository root:

```bash
python bioinformatics_dna_toolkit.py --input data/sample.fasta --outdir results
```

The program writes:

- `results/analysis.json`: detailed record-level statistics and ORF candidates
- `results/summary.csv`: spreadsheet-friendly summary

Change the minimum candidate ORF length (amino acids) with:

```bash
python bioinformatics_dna_toolkit.py --input data/sample.fasta --outdir results --min-orf-aa 10
```

## Optional: run NCBI BLAST

```bash
python bioinformatics_dna_toolkit.py --input data/sample.fasta --outdir results --blast
```

Optional settings:

```bash
python bioinformatics_dna_toolkit.py --input data/sample.fasta --outdir results --blast --blast-db nt --blast-program blastn --blast-hits 10
```

BLAST requires internet access and depends on NCBI service availability. Remote searches may take time; do not submit high-volume jobs. Review current NCBI BLAST usage guidance before repeated or large-scale queries. The generated XML is a raw result file; inspect it using Biopython's BLAST XML parser or the NCBI web interface before reporting any biological interpretation.

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Interpretation notes

- **GC%** is calculated as `(G + C) / (A + C + G + T) × 100`. Ambiguous `N` bases are excluded from the denominator.
- **AT%** uses the same canonical-base denominator.
- If a sequence has no canonical bases, GC% and AT% are reported as 0.
- **ORFs:** this implementation scans only the three forward frames and reports complete ATG-to-stop candidates meeting the chosen amino-acid threshold. It does not scan the reverse complement, predict promoters, or confirm genes.
- **BLAST:** sequence similarity depends on the database, parameters, query quality, and significance statistics. A high-scoring match alone does not prove identical function.

## Suggested next improvements

1. Add reverse-complement ORF scanning with strand-aware coordinates.
2. Add a `pandas` report and GC-content visualization.
3. Parse BLAST XML into a ranked tabular report with accession, identity, query coverage, E-value, and bit score.
4. Add GitHub Actions to run unit tests on every push.
5. Add a small, cited public reference dataset and record its source/access date.

## Reproducibility and attribution

- Python and Biopython versions should be recorded when reporting results.
- `data/sample.fasta` contains short synthetic teaching sequences and is not a biological reference dataset.
- Cite the source and version/date of any real sequence data added later.
