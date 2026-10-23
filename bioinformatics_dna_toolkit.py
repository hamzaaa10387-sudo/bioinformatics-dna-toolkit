#!/usr/bin/env python3
"""A small, reproducible DNA FASTA analysis toolkit.

Core analysis is local. Optional NCBI BLAST requires internet access and can take
time. This tool is educational and is not intended for clinical diagnosis.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

VALID_DNA = set("ACGTN")


def analyze_sequence(record: Any, min_orf_aa: int = 30) -> dict[str, Any]:
    """Calculate sequence statistics and find simple forward-strand ORFs.

    ORFs are reported as candidate open reading frames in the input orientation.
    This simple implementation scans all three forward reading frames only;
    it does not establish that an ORF is a real gene.
    """
    sequence = str(record.seq).upper().replace(" ", "").replace("\n", "")
    invalid = sorted(set(sequence) - VALID_DNA)
    if invalid:
        raise ValueError(
            f"Sequence '{record.id}' contains unsupported symbols: {', '.join(invalid)}"
        )
    length = len(sequence)
    counts = Counter(sequence)
    canonical_count = sum(counts[b] for b in "ACGT")
    gc_count = counts["G"] + counts["C"]
    gc_percent = (gc_count / canonical_count * 100) if canonical_count else 0.0
    at_percent = ((counts["A"] + counts["T"]) / canonical_count * 100) if canonical_count else 0.0

    return {
        "id": str(record.id),
        "description": str(record.description),
        "length_bp": length,
        "A": counts["A"],
        "C": counts["C"],
        "G": counts["G"],
        "T": counts["T"],
        "N": counts["N"],
        "ambiguous_bases": counts["N"],
        "canonical_bases": canonical_count,
        "gc_percent": round(gc_percent, 3),
        "at_percent": round(at_percent, 3),
        "orf_candidates": find_orfs(sequence, min_orf_aa=min_orf_aa),
    }


def find_orfs(sequence: str, min_orf_aa: int = 30) -> list[dict[str, Any]]:
    """Find ATG-to-stop ORFs in the three forward reading frames.

    Coordinates are 1-based inclusive. Only complete ORFs with a stop codon
    are returned. Ambiguous codons may translate to X and are not treated as
    stop codons.
    """
    if min_orf_aa < 1:
        raise ValueError("min_orf_aa must be at least 1")
    seq = sequence.upper()
    stops = {"TAA", "TAG", "TGA"}
    results: list[dict[str, Any]] = []
    for frame in range(3):
        starts: list[int] = []
        for pos in range(frame, len(seq) - 2, 3):
            codon = seq[pos:pos + 3]
            if codon == "ATG":
                starts.append(pos)
            elif codon in stops and starts:
                for start in starts:
                    coding_dna = seq[start:pos + 3]
                    peptide = str(Seq(coding_dna).translate(to_stop=True))
                    if len(peptide) >= min_orf_aa:
                        results.append({
                            "start_1based": start + 1,
                            "end_1based": pos + 3,
                            "frame": frame + 1,
                            "length_nt": len(coding_dna),
                            "length_aa": len(peptide),
                            "protein_sequence": peptide,
                        })
                starts = []
    return sorted(results, key=lambda x: (x["start_1based"], x["end_1based"]))


def load_fasta(path: str | Path) -> list[Any]:
    """Read FASTA records and reject empty files or empty sequences."""
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in: {fasta_path}")
    for record in records:
        if len(record.seq) == 0:
            raise ValueError(f"Record '{record.id}' has an empty sequence")
    return records


def analyze_fasta(input_path: str | Path, min_orf_aa: int = 30) -> list[dict[str, Any]]:
    """Analyze all records in a FASTA file."""
    return [analyze_sequence(record, min_orf_aa=min_orf_aa) for record in load_fasta(input_path)]


def write_json(results: list[dict[str, Any]], output_path: str | Path) -> None:
    """Write analysis results as JSON."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)


def write_csv(results: list[dict[str, Any]], output_path: str | Path) -> None:
    """Write summary statistics to CSV (ORFs summarized as a count)."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "id", "description", "length_bp", "A", "C", "G", "T", "N",
        "ambiguous_bases", "canonical_bases", "gc_percent", "at_percent",
        "orf_candidate_count",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            row = {key: result.get(key, "") for key in fields if key != "orf_candidate_count"}
            row["orf_candidate_count"] = len(result["orf_candidates"])
            writer.writerow(row)


def run_blast(input_path: str | Path, output_path: str | Path,
              database: str = "nt", program: str = "blastn",
              max_records: int = 10, email: str | None = None) -> None:
    """Run optional remote NCBI BLAST and save results in XML format.

    Use responsibly: NCBI remote services have usage policies and are not
    intended for high-volume batch jobs. Provide an email address when possible.
    """
    from Bio.Blast import NCBIWWW  # imported only if BLAST is requested

    records = load_fasta(input_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as output:
        for record in records:
            handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=str(record.seq),
                hitlist_size=max_records,
                entrez_query=None,
            )
            output.write(f"<!-- Query: {record.id} -->\n")
            output.write(handle.read())
            handle.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyze DNA FASTA sequences and optionally run NCBI BLAST."
    )
    parser.add_argument("--input", default="data/sample.fasta", help="Input FASTA file")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--min-orf-aa", type=int, default=30,
                        help="Minimum ORF peptide length in amino acids (default: 30)")
    parser.add_argument("--blast", action="store_true",
                        help="Also submit each sequence to remote NCBI BLAST")
    parser.add_argument("--blast-db", default="nt", help="NCBI BLAST database (default: nt)")
    parser.add_argument("--blast-program", default="blastn", help="BLAST program (default: blastn)")
    parser.add_argument("--blast-hits", type=int, default=10, help="Maximum BLAST hits per query")
    args = parser.parse_args(argv)

    try:
        results = analyze_fasta(args.input, min_orf_aa=args.min_orf_aa)
        outdir = Path(args.outdir)
        write_json(results, outdir / "analysis.json")
        write_csv(results, outdir / "summary.csv")
        print(f"Analyzed {len(results)} sequence(s).")
        print(f"JSON report: {outdir / 'analysis.json'}")
        print(f"CSV summary: {outdir / 'summary.csv'}")
        for item in results:
            print(
                f"- {item['id']}: {item['length_bp']} bp; "
                f"GC={item['gc_percent']:.3f}%; "
                f"ambiguous bases={item['ambiguous_bases']}; "
                f"ORF candidates={len(item['orf_candidates'])}"
            )
        if args.blast:
            run_blast(args.input, outdir / "blast_results.xml",
                      database=args.blast_db, program=args.blast_program,
                      max_records=args.blast_hits)
            print(f"BLAST XML: {outdir / 'blast_results.xml'}")
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        if args.blast:
            print(f"BLAST or analysis failed: {exc}", file=sys.stderr)
            return 3
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
