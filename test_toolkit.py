import tempfile
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bioinformatics_dna_toolkit import analyze_sequence, find_orfs, load_fasta, write_csv, write_json


class TestSequenceAnalysis(unittest.TestCase):
    def test_gc_content(self):
        record = SeqRecord(Seq("AACCGGTT"), id="test")
        result = analyze_sequence(record, min_orf_aa=1)
        self.assertEqual(result["length_bp"], 8)
        self.assertEqual(result["gc_percent"], 50.0)
        self.assertEqual(result["A"], 2)
        self.assertEqual(result["N"], 0)

    def test_ambiguous_bases_excluded_from_gc_denominator(self):
        record = SeqRecord(Seq("GGCCNN"), id="test")
        result = analyze_sequence(record, min_orf_aa=1)
        self.assertEqual(result["gc_percent"], 100.0)
        self.assertEqual(result["ambiguous_bases"], 2)

    def test_invalid_symbol_rejected(self):
        record = SeqRecord(Seq("ACGTX"), id="bad")
        with self.assertRaises(ValueError):
            analyze_sequence(record)

    def test_orf_coordinates(self):
        # ATG + 30 alanine codons + TAA gives a complete 31-aa ORF.
        sequence = "ATG" + "GCT" * 30 + "TAA"
        result = find_orfs(sequence, min_orf_aa=30)
        self.assertTrue(any(x["start_1based"] == 1 and x["end_1based"] == len(sequence)
                            for x in result))

    def test_outputs_written(self):
        record = SeqRecord(Seq("AACCGGTT"), id="test")
        result = [analyze_sequence(record, min_orf_aa=1)]
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / "nested" / "out.csv"
            json_path = Path(tmp) / "nested" / "out.json"
            write_csv(result, csv_path)
            write_json(result, json_path)
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())


if __name__ == "__main__":
    unittest.main()
