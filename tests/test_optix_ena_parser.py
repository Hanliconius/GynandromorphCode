#!/usr/bin/env python3
"""Offline contract tests for historical ERS and canonical ENA BioSample handling."""
import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/figures/figure3/optix_joint_candidate_review.py"
SPEC = importlib.util.spec_from_file_location("optix_ena", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def example(**kwargs):
    row = {
        "run_accession": "ERR123456",
        "sample_accession": "SAMEA123456",
        "secondary_sample_accession": "ERS123456",
        "library_strategy": "WGS",
        "library_layout": "PAIRED",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/ERR123456_1.fastq.gz;ftp.sra.ebi.ac.uk/ERR123456_2.fastq.gz",
        "fastq_md5": "a" * 32 + ";" + "b" * 32,
        "fastq_bytes": "123;456",
    }
    row.update(kwargs)
    return row


class TestENAParser(unittest.TestCase):
    def test_primary_and_secondary_sample_ids(self):
        parsed, reason = module.parse_pair(example())
        self.assertEqual(reason, "OK")
        self.assertEqual(parsed[1], "SAMEA123456")
        self.assertEqual(parsed[2], "ERS123456")

    def test_historical_ers_lookups(self):
        record = example()
        old = module.ena_record
        try:
            module.ena_record = lambda _: ([record], "")
            _, rows, error = module.lookup_one(
                ("cythera", "CJ2856", "ERS123456", "test"))
        finally:
            module.ena_record = old
        self.assertFalse(error)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["biosample"], "SAMEA123456")
        self.assertEqual(rows[0]["secondary_sample_accession"], "ERS123456")

    def test_extra_singleton_file(self):
        row = example(
            fastq_ftp=example()["fastq_ftp"] + ";ftp.sra.ebi.ac.uk/ERR123456.fastq.gz",
            fastq_md5=example()["fastq_md5"] + ";" + "c" * 32,
            fastq_bytes="123;456;90")
        parsed, reason = module.parse_pair(row)
        self.assertEqual(reason, "OK")
        self.assertEqual(len(parsed[3]), 2)

    def test_missing_data_is_not_accepted(self):
        parsed, reason = module.parse_pair(example(fastq_md5=None))
        self.assertIsNone(parsed)
        self.assertNotEqual(reason, "OK")

    def test_capture_record_is_not_wgs(self):
        parsed, reason = module.parse_pair(example(library_strategy="OTHER"))
        self.assertIsNone(parsed)
        self.assertEqual(reason, "NOT_WGS")

    def test_false_specimen_link_is_rejected(self):
        old = module.ena_record
        try:
            module.ena_record = lambda _: ([example()], "")
            _, rows, reason = module.lookup_one(
                ("cythera", "incorrect", "ERS987654", "test"))
        finally:
            module.ena_record = old
        self.assertEqual(rows, [])
        self.assertIn("HISTORICAL_SAMPLE_ID_MISMATCH", reason)


if __name__ == "__main__":
    unittest.main()
