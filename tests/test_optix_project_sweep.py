#!/usr/bin/env python3
"""Offline tests for independent cythera discovery and race-label audits."""
import importlib
import sys
import unittest
from pathlib import Path

FIG = Path(__file__).resolve().parents[1] / "scripts/figures/figure3"
sys.path.insert(0, str(FIG))
sweep = importlib.import_module("optix_joint_project_sweep")
prior = importlib.import_module("optix_joint_candidate_review")


def record(**changes):
    row = {
        "run_accession": "SRR5000001",
        "sample_accession": "SAMN5000001",
        "secondary_sample_accession": "SRS5000001",
        "library_strategy": "WGS",
        "library_layout": "PAIRED",
        "scientific_name": "Heliconius melpomene cythera",
        "sample_title": "cythera individual new",
        "sample_alias": "",
        "sample_description": "",
        "sub_species": "",
        "tax_id": "123",
        "study_accession": "PRJNA485723",
        "fastq_ftp": "archive/SRR5000001_1.fastq.gz;archive/SRR5000001_2.fastq.gz",
        "fastq_md5": "a"*32 + ";" + "b"*32,
        "fastq_bytes": "500;500",
    }
    row.update(changes)
    return row


class TestCytheraProjectSweep(unittest.TestCase):
    def test_confirmed_taxon(self):
        self.assertEqual(sweep.race_evidence(record())[:2],
                         ("cythera", "TAXON_CONFIRMED"))

    def test_no_cross_race_assignment(self):
        r = record(scientific_name="Heliconius melpomene vulcanus",
                   sample_title="vulcanus individual")
        self.assertEqual(sweep.race_evidence(r)[1], "NOT_CONFIRMED")

    def test_title_only_requires_review(self):
        r = record(scientific_name="Heliconius melpomene")
        self.assertEqual(sweep.race_evidence(r)[1], "MANUAL_REVIEW")

    def test_prior_cythera_not_counted_as_new(self):
        r = sweep.annotate(record(), "project:PRJNA485723",
                           {"SAMN5000001"})
        self.assertEqual(r["already_in_previous"], "yes")
        self.assertEqual(r["wgs_pair"], "yes")

    def test_capture_record_rejected(self):
        r = sweep.annotate(record(library_strategy="OTHER"),
                           "project:PRJNA485723", set())
        self.assertEqual(r["wgs_pair"], "no")
        self.assertEqual(r["reason"], "NOT_WGS")

    def test_vulcanus_mismatch_flagged(self):
        r = sweep.annotate(record(), "audit:vulcanus:SRR5000001",
                           {"SAMN5000001"})
        old = [{"subspecies": "vulcanus", "biosample": "SAMN5000001",
                "run_accession": "SRR5000001"}]
        audit = sweep.audit_previous(old, [r])
        self.assertEqual(audit[0]["status"], "TAXON_CONFLICT")


if __name__ == "__main__":
    unittest.main()
