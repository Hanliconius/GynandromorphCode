#!/usr/bin/env python3
"""Offline safety tests for full-genome Rosina/Yb/Optix joint-calling plan."""
import csv
import importlib.util
import inspect
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/figures/figure3/rosina_whole_genome_bcftools.py"
spec=importlib.util.spec_from_file_location("whole_rosina",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class GenomeWidePlanTests(unittest.TestCase):
    def setUp(self):
        work=tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.base=Path(work.name)
        self.source=self.base/"oct14"
        self.source.mkdir()
        rows=[]
        for i in range(44):
            sample="SRR11243351" if i==14 else "C"+str(i)
            rows.append(dict(sample=sample,specimen="BIO"+str(i),
                 group="cythera" if i<5 else ("vulcanus" if i<10 else "PK"),
                 role="comparative",side="",label="P4003" if i==14 else "P"+str(i)))
        for i in range(4):
            rows.append(dict(sample="F"+str(i),specimen="FBIO"+str(i),
                 group="focal",role="focal",side="Left" if i<2 else "Right",
                 label="T"+str(i)))
        self.rows=rows
        with (self.source/"manifest.tsv").open("w") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        self.fasta=self.base/"Rosina.fna"
        self.fasta.write_text("ACTG\n")
        self.refs=[
            ("NC_140403.1",4_000_000),("NC_140406.1",3_000_000),
            *[("scaffold"+str(i),i*100_000+1) for i in range(1,11)]
        ]
        Path(str(self.fasta)+".fai").write_text("".join(
            name+"\t"+str(length)+"\t0\t0\t0\n" for name,length in self.refs))
        (self.source/"run.json").write_text(json.dumps(
            dict(original_fasta=str(self.fasta),samples=rows,
                 comparative=[r["sample"] for r in rows[:44]],
                 focal=[r["sample"] for r in rows[44:]])))
        (self.source/"jobs.json").write_text(json.dumps(dict(align=m.ALIGN_ARRAY)))
        old=self.base/"old.dedup.bam"
        old.write_bytes(b"BAM_EXAMPLE_NOT_OPENED_BY_PLAN")
        self.audit=self.base/"existing_audit.tsv"
        with self.audit.open("w") as f:
            fields=["sample","candidate","status",
                "source_ref_explicit_match","long_contigs_mapped",
                "picard_metrics_present","has_picard_pg"]
            writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
            writer.writeheader()
            writer.writerow(dict(sample="C0",candidate=str(old),
                status="POTENTIAL_REUSE",source_ref_explicit_match="YES",
                long_contigs_mapped="18",picard_metrics_present="YES",
                has_picard_pg="NO"))
            writer.writerow(dict(sample="C1",candidate=str(old),
                status="POTENTIAL_REUSE",source_ref_explicit_match="NO",
                long_contigs_mapped="18",picard_metrics_present="YES",
                has_picard_pg="NO"))

    def test_every_reference_base_covered_once(self):
        refs=m.reference_contigs(self.fasta)
        win=m.windows_for(refs,width=1_000_000)
        self.assertEqual(len(refs),12)
        for i,r in enumerate(refs):
            pieces=[x for x in win if x["reference_index"]==i]
            self.assertEqual(pieces[0]["start"],1)
            self.assertEqual(pieces[-1]["end"],r["length"])
            self.assertEqual(sum(x["end"]-x["start"]+1 for x in pieces),
                             r["length"])
            for x,y in zip(pieces,pieces[1:]):
                self.assertEqual(x["end"]+1,y["start"])
        self.assertEqual(len({w["task"] for w in win}),len(win))

    def test_one_genuine_prior_bam_approved_without_unverified_candidates(self):
        p=m.build_plan(self.source,self.audit,1_000_000)
        self.assertEqual(len(p["comparative"]),43)
        self.assertEqual(len(p["focal"]),4)
        self.assertEqual(p["verified_reuse"],["C0"])
        self.assertNotIn(14,p["pending_alignment_ids"])
        self.assertEqual(len(p["pending_alignment_ids"]),46)
        self.assertEqual({r["sample"] for r in p["samples"]},
                         {r["sample"] for r in self.rows if r["sample"]!="SRR11243351"})
        self.assertTrue(p["no_hmel25"])

    def test_stale_candidate_audit_is_rejected(self):
        self.audit.write_text("sample\tcandidate\nC0\tfake.bam\n")
        with self.assertRaises(RuntimeError):
            m.build_plan(self.source,self.audit,1_000_000)

    def test_no_individual_caller_or_other_reference(self):
        fn=inspect.getsource(m.call_window)
        self.assertIn('"-b",str(output/"comparative_bams.list")',fn)
        self.assertIn('"mpileup"',fn)
        self.assertIn('"call"',fn)
        self.assertIn('"POS>="',fn)
        self.assertEqual(set(m.STAGES),
            {"validate","call","chrom","combine","extract","focal","publish"})
        self.assertNotIn("haplotype",m.STAGES)

    def test_both_yb_and_optix_handoffs_without_recalls(self):
        text=inspect.getsource(m.export_locus)
        self.assertIn('"Yb"',text)
        self.assertIn('"Optix"',text)
        self.assertIn('"view"',text)
        self.assertNotIn('"mpileup"',text)


if __name__=="__main__":
    unittest.main()
