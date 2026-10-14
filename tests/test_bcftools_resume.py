#!/usr/bin/env python3
"""Offline contract checks: restore bcftools caller, exclude P4003, no GATK stages."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/figures/figure3/resume_bcftools_optix_joint.py"
spec=importlib.util.spec_from_file_location("bcftools_resume",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestBcftoolsResume(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)
        originals=[]
        # Precisely 44 comparative specimens and four focal libraries.
        for i in range(44):
            sample="SRR11243351" if i==14 else "SAMPLE"+str(i)
            group=("cythera" if i in (0,1,2,3,4) else
                   "vulcanus" if i in (5,6,7,8,9) else "PK_pale")
            originals.append(dict(sample=sample,specimen="BIO"+str(i),
                                  label="P4003" if i==14 else "label"+str(i),
                                  role="comparative",group=group,side=""))
        for i in range(4):
            originals.append(dict(sample="FOCAL"+str(i),
                specimen="FOCALBIO"+str(i),label="focal"+str(i),
                role="focal",group="focal",
                side="Left" if i<2 else "Right"))
        with (self.path/"manifest.tsv").open("w") as f:
            import csv
            w=csv.DictWriter(f,fieldnames=list(originals[0]),delimiter="\t")
            w.writeheader()
            w.writerows(originals)
        (self.path/"run.json").write_text(json.dumps(
            dict(comparative=[r["sample"] for r in originals[:44]],
                 focal=[r["sample"] for r in originals[44:]])))

    def test_only_p4003_is_excluded(self):
        selected,comparisons,focal,indices=mod.cohort(self.path)
        self.assertEqual((len(selected),len(comparisons),len(focal)),(47,43,4))
        self.assertEqual(indices,list(range(14))+list(range(15,48)))
        self.assertNotIn("SRR11243351",[r["sample"] for r in selected])

    def test_established_caller_parameters(self):
        self.assertEqual(mod.MP_OPTIONS,[
            "-Ou","-a","FORMAT/DP,FORMAT/AD","-q","30","-Q","20"])
        self.assertEqual(mod.CALL_OPTIONS,["-m","-A","-f","GQ","-Oz"])
        self.assertEqual(set(mod.STAGES),
                         {"validate","shard","concat","focal","audit"})
        self.assertNotIn("haplotype",mod.STAGES)
        self.assertNotIn("joint",mod.STAGES)

    def test_full_loci_sharded_without_gaps_or_overlap(self):
        shards=mod.make_shards()
        self.assertEqual(len(shards),56)
        from collections import Counter
        self.assertEqual(Counter(s["locus"] for s in shards),
                         {"Optix":36,"Yb":20})
        for locus,(chrom,left,right) in mod.LOCI.items():
            chunk=[s for s in shards if s["locus"]==locus]
            self.assertEqual(chunk[0]["start"],left)
            self.assertEqual(chunk[-1]["end"],right)
            self.assertTrue(all(s["chrom"]==chrom for s in chunk))
            for a,b in zip(chunk,chunk[1:]):
                self.assertEqual(a["end"]+1,b["start"])
            self.assertEqual(sum(s["end"]-s["start"]+1 for s in chunk),
                             right-left+1)

    def test_array_shards_not_per_sample(self):
        # Every shard calls all comparative BAMs in the single cohort-wide list.
        import inspect
        source=inspect.getsource(mod.shard)
        self.assertIn('"-b",cfg["bam_list"]',source)
        self.assertIn('"mpileup"',source)
        self.assertIn('"call"',source)

    def test_new_pipeline_has_explicit_per_locus_concat(self):
        import inspect
        source=inspect.getsource(mod.concat)
        self.assertIn('"concat"',source)
        self.assertIn('require(str(tmp)+".tbi")',source)

    def test_does_not_silently_drop_other_samples(self):
        data=(self.path/"manifest.tsv").read_text().replace("SAMPLE13","SRR11243351")
        (self.path/"manifest.tsv").write_text(data)
        with self.assertRaises(RuntimeError):
            mod.cohort(self.path)


if __name__=="__main__":
    unittest.main()
