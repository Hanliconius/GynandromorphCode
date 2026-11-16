# NDJ chr15: targeted junction follow-up

First pass: run ec8e51556aa7, chr15 [2050000,2070000), Rosina GCF_056956725.1.
No feature survived the exploratory multiple-testing adjustment. The AT-only tract
screen is a lead, not an established association. At 2056847 the largest exact clip
cluster in the bracket had 10 right-side fragments and one left-side fragment.
The nearby inverted-repeat screen was selected post hoc and does not establish causality.

## Run

In the usual gynandrocalling environment on Pegasus:

```bash
python scripts/ndj_breakpoint/followup_junctions.py --submit
```

Default input is the completed ec8e51556aa7 cluster run; override --source if moved.
One nano job, one CPU, 4 GB, 30 minutes. Exact source BAM stat signatures, indexes
and contig lengths are rechecked; no download, realignment or new variant calls.
The job freezes its code and writes per-sample reuse decisions before submission.
Completed outputs are checksum-checked for reuse. Inspect incomplete jobs before
retrying in a fresh --output directory. Output resources are registered in the
existing project ledger under a file lock with a pre-update snapshot.

## What is compared

All exact clip sites in the bracket with >=3 fragments on either side, plus the
20 strongest such sites in the 100-kb flanks from the original scan. These are
selected comparison sites, not independent statistical controls. Sites within
5 bp can share supporting reads and are deliberately retained as separate rows.

For each site and sample, denominator = unique RG+QNAME fragments whose primary
alignment intersects the +/-5-bp anchor neighbourhood, MAPQ >=30. The numerator
has a >=20-bp soft clip at the selected edge within 5 bp. A quality subset requires
mean clipped-base quality >=20. Fractions are descriptive, not independent allele
fractions or proof of mosaicism. Different template start/end/reverse geometries
are also counted for each partner group; this is not UMI-based molecule validation.

SA CIGARs are normalized to original-read query orientation, including hard clips.
Retain SA alignments with MAPQ >=30 that cover >=20 bp and >=half of the clipped
query interval. Report exact partner interbase coordinate, both orientations,
query overlap/gap and supporting fragment count. Raw SA positions are converted
from 1-based to 0-based. Negative query gaps must NOT be labelled microhomology.
The script does not remap clips without SA tags, assemble an alternate junction,
or infer a repair mechanism. These remain limitations of a negative screen.

## Review and stopping rule

Inspect REPORT.md, candidate_summary.tsv and SA_partner_groups.tsv. A descriptive
flag requires >=3 fragments with >=2 distinct primary anchor geometries pointing
to the same exact partner and orientation. It is not a validated SV threshold.
Also inspect nearby partner endpoints: alignment ambiguity can distribute support
over adjacent coordinates. Consistent local partners can be ordinary small indels.

If no coherent partner evidence emerges, report **nothing obvious in the tested
sequence features and short-read junction screen**. Do not keep relaxing motif,
distance or support thresholds to obtain a positive result. Lack of evidence can
also reflect unresolved repetitive sequence or short-read limitations.

If a coherent group does emerge, inspect primary and supplementary alignments and
consider local assembly/read remapping for that specific junction. The existence
of a motif nearby is secondary to this validation. Existing homology_search.py
remains an optional reference-homology search; it does not resolve the specimen's
junction, and is not necessary to force a mechanistic explanation.

Only small summary tables, report, plan and completion marker need to be pushed
for review. No raw sequences or read names appear in these summary outputs.

## Validation status

The development executor was unavailable for this follow-up. Python execution
has NOT been tested here. Included CIGAR/orientation/filter assertions run before
submission and again in the worker; they cover forward/reverse partners, hard
clips, non-clip SA alignments and low-MAPQ exclusion. Actual BAM behavior remains
to be checked in the Pegasus run.
