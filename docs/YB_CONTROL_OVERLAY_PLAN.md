# Yb control overlay and annotation audit

Sources:
- Livraghi et al. 2021 eLife 68549 Figure 8 source data 1:
  https://cdn.elifesciences.org/articles/68549/elife-68549-fig8-data1-v1.docx
- Figure 8 source data 2:
  https://cdn.elifesciences.org/articles/68549/elife-68549-fig8-data2-v1.docx
- Krista-Pipho/gynandromorph_2024 config.yaml candidate controls.

The source table has 78 rows, versus 79 individuals reported in the methods.
Its melpomene/timareta subset has 44 rows, with ERR2298207 assigned twice to
opposite phenotype groups (rosina and melpomene). Both rows are retained but
excluded pending metadata resolution. The four Krista candidates remain
phenotype-unresolved. ERS sample accessions need run resolution, not treating
them as downloadable FASTQ run IDs.

First run audit_Yb_control_candidates.py on Pegasus. It searches the existing
Joe scratch area for sequence files and small manifest/sample tables, resolving
accession/run/library aliases before judging local availability. This prevents
the narrow filename-search mistake from the first rebuild. File candidates
are not certified BAM reference matches until headers and provenance are checked.

Next resolve source-table conflicts and missing run IDs, verify Krista's figure
sample identities, and choose available confirmed yellow-bar-present/absent
controls. Download only genuinely missing required reads. Align controls to
both Hmel2.5 and HmRos with matching processing.

The five published Sanger CRE-neighborhood sequences are included verbatim as
FASTA with IUPAC ambiguity retained. Their alignment anchors the region on both
references; their full alignment span is not the exact CRISPR-tested CRE.
Recover the paper's tested ATAC interval/guide coordinates before labelling
an exact tested CRE. Current Yb(d)/Yb(v) annotations alone do not establish that.

For comparable overlays, calculate 50-bp MAPQ>=30 window coverage and divide
each individual's local values by their full target scaffold/chromosome mean.
The paper describes median window depths; implement/review that detail rather
than silently replacing it with the existing unfiltered bamCoverage means.
Keep raw coverage as a companion. Plot independent T2/T3 traces plus Left/Right
summaries, and confirmed Yb-present/absent controls with individual counts.
Use native reference coordinates and annotation provenance on both references.
