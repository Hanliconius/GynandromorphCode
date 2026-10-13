# Post hoc dark and pale piano-key controls

Start raw-read acquisition immediately with `extend_PK_optix_Yb.py --download-new`.
This submits one download task per PK accession absent from the original optix
VCF, without needing either reference genome. The full `--submit --optix-run RUN`
workflow automatically discovers active download array elements from these
download runs and waits on them rather than duplicating downloads. Completed
downloads are revalidated by the download worker. Failed elements are retried.
Supply `--optix-reference FASTA` when the matching indexed whole genome cannot
be discovered automatically; its frozen optix SNP REF bases must agree.

This extension uses the 21 individuals explicitly labelled Piano Keys, dark
or Piano Keys, pale in the supplied sample table: 12 dark, 9 pale. It excludes
cydno, all cydno crosses and backcrosses, and ivory-only individuals. All 21
are B-positive; six are Dennis-positive and four are Ray-positive. These are
recorded phenotype annotations, not genotype calls. Inferred pattern genotype
and ivory annotations are stored separately in the manifest.

The existing comparative SNP scores, allele orientation and focal reciprocal
mixture calls are frozen. PK individuals are appended after ascertainment.
Optix displays retain the BA >=0.80 and >=0.90 labelled-interval selections,
remove missing focal pairs and retain conflicting focal genotype states.

The workflow inspects the existing optix VCF and reuses its available PK calls.
Missing individuals require a full reference genome with Hmr18 and exact REF
agreement at every frozen SNP. An unambiguous indexed genome is discovered
locally, or supplied with --optix-reference. No coordinate transfer is guessed.
New optix calls use known-allele-constrained diploid bcftools calling with
MAPQ >=30, BQ >=20, DP >=2 and GQ >=20. Their call source is recorded separately
from the existing VCF calls.

Yb uses both prepared whole genomes (Hmel2.5 and HmRos). Existing or queued
five-PK alignments are reused by individual task dependencies. New data are
aligned competitively against the full genome; duplicates are marked before
extraction. Only the whole Yb contig and optix region BAMs are retained from
new alignments; temporary full-genome BAMs are removed after successful local
extraction. Whole-contig Yb data are retained to normalize coverage by each
sample's own contig mean, including zero-depth positions. MAPQ >=30, duplicate
exclusion and paired overlap suppression apply. Local Yb VCFs are retained
for subsequent allele-evidence analyses; no Yb phenotype is inferred from
ivory, Band, Dennis or Ray states.

Download tasks validate ENA checksums, byte counts and BioSample identities.
The existing five-PK download pipeline is not duplicated. New task arrays are
unthrottled; the scheduler assigns resources. New profile tasks use aftercorr
dependencies on their matching alignment tasks. Reused profile tasks depend
only on the relevant original alignment array elements. A final nano report
depends on all profile arrays, and assembles individual optix genoplots and
Yb coverage pages into PK_dark_pale_optix_Yb_posthoc.pdf.

Executable workers, report dependencies, sample annotations and the existing
optix values are frozen per run. A complete.json is written only after the
report succeeds. Individual stage completion markers support safe task reruns.
