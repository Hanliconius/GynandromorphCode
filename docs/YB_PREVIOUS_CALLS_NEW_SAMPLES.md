# Additional samples displayed on the previous Yb calls

2026-10-16 correction to the restoration: preserve the old calling method
and show new samples without using them for ascertainment.

`scripts/yb/yb_previous_calls_new_samples.py --submit` takes the completed
restoration log, the original interval count run, and the new shared Rosina
run. It freezes the actual restored genotype tables, original count tables,
saved fit, archived caller, previous PK calls and race-frequency tables.
It also copies the original fragment counter from the pinned restoration.

The original compatible SNP positions and allele orientation are read from
the restored `*_genotype_pairs.tsv.gz` tables. Left/Right dosages are copied
directly, including missing values. There is no new focal recounting, genotype
inference, parameter fitting, balanced-accuracy filtering, or SNP selection.
Both original control frequency rows are reused from the original counts.
The worker asserts exact equality of plotted Left/Right values to the frozen
tables and verifies input hashes.

Additional comparative samples come from the retained shared-run manifest,
not filesystem wildcard guesses. Preflight checks BAM indexes, read-group
sample identities and the entire reference sequence dictionary against the
previous Rosina FASTA. Cythera and vulcanus must be present. These BAMs are
Rosina-only; they are never used at Hmel2.5 coordinates. Both Hmel2.5 panels
retain the original display cohort and explicitly say so.

Existing five-PK dosage tables and earlier other-race frequency values are
reused. Remaining PKs and other comparative samples are counted at exactly
the historical SNPs, using the original fragment counter: MAPQ30, baseQ30,
duplicate/secondary/supplementary/QC-fail exclusion, overlapping mate agreement
counted once per fragment per position. Additional races are mean measured
read-frequency strips, with each individual's selected-biallelic depth>=8;
as before, a race cell can use one or more callable individuals. Per-site
individual counts and values are written out. The original rosina/melpomene
control rows stay unchanged; extra rosina/melpomene individuals get separately
labelled display strips.

New PKs use the historical unmixed diploid three-state beta-binomial model
with the original saved error, bias and dispersion and posterior>=0.95,
minimum selected-biallelic depth2. Original five PK dosage values are kept
unchanged. All retained shared-manifest PKs must be included or the job fails.
PK labels show retained specimen IDs, dark/pale group, and user-confirmed
dorsal+/ventral+ phenotype. Phenotype labels do not force homozygous calls.
Cythera is dorsal−/ventral−; vulcanus dorsal−/ventral+, per user confirmation.

The PDF contains four panels (d/v on each reference), the same palette and
locus connections. Added races and PKs appear on the Rosina panels. Exact
plot matrices, frozen focal call copies, new-sample count tables, new PK
posterior calls, membership and an audit are saved in the new output folder.
The original restored outputs are never overwritten. No new SNP recall,
downloads, or alignments are submitted. One nano job (30min, 2 CPUs, 12G)
handles the small local BAM regions and plotting on Pegasus.
