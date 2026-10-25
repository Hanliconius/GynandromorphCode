# Reviewed Figure 3 package

The input authority is `data/figure3_review/final_specimens.tsv`: a 550-row snapshot of the author-reviewed Google Sheet taken on 22 October 2026. Later online edits do not alter a running job. Every run freezes the input bytes and records SHA-256 hashes.

Run on Pegasus:

```bash
git -C /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo pull --quiet --ff-only origin main &&
micromamba run -n gynandrocalling python \
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/figures/figure3/build_Figure3_final_package.py \
--submit --overwrite --yb-cohort lower-panel
```

Two nano array tasks independently calculate Dennis and Yb(v) statistics. A dependent nano job creates the PDFs and share archive. All jobs request 30 minutes and 8 GB. No download, alignment, SNP calling, or focal mixture fitting is performed.

The stable output is `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Figure3_for_Arnaud`. A successful run writes `complete.json` and `Figure3_for_Arnaud.zip` only after all eight PDFs and the archive succeed. Failed reruns overwrite generated files in this same owned directory, and active jobs are protected.

## Figures

1. `01_Figure3_main_BA_dots.pdf`: blue/yellow main panels, red BA >=0.90 dots on the annotations, unchanged corrected Ho, retained nucleotide focus panels.
2. `02_Supplement_full_BA.pdf`: same layout with full BA scatter plots.
3. `03_Dennis_all_individuals.pdf`: all reviewed specimens with existing Optix focus calls, excluding global exclusions. Display-only specimens are labelled BA=no.
4. `04_Dennis_race_frequencies.pdf`: final reviewed Dennis cohort pooled by the exact population and phenotype fields.
5. `05_Ybv_all_individuals.pdf`: all reviewed specimens with existing Yb focus calls, including available references excluded from BA, labelled BA=no.
6. `06_Ybv_race_frequencies.pdf`: final reviewed Yb comparison cohort, pooled by geography, population label and ventral phenotype, merging the cythera ssp/spp spelling variants.
7. `07_BA_selection_individuals.pdf`: specimens actually contributing to the main BA, separately for Dennis and restricted Yb(v).
8. `08_Phenotype_frequencies.pdf`: phenotype-level frequencies for the main cohorts.

The individual plots are paginated and labels stay on one line. Supplement palette: pale blue `#D7EDF5`, heterozygote gold `#ECCB63`, deep purple `#662A83`; grey missing. These supplements use raw ALT dosage. Frequency midpoints do not imply that every individual is heterozygous.

## Statistical decisions

- Use the reviewed `D`, `Yb_ventral`, population names, geography and inclusion flags, not the historical `*_used_Fig3*` columns.
- Deny numata and pardalinus by names, known accessions, or the global exclusion flag. They never appear in genotype plots or contribute to statistics; their metadata remain in the audit snapshot.
- One reviewed specimen contributes once. Multiple qualified run calls are united if concordant; conflicts are missing and counted in audit tables. No run is selected for its phenotype agreement.
- Dennis uses the carrier rule over the final reviewed Dennis cohort. Yb(v) uses the homozygote rule within the exact frozen eight-population cohort (63 positive, 107 negative). Allele orientation is selected within that restricted cohort. The input membership and 27-marker reference scores come from the completed `Ybv_BA_rule_comparison/lower_panel` diagnostic.
- Both require >=80% called within each tested phenotype class. Qualified Optix calls require DP>=2 and reported GQ>=20. Yb uses the unchanged GQ/PL quality engine. Existing legacy/new calling histories are retained; identical calling provenance is not claimed.
- Keep the established 24 Dennis and 27 Yb focus SNP identities. Recalculate their BA, frequencies and class denominators without reselection. Export whether each still passes BA 0.90. Annotation dots show all newly passing SNPs across the full region, not only the lower focus markers.
- Main Dennis consensus and example rows use the final Dennis cohort. Main Yb(v) consensus, full-span BA and the eight population rows use exactly the same 170 specimens. Membership is frozen from the completed lower-panel diagnostic, not selected again by group size or genotype. With elevatus replacing linaresi: positive examples are thelxinoe, cythera, besckei and vulcanus; negative examples are malleti, timareta nov. ECU, heurippa and elevatus. Main colours retain the previous allele polarity: Dennis blue is the original Dennis-associated allele, yellow the other; Yb blue is ALT, yellow REF. Diploid genotypes use one row and unphased top/bottom allele boxes.
- Keep the existing corrected L/R calls and Ho estimator: 10-kb windows, 2.5-kb slide, common callable denominator, posterior 0.95 and 0.99.
- Set explicit Ho y-ticks every 0.05. Window values are plotted at their midpoints, so the first point is approximately 1.065 Mb for Optix and 1.635 Mb for Yb, while axis origins remain 1.060 and 1.630 Mb.
- Keep fixed spans: Optix NC_140406.1:1060000-1240000, Yb NC_140403.1:1630000-1692000. Verified cortex and partial ivory annotations are preserved.

The package includes the metadata snapshot, BA tables, per-site frequency denominators, specimen membership, run-conflict audit, current focus-marker scores and consensus frequencies. The wide genotype matrices and frozen worker code remain in the run's `data` and `code` subdirectories but are omitted from the share ZIP.

## Coordinate and label checks

`tables/Ybv_BA90_coordinate_audit.tsv` records every Yb BA >=0.90 dot with its exact REF/ALT identity, whether it belongs to the 27 displayed markers, its offset from the nearest displayed marker, and whether it lies to the right of the Ybv interval or first mapped ivory exon. Dots retain exact genomic positions and are never snapped to a nearby focus marker. The interval is 1,658,079–1,659,072 and the first mapped ivory exon is 1,659,008–1,659,348.

`tables/population_label_normalizations.tsv` records the cythera merge. The raw reviewed snapshot is preserved; its one `Heliconius melpomene ssp. cythera` population label is pooled with the other 13 cythera specimens. No inclusion flags, phenotypes, geography or calls change. Per-panel JSON records x limits, Ho tick spacing, first window midpoint, lower-panel labels and class counts.

## Restricted-cohort remake

The current default is `--yb-cohort lower-panel`. It requires the completed restricted BA diagnostic and freezes its specimen IDs, recorded population rows and 27-marker homozygote BA reference table before submission. The main SNP rows must cover exactly the main BA membership with no duplicates, additions or omissions. Cythera is one row, and elevatus stays in the negative rows; linaresi is not substituted back in. Consensus totals are 63 v+ and 107 v− for this snapshot.

The full 62-kb BA scan is recalculated from existing calls for these same specimens. At retained markers the exact diagnostic call cache is used; shared called genotypes must agree with the broad cache, and missingness differences or appended retained-marker records are audited. Each retained-marker BA must reproduce the completed restricted diagnostic before plotting. The main rule remains homozygote; carrier results remain in the sensitivity diagnostic. Both main red-dot and supplemental full-BA PDFs are rebuilt in place, together with the package supplements and ZIP. Existing calls, focal mixture estimates, Ho, fixed spans, Dennis cohort and marker identities are not changed.

`--yb-cohort west` retains the earlier western-only sensitivity workflow, but is not the current main figure choice.
