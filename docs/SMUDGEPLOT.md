# Smudgeplot workflow and calibration

## Purpose

Smudgeplot provides an alignment-free view of covariant k-mer pairs. Here it is used to compare:

- TN Left;
- TN Right;
- T2 cross-side merge;
- T3 cross-side merge;
- all four focal libraries;
- an empirical diploid control;
- a pseudo-tetraploid benchmark made by merging two diploid individuals.

Smudgeplot is treated as a diagnostic of k-mer-pair structure, not as a single-button ploidy caller.

## Upstream baseline

The original workflow inherited from Krista Pipho's gynandromorph_2024 repository used:

FastK -v -t4 -k22 -M100 -T60 <FASTQ inputs> -NFastK_Table

smudgeplot hetmers -L 12 -t 60 -o kmerpairs --verbose FastK_Table

smudgeplot all -o output kmerpairs.smu

Recent reruns used Smudgeplot v0.5.3.

See UPSTREAM_KRISTA.md for the exact upstream commit.

## Arrangement set

Historical k-mer-pair inputs include:

analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu

analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu

analysis/smudgeplot_arrangements/TN_T2_merge/kmerpairs.smu

analysis/smudgeplot_arrangements/TN_T3_merge/kmerpairs.smu

analysis/smudgeplot_arrangements/ERR1051012_krista/kmerpairs_L4.smu

analysis/smudgeplot_arrangements/pseudoTetraploid_ERR1307111_ERR1307112/kmerpairs.smu

analysis/smudgeplot_4legs_combined/kmerpairs.smu

These are provenance paths, not portable script requirements.

## Control calibration

### Earlier controls

ERR2298208 and several controls chosen from mapped BAM depth performed poorly as Smudgeplot single-sample controls.

This established an important point:

Mapped BAM depth is not sufficient to predict useful heterozygous k-mer-pair structure.

### Krista control set

The older upstream config listed:

- ERR1050907
- ERR1050958
- ERR1050994
- ERR1051012

ERR1051012 became the useful empirical diploid control.

### Coverage-harmonic problem

At L=4, ERR1051012 retained a large coherent balanced cloud. Automatic coverage selection could choose a half-coverage harmonic, causing the same balanced cloud to be labelled AABB instead of AB.

Forcing approximately 1n=8.01 on the same L=4 k-mer-pair set assigned about 53.17 million in-smudge k-mers to AB, with an error fraction near 12.4%.

Therefore:

Changing the assumed 1n coverage can relabel an existing cloud without changing the underlying k-mer-pair data.

An apparent AABB call must not be interpreted as copy-number evidence until the coverage harmonic is checked.

## Lower-cutoff sweep

A later systematic sweep tested L values including 6, 8, and 12.

For ERR1051012:

- L=4 was diagnostically useful but admitted additional low-count structure;
- L=6 gave the cleanest final diploid control and automatic AB classification;
- L=8 became too restrictive and lost the useful smudge signal.

The current preferred diploid control is therefore ERR1051012 at L=6.

## Why one universal L is inappropriate

A strongly imbalanced pair can have high total coverage while its minor arm lies below the lower cutoff.

This matters for candidate 1:3-like structure.

At L=12, exploratory diagnostics showed important imbalanced components near:

- pseudo-tetraploid: 12:36;
- T2 merge: 12:29;
- T3 merge: 12:30.

Thus L=12 preferentially censors the low-count arm of some imbalanced pairs.

The lower-cutoff sweep was designed to test this bias explicitly.

## Current preferred comparison

| dataset | role | preferred L |
|---|---|---:|
| ERR1051012 | empirical diploid | 6 |
| ERR1307111 + ERR1307112 | pseudo-tetraploid benchmark | 6 |
| TN Left | focal side | 12 |
| TN Right | focal side | 12 |
| TN T2 merge | cross-side merge | 6 |
| TN T3 merge | cross-side merge | 6 |
| all four TN legs | combined focal sample | 12 |

The final lower-L downloadable set was collected historically under:

analysis/smudgeplot_lowerL_20261002/L6_download/

## Original L=12 TN reference

### TN Left

- total pairs: 53,496,929
- error fraction: about 9.0%
- inferred 1n: 20
- AB: 44,600,980
- AABB: 2,356,462
- tiny higher-order class: 2,084

### TN Right

- total pairs: 22,222,743
- error fraction: about 7.2%
- inferred 1n: 20
- AB: 16,529,698
- AABB: 2,825,354

### T2 merge

- total pairs: 52,835,318
- error fraction: about 55.1%
- inferred 1n: 20
- AB: 20,312,886
- AABB: 1,702,456

### T3 merge

- total pairs: 56,207,120
- error fraction: about 58.1%
- inferred 1n: 22
- AB: 21,858,988

The much higher masked/error fraction in the cross-side merges is itself informative, but the called class proportions must be interpreted in light of cutoff bias.

## Pegasus operational lessons

- Give concurrent FastK jobs unique temporary directories.
- Directly prepend the known environment bin directory when concurrent micromamba invocations would contend for locks.
- FastK histogram files are binary; convert them with Histex before text processing.
- Do not infer Smudgeplot suitability from BAM mean depth alone.
- Automatic 1n coverage can choose a harmonic or zero; inspect the grid and the raw density before interpreting class labels.

## Plotting

The portable viridis renderer reads a .smu file directly and uses:

x = min(A,B) / (A+B)

y = A+B

weight = k-mer-pair count

This raw-density renderer is deliberately separate from Smudgeplot's class caller so that geometry can be inspected without silently inheriting a possibly incorrect ploidy/harmonic label.
