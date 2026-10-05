# Figure 2 panel scripts

This directory contains one plotting script per manuscript Figure 2 panel. Each script writes a vector PDF by default to:

`analysis/figure2/panels/Fig2<LETTER>.pdf`

All scripts accept `--project-root`; otherwise they use `$GYNANDRO_ROOT` if set, then fall back to the Pegasus project path.

## Panel map

| panel | script | purpose | principal inputs |
|---|---|---|---|
| A | `Fig2A_side_specific_sequencing.py` | side-specific sequencing schematic | none |
| B | `Fig2B_chromosome_coverage.py` | mirrored chromosome relative depth | chromosome coverage table |
| C | `Fig2C_smudge_controls.py` | diploid and pseudo-tetraploid Smudgeplot controls | final L=6 `.smu` files |
| D | `Fig2D_smudge_gynandromorph.py` | Left, Right, T2-merge and T3-merge Smudgeplots | Left/Right L=12; T2/T3 L=6 |
| E | `Fig2E_AF_marginal.py` | 1-D AF distribution by original GT, focal versus diploid control | ERR2298208/TN histogram cache |
| F | `Fig2F_founder_scenarios.py` | founder/haplotype cartoons with 20% reciprocal mixing | none |
| G | `Fig2G_simulated_AF.py` | four simulated paired-AF founder configurations | empirical-depth simulation caches |
| H | `Fig2H_observed_AF.py` | observed paired AF with 20% model overlay and original-GT marginals | 31,126,758-site hex table + histogram cache |

`common.py` contains shared plotting conventions, the Luka-Bleute palette loader, Smudgeplot density utilities, and the Figure 2 predictive-ring drawing convention.

## Locked Figure 2 choices

- C/D use ERR1051012 L=6, pseudo-tetraploid L=6, TN Left L=12, TN Right L=12, T2 merge L=6, T3 merge L=6.
- G/H use the same Luka-Bleute density palette.
- G compares `AB/A`, `AA/BC`, `AB/AC`, and `AB/CD` at 20% reciprocal read mixing.
- G/H show magenta expected centroids and thin 50% predictive outlines. These are predictive-region annotations, not parameter-confidence intervals.
- H requires exactly 31,126,758 retained genome-wide sites.
- H uses original GT only for marginal colour strata; observed AF itself comes from AD.
- AF and Smudgeplot density artists are kept vector in the PDFs.

## Generate every panel

```bash
GYNANDRO_ROOT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling \
  bash scripts/figures/figure2/run_Figure2_panels.sh
```

The runner stops on the first missing input or plotting error.
