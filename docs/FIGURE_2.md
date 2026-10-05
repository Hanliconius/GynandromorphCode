# Manuscript Figure 2

Figure 2 is maintained as eight independent vector-PDF panels rather than one monolithic plotting script. Canonical plotting code is under `scripts/figures/figure2/` and is named directly by panel (`Fig2A_...` through `Fig2H_...`).

## Evidence layers

A introduces the four-library Left/Right design. B is the chromosome-scale coverage diagnostic. C and D are alignment-free k-mer-pair diagnostics using the final lower-cutoff calibration. E summarizes one-dimensional observed AF behaviour relative to a diploid control. F defines the founder/haplotype configurations. G shows their simulated paired-AF distributions under 20% reciprocal mixing. H shows the observed paired-AF distribution with the same palette and overlay convention.

The final composite should preserve the distinction between direct observations (B, C, D, E, H) and model illustrations/simulations (F, G and the overlay annotations on H).

## Observed AF scope

The genome-wide paired-AF analysis starts from 35,195,121 control-discovered biallelic SNP records. It retains 31,126,758 positions with usable AD >= 10 in both `TN_Left_merged` and `TN_Right_merged`. Original focal GT is not used to define observed AF; it is used only to colour/category the marginal distributions.

## Smudgeplot scope

Final Figure 2 cutoffs:

- ERR1051012 diploid control: L=6
- ERR1307111 + ERR1307112 pseudo-tetraploid benchmark: L=6
- TN Left: L=12
- TN Right: L=12
- TN T2 merge: L=6
- TN T3 merge: L=6

The pseudo-tetraploid is a merged two-diploid benchmark, not an organismal tetraploid.

## Predictive-ring convention

Panels G/H use thin magenta 50% predictive-region annotations at the expected mixed-AF centroids. The canonical Figure 2 helper draws these from the beta-binomial spread convention used in the figure (illustrative depth 40 reads/side, concentration 22). They are predictive regions, not confidence intervals on fitted parameters. The panels should be visually checked against the historical Pegasus analytical-overlay implementation before final submission.
