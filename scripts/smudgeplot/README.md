# Smudgeplot scripts

## plot_smudge_smu_viridis.py

Portable renderer for Smudgeplot .smu files.

The script reads three columns:

- A k-mer coverage
- B k-mer coverage
- number of k-mer pairs

and converts them to:

x = min(A,B) / (A+B)

y = A+B

weight = number of k-mer pairs

It then draws a weighted viridis density map with logarithmic density normalization.

Example:

python scripts/smudgeplot/plot_smudge_smu_viridis.py \
    --input path/to/kmerpairs.smu \
    --output-prefix results/sample_smudge \
    --title "Sample"

Optional arguments:

--ymax 120

--no-labels

Outputs are PDF, SVG, and PNG. PDF/SVG are preferred for archival figures.

## Interpretation warning

The renderer can place theoretical AB/AAB/AAAB/AABB reference positions from an estimated balanced-cloud coverage. Those labels are visual guides, not an independent ploidy fit.

The project demonstrated two failure modes that must be checked before biological interpretation:

1. automatic 1n coverage can choose a harmonic;
2. a high lower-count cutoff can remove the minor arm of imbalanced k-mer pairs.

See ../../docs/SMUDGEPLOT.md.
