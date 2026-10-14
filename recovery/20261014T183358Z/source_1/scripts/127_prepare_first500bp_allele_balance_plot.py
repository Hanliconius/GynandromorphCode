#!/usr/bin/env python3

from pathlib import Path
import re


source = Path(
    "scripts/125_plot_focal_vs_high_depth_WGS_allele_balance.R"
)

target = Path(
    "scripts/127_plot_first500bp_focal_vs_WGS_allele_balance.R"
)

if not source.is_file():
    raise RuntimeError(f"Missing source plotter: {source}")

text = source.read_text()

# Specific filenames must be replaced before the shared stem.
replacements = [
    (
        "optix_selected_region_focal_and_WGS_controls.long.tsv",
        "optix_first500bp_focal_and_WGS_controls.long.tsv",
    ),
    (
        "optix_selected_mixed_genotype_region.tsv",
        "optix_first500bp_selection.tsv",
    ),
    (
        "optix_focal_vs_high_depth_WGS_allele_balance.pdf",
        "optix_first500bp_focal_vs_high_depth_WGS_allele_balance.pdf",
    ),
    (
        "optix_focal_vs_high_depth_WGS_allele_balance.png",
        "optix_first500bp_focal_vs_high_depth_WGS_allele_balance.png",
    ),
    (
        "optix_focal_vs_high_depth_WGS_allele_balance",
        "optix_first500bp_focal_vs_high_depth_WGS_allele_balance",
    ),
    (
        '" biallelic SNPs in the original regional VCF"',
        '" biallelic SNPs in the first 500-bp span of the broad VCF"',
    ),
]

for old, new in replacements:
    if old not in text:
        raise RuntimeError(
            f"Expected plotter text was not found: {old}"
        )

    text = text.replace(old, new, 1)

caption_pattern = re.compile(
    r'caption = paste\(\s*'
    r'"The three comparative controls are the highest-median-depth".*?'
    r'\)\s*'
    r'\)\s*\+\s*'
    r'theme_minimal',
    flags=re.DOTALL,
)

caption_replacement = '''caption = paste(
            "The three controls are the same highest-median-depth",
            "non-CS WGS samples used in the regional comparison.",
            "Every biallelic SNP in the first 500-bp span of the",
            "broad VCF is displayed; no BA or focal-genotype filter",
            "was applied.",
            "Colored points pass DP >= 2 and GQ >= 20;",
            "gray crosses retain allele-depth information but fail QC."
        )
    ) +

    theme_minimal'''

text, caption_replacements = caption_pattern.subn(
    caption_replacement,
    text,
    count=1,
)

if caption_replacements != 1:
    raise RuntimeError(
        "Could not replace the original figure caption."
    )

target.write_text(text)

print(f"Created: {target}")
