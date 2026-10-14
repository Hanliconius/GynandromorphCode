#!/usr/bin/env python3

from pathlib import Path
import re


source_python = Path(
    "scripts/120a_extract_all_comparative_R2_site_genotypes.py"
)

target_python = Path(
    "scripts/121a_extract_all_R_BA90_genotypes.py"
)

source_r = Path(
    "scripts/120b_plot_all_comparative_R2_site_genotypes.R"
)

target_r = Path(
    "scripts/121b_plot_all_R_BA90_genotypes.R"
)

for required in [source_python, source_r]:
    if not required.is_file():
        raise RuntimeError(f"Missing source script: {required}")


# ------------------------------------------------------------
# Python extraction script
# ------------------------------------------------------------

text = source_python.read_text()

text = text.replace(
    '"ray2-1174547-1187833.vcf"',
    '"Optix_Hmr18-938195-1239421_biallelic.vcf"',
)

text = text.replace(
    "optix_R2_all_comparative_selected",
    "optix_R_all_BA90_comparative_selected",
)

selection_pattern = re.compile(
    r"sites = strict\.loc\[\s*"
    r".*?"
    r"\]\.copy\(\)\s*"
    r"sites = sites\.drop_duplicates",
    flags=re.DOTALL,
)

selection_replacement = '''
if "minimum_BA_all_analyses" not in strict.columns:
    raise RuntimeError(
        "Strict-site table lacks minimum_BA_all_analyses."
    )

strict["minimum_BA_all_analyses"] = pd.to_numeric(
    strict["minimum_BA_all_analyses"],
    errors="coerce",
)

sites = strict.loc[
    (strict["module"] == "Ray")
    & (strict["chrom"] == "Hmr18")
    & (strict["minimum_BA_all_analyses"] >= 0.90)
].copy()

sites = sites.drop_duplicates'''

text, replacements = selection_pattern.subn(
    selection_replacement,
    text,
    count=1,
)

if replacements != 1:
    raise RuntimeError(
        "Could not replace the R2 site-selection block."
    )

text = text.replace(
    '"region": "R2"',
    '"region": "broad_optix_BA90"',
)

text = text.replace(
    '"R2",\n            "comparative_selected_SNPs"',
    '"broad_optix_BA90",\n'
    '            "comparative_selected_SNPs"',
)

text = text.replace(
    'print("ALL COMPARATIVE-SELECTED R2 GENOTYPES")',
    'print("ALL COMPARATIVE-SELECTED RAY BA90 GENOTYPES")',
)

text = text.replace(
    '"No comparative-selected Ray sites were found in R2."',
    '"No comparative-selected Ray sites with BA >= 0.90 were found."',
)

compile(text, str(target_python), "exec")
target_python.write_text(text)


# ------------------------------------------------------------
# R plotting script
# ------------------------------------------------------------

r_text = source_r.read_text()

r_text = r_text.replace(
    "optix_R2_all_comparative_selected",
    "optix_R_all_BA90_comparative_selected",
)

r_text = r_text.replace(
    "optix_R2_all_comparative_selected_genoplot",
    "optix_R_all_BA90_comparative_selected_genoplot",
)

r_text = r_text.replace(
    "Ray2 comparative-associated genotype plot",
    "Broad optix Ray-associated BA90 genotype plot",
)

subtitle_pattern = re.compile(
    r"subtitle = paste0\(\s*"
    r".*?"
    r"\),\s*"
    r'x = "Hmr18 position"',
    flags=re.DOTALL,
)

subtitle_replacement = '''subtitle = paste0(
            site_count,
            " comparative-selected Ray SNPs with minimum BA >= 0.90; ",
            sample_count,
            " comparative individuals. ",
            "Focal L/R agreement was not required."
        ),
        x = "Hmr18 position"'''

r_text, replacements = subtitle_pattern.subn(
    subtitle_replacement,
    r_text,
    count=1,
)

if replacements != 1:
    raise RuntimeError(
        "Could not replace the R plot subtitle."
    )

r_text = r_text.replace(
    "optix_R2_all_comparative_selected_genotypes.pdf",
    "optix_R_all_BA90_comparative_selected_genotypes.pdf",
)

r_text = r_text.replace(
    "R2 ALL COMPARATIVE-SELECTED GENOTYPE PLOT",
    "BROAD RAY BA90 COMPARATIVE-SELECTED GENOTYPE PLOT",
)

target_r.write_text(r_text)

print(f"Created and syntax-checked: {target_python}")
print(f"Created: {target_r}")
