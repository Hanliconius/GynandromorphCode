#!/usr/bin/env python3

from pathlib import Path
import re
import textwrap


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

SOURCE = (
    PROJECT
    / "scripts"
    / "150_plot_key_region_population_context_genoplots.R"
)

TARGET = (
    PROJECT
    / "scripts"
    / "151d_plot_key_region_population_context_genoplots_exact_BA.R"
)


if not SOURCE.is_file():
    raise RuntimeError(
        f"Missing source plotter: {SOURCE}"
    )

text = SOURCE.read_text()


# ============================================================
# 1. GIVE THE DERIVED OUTPUTS THEIR OWN NAMES
# ============================================================
#
# Do not overwrite the earlier population-context plots.
# This replaces only the output-name stem.
# ============================================================

old_stem = (
    "Optix_key_regions_population_context_genoplots"
)

new_stem = (
    "Optix_key_regions_population_context_genoplots."
    "exact_BA_cohorts"
)

stem_count = text.count(old_stem)

if stem_count < 1:
    raise RuntimeError(
        "Could not find the expected output-name stem in script 150."
    )

text = text.replace(
    old_stem,
    new_stem,
)


# ============================================================
# 2. LOAD THE EXACT BA SAMPLE MANIFEST
# ============================================================
#
# Insert this before the existing MAIN LOOP.
#
# This manifest was independently built from:
#
#   Band:
#     final Analysis-142 WGS+CS retained cohort
#
#   Ray:
#     final Ray all-comparative BA cohort
#
#   Dennis:
#     final Dennis all-comparative BA cohort
#
# Hard expected counts:
#
#   Band    163
#   Ray     368
#   Dennis  369
# ============================================================

main_loop_marker = "# MAIN LOOP"

main_loop_pos = text.find(
    main_loop_marker
)

if main_loop_pos < 0:
    raise RuntimeError(
        "Could not find '# MAIN LOOP' in script 150."
    )

manifest_block = r'''
# ============================================================
# EXACT COMPARATIVE BA COHORTS FOR POPULATION CONTEXT
# ============================================================
#
# IMPORTANT:
#
# Population frequencies below must be calculated from exactly
# the butterflies used in the corresponding comparative BA
# analysis.
#
# The focal samples and BWPK individuals are NOT members of
# these comparative population cohorts.
# ============================================================

exact_ba_manifest_file <- file.path(
    output_dir,
    "exact_BA_comparative_sample_manifest.tsv"
)

if (!file.exists(exact_ba_manifest_file)) {
    stop(
        paste0(
            "Missing exact BA comparative sample manifest: ",
            exact_ba_manifest_file
        )
    )
}

exact_ba_manifest <- read.delim(
    exact_ba_manifest_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
)

required_exact_ba_columns <- c(
    "BA_pattern",
    "sample"
)

missing_exact_ba_columns <- setdiff(
    required_exact_ba_columns,
    names(exact_ba_manifest)
)

if (length(missing_exact_ba_columns) > 0) {
    stop(
        paste(
            "Exact BA manifest is missing columns:",
            paste(
                missing_exact_ba_columns,
                collapse = ", "
            )
        )
    )
}

exact_ba_manifest$BA_pattern <- trimws(
    as.character(
        exact_ba_manifest$BA_pattern
    )
)

exact_ba_manifest$sample <- trimws(
    as.character(
        exact_ba_manifest$sample
    )
)

exact_ba_expected_n <- c(
    Band = 163L,
    Ray = 368L,
    Dennis = 369L
)

exact_ba_observed_n <- vapply(
    names(exact_ba_expected_n),
    function(pattern_name) {
        length(
            unique(
                exact_ba_manifest$sample[
                    exact_ba_manifest$BA_pattern ==
                        pattern_name
                ]
            )
        )
    },
    integer(1)
)

if (
    !identical(
        unname(exact_ba_observed_n),
        unname(exact_ba_expected_n)
    )
) {
    stop(
        paste(
            "Exact BA manifest count validation failed.",
            paste(
                names(exact_ba_expected_n),
                exact_ba_observed_n,
                sep = "=",
                collapse = ", "
            )
        )
    )
}

cat(
    "\nExact BA population-context cohorts loaded:\n"
)

for (
    pattern_name
    in names(exact_ba_expected_n)
) {
    cat(
        sprintf(
            "  %-6s %d butterflies\n",
            pattern_name,
            exact_ba_observed_n[
                pattern_name
            ]
        )
    )
}

'''

text = (
    text[:main_loop_pos]
    + manifest_block
    + "\n"
    + text[main_loop_pos:]
)


# ============================================================
# 3. FILTER COMPARATIVE_QC IMMEDIATELY BEFORE THE EXISTING
#    POPULATION-FREQUENCY AGGREGATION
# ============================================================
#
# The original script contains the calculation:
#
# aggregate(
#     plot_value ~ population_use + site_id + pos,
#     data = comparative_qc,
#     FUN = mean
# )
#
# We leave that calculation unchanged.
#
# We only restrict comparative_qc to the sample IDs belonging
# to the BA cohort for the current pattern before it happens.
# ============================================================

aggregate_pattern = re.compile(
    r'''
    ^
    (?P<indent>[ \t]*)
    (?P<lhs>
        [A-Za-z][A-Za-z0-9_.]*
        [ \t]*<-[ \t]*
    )
    aggregate
    \(
        \s*
        plot_value
        \s*~\s*
        population_use
        \s*\+\s*
        site_id
        \s*\+\s*
        pos
        \s*,\s*
        data
        \s*=\s*
        comparative_qc
        \s*,\s*
        FUN
        \s*=\s*
        mean
        \s*
    \)
    ''',
    flags=(
        re.MULTILINE
        | re.DOTALL
        | re.VERBOSE
    ),
)

matches = list(
    aggregate_pattern.finditer(
        text
    )
)

if len(matches) != 1:
    raise RuntimeError(
        "Expected exactly one population-frequency aggregate "
        f"statement, found {len(matches)}."
    )

match = matches[0]

indent = match.group(
    "indent"
)

filter_block_raw = r'''
# ------------------------------------------------------------
# Restrict population context to the exact comparative sample
# cohort used to calculate BA for this pattern.
#
# Band        -> Band BA cohort
# Ray1/Ray2   -> Ray BA cohort
# Dennis1/2   -> Dennis BA cohort
#
# Filtering occurs AFTER genotype QC but BEFORE population
# aggregation, so unavailable genotypes remain unavailable;
# they are not converted into reference calls.
# ------------------------------------------------------------

ba_pattern_for_region <- if (
    identical(
        region_name,
        "Band"
    )
) {
    "Band"
} else if (
    grepl(
        "^Ray",
        region_name
    )
) {
    "Ray"
} else if (
    grepl(
        "^Dennis",
        region_name
    )
) {
    "Dennis"
} else {
    NA_character_
}

if (is.na(ba_pattern_for_region)) {
    stop(
        paste0(
            "No BA-cohort mapping defined for region: ",
            region_name
        )
    )
}

allowed_ba_samples <- unique(
    exact_ba_manifest$sample[
        exact_ba_manifest$BA_pattern ==
            ba_pattern_for_region
    ]
)

expected_ba_sample_n <- as.integer(
    exact_ba_expected_n[
        ba_pattern_for_region
    ]
)

if (
    length(
        allowed_ba_samples
    ) !=
        expected_ba_sample_n
) {
    stop(
        paste0(
            "Wrong BA whitelist size for ",
            region_name,
            " / ",
            ba_pattern_for_region,
            ": expected ",
            expected_ba_sample_n,
            ", observed ",
            length(
                allowed_ba_samples
            )
        )
    )
}

comparative_samples_before_BA_filter <- unique(
    as.character(
        comparative_qc$sample_id
    )
)

comparative_qc <- comparative_qc[
    as.character(
        comparative_qc$sample_id
    ) %in%
        allowed_ba_samples,
    ,
    drop = FALSE
]

comparative_samples_after_BA_filter <- unique(
    as.character(
        comparative_qc$sample_id
    )
)

cat(
    sprintf(
        paste0(
            "\n%-8s -> %-6s BA cohort: ",
            "%d allowed butterflies; ",
            "%d with >=1 QC-passing genotype ",
            "represented in comparative_qc\n"
        ),
        region_name,
        ba_pattern_for_region,
        length(
            allowed_ba_samples
        ),
        length(
            comparative_samples_after_BA_filter
        )
    )
)

'''

filter_block = textwrap.indent(
    textwrap.dedent(
        filter_block_raw
    ).strip("\n"),
    indent,
)

replacement = (
    filter_block
    + "\n\n"
    + match.group(0)
)

text = (
    text[:match.start()]
    + replacement
    + text[match.end():]
)


# ============================================================
# 4. FINAL SAFETY CHECKS
# ============================================================

if (
    "exact_BA_comparative_sample_manifest.tsv"
    not in text
):
    raise RuntimeError(
        "Manifest-loading block was not inserted."
    )

if text.count(
    "ba_pattern_for_region"
) < 2:
    raise RuntimeError(
        "BA cohort filter block was not inserted correctly."
    )

if old_stem in text:
    raise RuntimeError(
        "Old output stem remains in derived script."
    )

TARGET.write_text(
    text
)

print(
    "Created derived plotter:"
)

print(
    TARGET
)

print()

print(
    "Source script left unchanged:"
)

print(
    SOURCE
)

print()

print(
    "Output stem:"
)

print(
    new_stem
)

print()

print(
    "Inserted exact BA cohort mapping:"
)

print(
    "  Band          -> Band   (163)"
)

print(
    "  Ray1 / Ray2   -> Ray    (368)"
)

print(
    "  Dennis1 / 2   -> Dennis (369)"
)
