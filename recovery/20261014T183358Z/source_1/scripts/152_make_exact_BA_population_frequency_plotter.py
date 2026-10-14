#!/usr/bin/env python3

from pathlib import Path


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

SOURCE = (
    PROJECT
    / "scripts"
    / "151_plot_population_frequency_plus_LR_BWPK.R"
)

TARGET = (
    PROJECT
    / "scripts"
    / "152_plot_population_frequency_plus_LR_BWPK_exact_BA.R"
)


if not SOURCE.is_file():
    raise RuntimeError(
        f"Missing working source plotter: {SOURCE}"
    )

text = SOURCE.read_text()


# ============================================================
# 1. PUT THIS DERIVATIVE IN ITS OWN OUTPUT DIRECTORY
# ============================================================

old_output = '''    "population_frequency_plus_LR_BWPK"
)'''

new_output = '''    "population_frequency_plus_LR_BWPK_exact_BA"
)'''

if text.count(old_output) != 1:
    raise RuntimeError(
        "Could not uniquely identify the 151 output directory."
    )

text = text.replace(
    old_output,
    new_output,
    1
)


# ============================================================
# 2. LOAD THE EXACT, VALIDATED BA MANIFEST
# ============================================================

one_region_marker = '''# ============================================================
# ONE REGION
# ============================================================'''

if text.count(one_region_marker) != 1:
    raise RuntimeError(
        "Could not uniquely find the ONE REGION section."
    )

manifest_block = r'''
# ============================================================
# EXACT BA COMPARATIVE COHORTS
#
# These are the validated final comparative sample sets:
#
#   Band    163
#   Ray     368
#   Dennis  369
#
# Population identity is also taken from this manifest.
# Metadata therefore cannot introduce samples or split the
# exact BA populations into arbitrary additional rows.
# ============================================================

exact_ba_manifest_file <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.tsv"
)

if (!file.exists(exact_ba_manifest_file)) {
    stop(
        paste(
            "Missing validated exact BA manifest:",
            exact_ba_manifest_file
        )
    )
}

exact_ba_manifest <- read.delim(
    exact_ba_manifest_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_exact_columns <- c(
    "BA_pattern",
    "sample",
    "population",
    "taxon",
    "country",
    "data_type"
)

missing_exact_columns <- setdiff(
    required_exact_columns,
    names(exact_ba_manifest)
)

if (length(missing_exact_columns) > 0) {
    stop(
        paste(
            "Exact BA manifest is missing:",
            paste(
                missing_exact_columns,
                collapse = ", "
            )
        )
    )
}

for (
    column
    in required_exact_columns
) {
    exact_ba_manifest[[column]] <- trimws(
        as.character(
            exact_ba_manifest[[column]]
        )
    )
}

if (
    anyDuplicated(
        paste(
            exact_ba_manifest$BA_pattern,
            exact_ba_manifest$sample,
            sep = "||"
        )
    )
) {
    stop(
        paste(
            "Exact BA manifest contains duplicate",
            "pattern/sample combinations."
        )
    )
}

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
    any(
        exact_ba_observed_n !=
            exact_ba_expected_n
    )
) {
    stop(
        paste(
            "Exact BA cohort-size validation failed:",
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
    "\nVALIDATED EXACT BA COHORTS\n"
)

for (
    pattern_name
    in names(exact_ba_expected_n)
) {

    cat(
        sprintf(
            "  %-6s %d butterflies; %d populations\n",
            pattern_name,
            exact_ba_observed_n[
                pattern_name
            ],
            length(
                unique(
                    exact_ba_manifest$population[
                        exact_ba_manifest$BA_pattern ==
                            pattern_name
                    ]
                )
            )
        )
    )
}

cat("\n")
'''

text = text.replace(
    one_region_marker,
    manifest_block
    + "\n\n"
    + one_region_marker,
    1
)


# ============================================================
# 3. REPLACE THE OLD METADATA-DRIVEN COMPARATIVE POPULATION
#    ASSIGNMENT WITH EXACT BA-COHORT FILTERING
# ============================================================

start_marker = '''    comparative <- comparative[
        comparative$sample_id %in% meta$sample,
        ,
        drop = FALSE
    ]'''

end_marker = '''    comparative$GT_frequency <- genotype_frequency(
        comparative$display_GT
    )'''

start = text.find(
    start_marker
)

if start < 0:
    raise RuntimeError(
        "Could not find the start of the old comparative "
        "metadata-assignment block."
    )

end_start = text.find(
    end_marker,
    start
)

if end_start < 0:
    raise RuntimeError(
        "Could not find the end of the old comparative "
        "metadata-assignment block."
    )

end = (
    end_start
    + len(end_marker)
)


replacement = r'''    # --------------------------------------------------------
    # Exact BA cohort for this key region.
    #
    # Ray1 and Ray2 use the same Ray BA sample set.
    # Dennis1 and Dennis2 use the same Dennis BA sample set.
    # Band uses the final Analysis-142 Band BA sample set.
    # --------------------------------------------------------

    ba_pattern <- if (
        identical(
            panel_name,
            "Band"
        )
    ) {

        "Band"

    } else if (
        grepl(
            "^Ray",
            panel_name
        )
    ) {

        "Ray"

    } else if (
        grepl(
            "^Dennis",
            panel_name
        )
    ) {

        "Dennis"

    } else {

        stop(
            paste(
                "No exact BA cohort mapping for panel:",
                panel_name
            )
        )
    }


    cohort_manifest <- exact_ba_manifest[
        exact_ba_manifest$BA_pattern ==
            ba_pattern,
        ,
        drop = FALSE
    ]

    allowed_samples <- unique(
        cohort_manifest$sample
    )

    expected_sample_n <- as.integer(
        exact_ba_expected_n[
            ba_pattern
        ]
    )


    # --------------------------------------------------------
    # The established genotype-plot table must contain every
    # butterfly that actually entered the corresponding BA.
    #
    # Fail rather than silently calculating a frequency from a
    # partial BA cohort.
    # --------------------------------------------------------

    comparative_samples_in_genoplot <- unique(
        comparative$sample_id
    )

    missing_BA_samples <- setdiff(
        allowed_samples,
        comparative_samples_in_genoplot
    )

    if (length(missing_BA_samples) > 0) {

        stop(
            paste0(
                panel_name,
                ": ",
                length(missing_BA_samples),
                " exact ",
                ba_pattern,
                " BA samples are absent from the established ",
                "genotype-plot input. First missing IDs: ",
                paste(
                    head(
                        missing_BA_samples,
                        20
                    ),
                    collapse = ", "
                )
            )
        )
    }


    # --------------------------------------------------------
    # Now discard every comparative butterfly that was NOT in
    # the final BA cohort for this pattern.
    # --------------------------------------------------------

    comparative <- comparative[
        comparative$sample_id %in%
            allowed_samples,
        ,
        drop = FALSE
    ]

    recovered_sample_n <- length(
        unique(
            comparative$sample_id
        )
    )

    if (
        recovered_sample_n !=
            expected_sample_n
    ) {

        stop(
            paste0(
                panel_name,
                ": expected exactly ",
                expected_sample_n,
                " comparative ",
                ba_pattern,
                " BA butterflies after filtering, but recovered ",
                recovered_sample_n,
                "."
            )
        )
    }


    # --------------------------------------------------------
    # Assign population directly from the exact BA manifest.
    #
    # This is deliberate:
    #
    #   * no metadata-driven inclusion
    #   * no within-species subgroup splitting
    #   * no extra populations
    #
    # Population rows therefore correspond exactly to the
    # population list in exact_BA_population_summary.tsv.
    # --------------------------------------------------------

    cohort_index <- match(
        comparative$sample_id,
        cohort_manifest$sample
    )

    if (any(is.na(cohort_index))) {
        stop(
            paste(
                panel_name,
                ": failed to map retained comparative samples",
                "back to exact BA manifest."
            )
        )
    }

    comparative$population_key <- (
        cohort_manifest$population[
            cohort_index
        ]
    )

    comparative$population_label <- (
        cohort_manifest$population[
            cohort_index
        ]
    )

    comparative$species_group <- ""

    comparative$taxon <- (
        cohort_manifest$taxon[
            cohort_index
        ]
    )

    comparative$country <- (
        cohort_manifest$country[
            cohort_index
        ]
    )


    # --------------------------------------------------------
    # Hard population-set validation.
    # --------------------------------------------------------

    expected_populations <- sort(
        unique(
            cohort_manifest$population
        )
    )

    observed_populations <- sort(
        unique(
            comparative$population_label
        )
    )

    if (
        !identical(
            expected_populations,
            observed_populations
        )
    ) {

        stop(
            paste0(
                panel_name,
                ": population-set validation failed.\n",
                "Expected but absent: ",
                paste(
                    setdiff(
                        expected_populations,
                        observed_populations
                    ),
                    collapse = ", "
                ),
                "\nUnexpected: ",
                paste(
                    setdiff(
                        observed_populations,
                        expected_populations
                    ),
                    collapse = ", "
                )
            )
        )
    }


    comparative$GT_frequency <- genotype_frequency(
        comparative$display_GT
    )


    cat(
        sprintf(
            paste0(
                "%-8s | BA=%-6s | ",
                "%d comparative butterflies | ",
                "%d populations\n"
            ),
            panel_name,
            ba_pattern,
            recovered_sample_n,
            length(
                observed_populations
            )
        )
    )'''

text = (
    text[:start]
    + replacement
    + text[end:]
)


# ============================================================
# 4. UPDATE THE CAPTION SO THE PLOT STATES THE RULE
# ============================================================

old_caption_piece = '''                "Population rows show the mean pattern-associated allele dosage",
                "among callable individuals: A/A=0, P/A=0.5, P/P=1.",'''

new_caption_piece = '''                "Population rows contain only butterflies used in the corresponding final comparative BA cohort",
                "and show mean pattern-associated allele dosage among callable individuals: A/A=0, P/A=0.5, P/P=1.",'''

if text.count(old_caption_piece) != 1:
    raise RuntimeError(
        "Could not uniquely identify the existing population caption."
    )

text = text.replace(
    old_caption_piece,
    new_caption_piece,
    1
)


# ============================================================
# 5. VERSION LABEL
# ============================================================

text = text.replace(
    "ANALYSIS 151 COMPLETE",
    "ANALYSIS 152 EXACT-BA POPULATION FIGURE COMPLETE",
    1
)


TARGET.write_text(
    text
)

print(
    "Created:"
)

print(
    TARGET
)

print()

print(
    "Source left unchanged:"
)

print(
    SOURCE
)

print()

print(
    "Derived output directory:"
)

print(
    PROJECT
    / "analysis/figures"
    / "population_frequency_plus_LR_BWPK_exact_BA"
)
