#!/usr/bin/env python3

from pathlib import Path


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

SOURCE = (
    PROJECT
    / "scripts"
    / "152_plot_population_frequency_plus_LR_BWPK_exact_BA.R"
)

TARGET = (
    PROJECT
    / "scripts"
    / "153_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref.R"
)


if not SOURCE.is_file():
    raise RuntimeError(
        f"Missing working Analysis-152 plotter: {SOURCE}"
    )

text = SOURCE.read_text()


# ============================================================
# 1. NEW OUTPUT DIRECTORY
# ============================================================

old_output = '''    "population_frequency_plus_LR_BWPK_exact_BA"
)'''

new_output = '''    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref"
)'''

if text.count(old_output) != 1:
    raise RuntimeError(
        "Could not uniquely identify the Analysis-152 "
        "output directory."
    )

text = text.replace(
    old_output,
    new_output,
    1
)


# ============================================================
# 2. INSERT SITE FILTER IMMEDIATELY BEFORE THE EXISTING
#    POPULATION + INDIVIDUAL COMBINATION.
#
# By this point Analysis 152 has already:
#
#   * selected the exact BA cohort
#   * calculated population frequencies
#   * recovered TN_Left / TN_Right reciprocal rows
#   * recovered BWPK individual rows
#
# We therefore use raw_GT from exactly those underlying
# displayed butterflies to decide whether a physical SNP is
# monomorphic reference.
# ============================================================

marker = '''    # --------------------------------------------------------
    # Combine population-frequency and individual rows.
    # --------------------------------------------------------'''

if text.count(marker) != 1:
    raise RuntimeError(
        "Could not uniquely find the population/individual "
        "combination marker in Analysis 152."
    )


filter_block = r'''
    # ========================================================
    # REMOVE SITES MONOMORPHIC FOR THE RAW REFERENCE ALLELE
    #
    # Display scope:
    #
    #   1. exact comparative BA cohort for this pattern
    #   2. TN_Left reciprocal-mixture row
    #   3. TN_Right reciprocal-mixture row
    #   4. all Butterfly World piano-key individuals
    #
    # A site is removed when:
    #
    #   * at least one displayed sample has a callable raw GT
    #   * every callable raw GT is 0/0
    #
    # Missing/uncertain calls neither create nor erase evidence
    # of a non-reference allele.
    #
    # All-missing sites are NOT classified as fixed reference.
    # ========================================================

    if (!("raw_GT" %in% names(d))) {
        stop(
            paste0(
                panel_name,
                ": established genotype-plot input does not ",
                "contain raw_GT; refusing to infer reference ",
                "state from biological display genotype."
            )
        )
    }


    # --------------------------------------------------------
    # Exactly the rows contributing to the displayed figure.
    #
    # For comparative butterflies, allowed_samples has already
    # been frozen to the exact Band/Ray/Dennis BA cohort by
    # Analysis 152.
    #
    # For focal samples we deliberately use ONLY the reciprocal
    # rows, matching the figure itself.
    # --------------------------------------------------------

    in_exact_comparative_cohort <- (
        !is_bwpk &
        !is_focal_group &
        !is_left_recip &
        !is_right_recip &
        d$sample_id %in% allowed_samples
    )

    in_display_scope <- (
        in_exact_comparative_cohort |
        is_bwpk |
        is_left_recip |
        is_right_recip
    )

    raw_scope <- d[
        in_display_scope,
        c(
            "site_key",
            "raw_GT"
        ),
        drop = FALSE
    ]


    # --------------------------------------------------------
    # Canonicalize raw GT.
    # --------------------------------------------------------

    raw_scope$raw_GT_canonical <- gsub(
        "\\|",
        "/",
        trimws(
            as.character(
                raw_scope$raw_GT
            )
        )
    )

    raw_scope$raw_GT_canonical[
        raw_scope$raw_GT_canonical == "1/0"
    ] <- "0/1"

    raw_scope$raw_GT_canonical[
        !raw_scope$raw_GT_canonical %in% c(
            "0/0",
            "0/1",
            "1/1"
        )
    ] <- NA_character_


    # --------------------------------------------------------
    # One audit row per physical SNP.
    # --------------------------------------------------------

    raw_split <- split(
        raw_scope$raw_GT_canonical,
        raw_scope$site_key
    )

    fixed_ref_audit <- do.call(
        rbind,
        lapply(
            names(raw_split),
            function(site_key_value) {

                gt <- raw_split[[
                    site_key_value
                ]]

                callable <- !is.na(gt)

                callable_n <- sum(
                    callable
                )

                ref_n <- sum(
                    gt == "0/0",
                    na.rm = TRUE
                )

                het_n <- sum(
                    gt == "0/1",
                    na.rm = TRUE
                )

                alt_hom_n <- sum(
                    gt == "1/1",
                    na.rm = TRUE
                )

                nonref_n <- (
                    het_n +
                    alt_hom_n
                )

                data.frame(
                    site_key = site_key_value,
                    callable_n = callable_n,
                    raw_ref_hom_n = ref_n,
                    raw_het_n = het_n,
                    raw_alt_hom_n = alt_hom_n,
                    raw_nonref_n = nonref_n,

                    fixed_reference = (
                        callable_n > 0 &
                        nonref_n == 0
                    ),

                    stringsAsFactors = FALSE
                )
            }
        )
    )


    # Include every original site in the audit, even if somehow
    # no displayed row existed for it.
    site_audit <- merge(
        sites[
            ,
            c(
                "site_key",
                "chrom",
                "pos",
                "ref",
                "alt"
            ),
            drop = FALSE
        ],
        fixed_ref_audit,
        by = "site_key",
        all.x = TRUE,
        sort = FALSE
    )

    site_audit$fixed_reference[
        is.na(
            site_audit$fixed_reference
        )
    ] <- FALSE

    site_audit$callable_n[
        is.na(
            site_audit$callable_n
        )
    ] <- 0L

    for (
        count_column
        in c(
            "raw_ref_hom_n",
            "raw_het_n",
            "raw_alt_hom_n",
            "raw_nonref_n"
        )
    ) {

        site_audit[[count_column]][
            is.na(
                site_audit[[count_column]]
            )
        ] <- 0L
    }


    fixed_ref_site_keys <- site_audit$site_key[
        site_audit$fixed_reference
    ]

    keep_site_keys <- site_audit$site_key[
        !site_audit$fixed_reference
    ]

    original_site_n <- nrow(
        sites
    )

    fixed_ref_removed_n <- length(
        fixed_ref_site_keys
    )

    retained_site_n <- length(
        keep_site_keys
    )


    if (retained_site_n == 0) {
        stop(
            paste0(
                panel_name,
                ": filtering fixed-reference positions removed ",
                "every SNP."
            )
        )
    }


    # --------------------------------------------------------
    # Save explicit audit BEFORE changing the objects used for
    # plotting.
    # --------------------------------------------------------

    fixed_ref_audit_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            panel_name,
            "_fixed_reference_site_filter.tsv"
        )
    )

    site_audit <- site_audit[
        order(
            site_audit$pos,
            site_audit$ref,
            site_audit$alt
        ),
        ,
        drop = FALSE
    ]

    write.table(
        site_audit,
        file = fixed_ref_audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    # --------------------------------------------------------
    # Physically remove those columns from EVERY plot object.
    # --------------------------------------------------------

    sites <- sites[
        sites$site_key %in%
            keep_site_keys,
        ,
        drop = FALSE
    ]

    sites <- sites[
        order(
            sites$pos,
            sites$ref,
            sites$alt
        ),
        ,
        drop = FALSE
    ]

    # Re-number remaining SNPs consecutively so removed positions
    # do not leave visual gaps.
    sites$site_index <- seq_len(
        nrow(sites)
    )


    pop_grid <- pop_grid[
        pop_grid$site_key %in%
            keep_site_keys,
        ,
        drop = FALSE
    ]

    pop_grid$site_index <- sites$site_index[
        match(
            pop_grid$site_key,
            sites$site_key
        )
    ]


    focal <- focal[
        focal$site_key %in%
            keep_site_keys,
        ,
        drop = FALSE
    ]

    focal$site_index <- sites$site_index[
        match(
            focal$site_key,
            sites$site_key
        )
    ]


    bwpk <- bwpk[
        bwpk$site_key %in%
            keep_site_keys,
        ,
        drop = FALSE
    ]

    bwpk$site_index <- sites$site_index[
        match(
            bwpk$site_key,
            sites$site_key
        )
    ]


    cat(
        sprintf(
            paste0(
                "%-8s | SNP columns: ",
                "%d original -> %d retained; ",
                "%d fixed-reference positions removed\n"
            ),
            panel_name,
            original_site_n,
            retained_site_n,
            fixed_ref_removed_n
        )
    )


    # --------------------------------------------------------
    # Combine population-frequency and individual rows.
    # --------------------------------------------------------'''

text = text.replace(
    marker,
    filter_block,
    1
)


# ============================================================
# 3. UPDATE CAPTION
# ============================================================

old_caption = '''                "SNP columns are equal width and ordered by genomic position."
            )'''

new_caption = '''                "Positions at which every callable displayed raw genotype is 0/0 are omitted.",
                "SNP columns are equal width and ordered by genomic position."
            )'''

if text.count(old_caption) != 1:
    raise RuntimeError(
        "Could not uniquely identify the existing caption ending."
    )

text = text.replace(
    old_caption,
    new_caption,
    1
)


# ============================================================
# 4. VERSION LABEL
# ============================================================

old_version = (
    "ANALYSIS 152 EXACT-BA POPULATION FIGURE COMPLETE"
)

new_version = (
    "ANALYSIS 153 EXACT-BA NON-FIXED-REFERENCE "
    "POPULATION FIGURE COMPLETE"
)

if text.count(old_version) != 1:
    raise RuntimeError(
        "Could not uniquely identify Analysis-152 version label."
    )

text = text.replace(
    old_version,
    new_version,
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
    "New output directory:"
)

print(
    PROJECT
    / "analysis/figures"
    / "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref"
)

print()

print(
    "Filter rule:"
)

print(
    "  remove SNP if >=1 displayed raw GT is callable "
    "and all callable displayed raw GTs are 0/0"
)
