#!/usr/bin/env python3

from pathlib import Path

project = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

source = (
    project
    / "scripts"
    / "152_plot_population_frequency_plus_LR_BWPK_exact_BA.R"
)

target = (
    project
    / "scripts"
    / "153_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref.R"
)

if not source.is_file():
    raise RuntimeError(
        f"Missing working Analysis-152 source: {source}"
    )

text = source.read_text()


# ============================================================
# NEW OUTPUT DIRECTORY
# ============================================================

old = '"population_frequency_plus_LR_BWPK_exact_BA"'
new = '"population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref"'

if text.count(old) != 1:
    raise RuntimeError(
        f"Expected one Analysis-152 output-directory string; "
        f"found {text.count(old)}."
    )

text = text.replace(
    old,
    new,
    1
)


# ============================================================
# INSERT FILTER DIRECTLY BEFORE pop_grid/focal/BWPK ARE
# COMBINED FOR PLOTTING
# ============================================================

marker = """    # --------------------------------------------------------
    # Combine population-frequency and individual rows.
    # --------------------------------------------------------
"""

if text.count(marker) != 1:
    raise RuntimeError(
        "Could not uniquely find the existing "
        "'Combine population-frequency and individual rows' marker."
    )


filter_code = r'''    # ========================================================
    # DROP POSITIONS FIXED FOR RAW REFERENCE
    #
    # The scope is exactly what this figure displays:
    #
    #   * exact BA comparative cohort for this pattern
    #   * TN_Left reciprocal-mixture row
    #   * TN_Right reciprocal-mixture row
    #   * all BWPK individuals
    #
    # A position is removed only when:
    #
    #   * >=1 displayed genotype is callable
    #   * every callable displayed raw genotype is 0/0
    #
    # Missing / N / Uncertain calls are ignored.
    # ========================================================

    if (!("raw_GT" %in% names(d))) {
        stop(
            paste0(
                panel_name,
                ": raw_GT is absent from established ",
                "genotype-plot input."
            )
        )
    }


    # --------------------------------------------------------
    # Which rows belong to samples actually represented in the
    # plot?
    # --------------------------------------------------------

    in_display_scope <- (
        d$sample_id %in% allowed_samples |
        is_bwpk |
        is_left_recip |
        is_right_recip
    )


    site_calls <- d[
        in_display_scope,
        c(
            "site_key",
            "display_GT",
            "raw_GT"
        ),
        drop = FALSE
    ]


    # --------------------------------------------------------
    # Only genotypes that the established genotype plot regards
    # as callable may contribute to the fixed-reference test.
    #
    # This avoids allowing a low-quality raw 0/1 call hidden as
    # N/Uncertain to keep an otherwise fixed-reference SNP.
    # --------------------------------------------------------

    site_calls$callable <- (
        site_calls$display_GT %in%
            c(
                "A/A",
                "P/A",
                "P/P"
            )
    )


    site_calls$raw_GT_canonical <- gsub(
        "\\|",
        "/",
        trimws(
            as.character(
                site_calls$raw_GT
            )
        )
    )

    site_calls$raw_GT_canonical[
        site_calls$raw_GT_canonical == "1/0"
    ] <- "0/1"

    site_calls$raw_GT_canonical[
        !site_calls$raw_GT_canonical %in%
            c(
                "0/0",
                "0/1",
                "1/1"
            )
    ] <- NA_character_


    # A raw genotype is usable only if BOTH:
    #
    #   1. the established plot considers it callable
    #   2. raw_GT resolves to 0/0, 0/1 or 1/1
    #
    site_calls$usable <- (
        site_calls$callable &
        !is.na(
            site_calls$raw_GT_canonical
        )
    )


    calls_by_site <- split(
        site_calls,
        site_calls$site_key
    )


    fixed_ref_rows <- lapply(
        names(calls_by_site),
        function(site_key_value) {

            z <- calls_by_site[
                [site_key_value]
            ]

            usable_gt <- z$raw_GT_canonical[
                z$usable
            ]

            callable_n <- length(
                usable_gt
            )

            ref_n <- sum(
                usable_gt == "0/0"
            )

            het_n <- sum(
                usable_gt == "0/1"
            )

            alt_hom_n <- sum(
                usable_gt == "1/1"
            )

            data.frame(
                site_key = site_key_value,

                callable_n = callable_n,

                ref_hom_n = ref_n,

                het_n = het_n,

                alt_hom_n = alt_hom_n,

                nonref_n = (
                    het_n +
                    alt_hom_n
                ),

                fixed_reference = (
                    callable_n > 0 &&
                    het_n == 0 &&
                    alt_hom_n == 0
                ),

                stringsAsFactors = FALSE
            )
        }
    )


    fixed_ref_audit <- do.call(
        rbind,
        fixed_ref_rows
    )


    # --------------------------------------------------------
    # Join to the complete original site list so the audit also
    # records any all-missing position.
    # --------------------------------------------------------

    fixed_ref_audit <- merge(
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


    fixed_ref_audit$fixed_reference[
        is.na(
            fixed_ref_audit$fixed_reference
        )
    ] <- FALSE


    count_columns <- c(
        "callable_n",
        "ref_hom_n",
        "het_n",
        "alt_hom_n",
        "nonref_n"
    )

    for (
        count_column
        in count_columns
    ) {

        fixed_ref_audit[[count_column]][
            is.na(
                fixed_ref_audit[[count_column]]
            )
        ] <- 0L
    }


    fixed_ref_audit <- fixed_ref_audit[
        order(
            fixed_ref_audit$pos,
            fixed_ref_audit$ref,
            fixed_ref_audit$alt
        ),
        ,
        drop = FALSE
    ]


    remove_site_keys <- fixed_ref_audit$site_key[
        fixed_ref_audit$fixed_reference
    ]

    keep_site_keys <- fixed_ref_audit$site_key[
        !fixed_ref_audit$fixed_reference
    ]


    original_site_n <- nrow(
        sites
    )

    removed_site_n <- length(
        remove_site_keys
    )

    retained_site_n <- length(
        keep_site_keys
    )


    if (retained_site_n < 1) {
        stop(
            paste0(
                panel_name,
                ": fixed-reference filtering removed every site."
            )
        )
    }


    # --------------------------------------------------------
    # Audit exactly what was removed.
    # --------------------------------------------------------

    audit_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            panel_name,
            "_fixed_reference_site_filter.tsv"
        )
    )

    write.table(
        fixed_ref_audit,
        file = audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    # --------------------------------------------------------
    # Filter and then REINDEX retained SNP columns.
    #
    # Thus fixed-reference positions vanish rather than leaving
    # blank gaps in the genotype plot.
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
                "%-8s | SNPs: %d original -> %d retained; ",
                "%d fixed-reference positions removed\n"
            ),
            panel_name,
            original_site_n,
            retained_site_n,
            removed_site_n
        )
    )


    # --------------------------------------------------------
    # Combine population-frequency and individual rows.
    # --------------------------------------------------------
'''

text = text.replace(
    marker,
    filter_code,
    1
)


# ============================================================
# UPDATE FIGURE CAPTION
# ============================================================

caption_old = '''                "SNP columns are equal width and ordered by genomic position."
            )'''

caption_new = '''                "Sites fixed for raw REF (0/0) across all callable displayed samples are omitted.",
                "SNP columns are equal width and ordered by genomic position."
            )'''

if text.count(caption_old) != 1:
    raise RuntimeError(
        "Could not uniquely locate Analysis-152 caption ending."
    )

text = text.replace(
    caption_old,
    caption_new,
    1
)


# ============================================================
# UPDATE VERSION STRING
# ============================================================

text = text.replace(
    "ANALYSIS 152 EXACT-BA POPULATION FIGURE COMPLETE",
    "ANALYSIS 153 NON-FIXED-REF EXACT-BA POPULATION FIGURE COMPLETE",
    1
)


# ============================================================
# WRITE ONLY AFTER EVERY PATCH HAS SUCCEEDED
# ============================================================

target.write_text(
    text
)

print(
    f"Created: {target}"
)

print(
    "Filter: remove sites with >=1 callable displayed genotype "
    "and all callable raw genotypes = 0/0."
)
