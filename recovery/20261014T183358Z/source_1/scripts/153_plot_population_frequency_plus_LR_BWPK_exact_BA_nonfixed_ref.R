#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

dr_input_dir <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "krista_DR_regions_reciprocal_genoplots"
)

band_input <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_band_original_region_reciprocal_genoplot.long.tsv.gz"
)

metadata_file <- file.path(
    project,
    "data/derived",
    "region_panel_metadata.tsv"
)

output_dir <- file.path(
    project,
    "analysis/figures",
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PANELS
# ============================================================

panels <- data.frame(
    panel = c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    ),
    input = c(
        band_input,
        file.path(
            dr_input_dir,
            "optix_Ray1_original_region_reciprocal_genoplot.long.tsv.gz"
        ),
        file.path(
            dr_input_dir,
            "optix_Ray2_original_region_reciprocal_genoplot.long.tsv.gz"
        ),
        file.path(
            dr_input_dir,
            "optix_Dennis1_original_region_reciprocal_genoplot.long.tsv.gz"
        ),
        file.path(
            dr_input_dir,
            "optix_Dennis2_original_region_reciprocal_genoplot.long.tsv.gz"
        )
    ),
    stringsAsFactors = FALSE
)


# ============================================================
# PREVIOUS GENOTYPE-PLOT AESTHETIC
#
# Same palette family used in the true-vector Krista plots.
# Individual:
#   A/A = 0
#   P/A = 0.5
#   P/P = 1
#
# Population rows use the mean of that same dosage scale.
# ============================================================

AA_colour <- "#2C9BCB"
PA_colour <- "#F28E2B"
PP_colour <- "#B61F72"
missing_colour <- "#F2F2F2"
uncertain_colour <- "#777777"


# ============================================================
# HELPERS
# ============================================================

normalise_name <- function(
    x
) {
    tolower(
        gsub(
            "[^a-z0-9]+",
            "",
            as.character(x)
        )
    )
}


find_col <- function(
    d,
    candidates,
    required = TRUE
) {

    normal_names <- normalise_name(
        names(d)
    )

    for (candidate in candidates) {

        hit <- which(
            normal_names ==
                normalise_name(candidate)
        )

        if (length(hit) == 1) {
            return(
                names(d)[hit]
            )
        }
    }

    if (required) {
        stop(
            paste(
                "Could not resolve column; tried:",
                paste(
                    candidates,
                    collapse = ", "
                )
            )
        )
    }

    NA_character_
}


clean_text <- function(
    x
) {

    x <- gsub(
        "[\r\n\t]+",
        " ",
        as.character(x)
    )

    x <- gsub(
        " +",
        " ",
        x
    )

    trimws(x)
}


as_bool <- function(
    x
) {

    if (is.logical(x)) {
        return(x)
    }

    tolower(
        trimws(
            as.character(x)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
    )
}


genotype_frequency <- function(
    x
) {

    x <- as.character(x)

    out <- rep(
        NA_real_,
        length(x)
    )

    out[x == "A/A"] <- 0
    out[x == "P/A"] <- 0.5
    out[x == "P/P"] <- 1

    out
}


extract_sample_id <- function(
    d
) {

    if ("sample" %in% names(d)) {
        return(
            as.character(
                d$sample
            )
        )
    }

    # Previous plot labels begin with the accession/sample ID.
    sub(
        " *\\|.*$",
        "",
        as.character(
            d$sample_label
        )
    )
}


# ============================================================
# METADATA
# ============================================================

if (!file.exists(metadata_file)) {
    stop(
        paste(
            "Missing canonical metadata:",
            metadata_file
        )
    )
}

metadata <- read.delim(
    metadata_file,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

accession_col <- find_col(
    metadata,
    c(
        "Accession",
        "sample",
        "sample_id"
    )
)

taxon_col <- find_col(
    metadata,
    c(
        "Taxon Name",
        "taxon"
    )
)

country_col <- find_col(
    metadata,
    c(
        "Country",
        "country"
    )
)

within_col <- find_col(
    metadata,
    c(
        "Grouping-witihin-species",
        "Grouping-within-species",
        "population",
        "population_group"
    ),
    required = FALSE
)

species_group_col <- find_col(
    metadata,
    c(
        "Species_Group",
        "Species Group"
    ),
    required = FALSE
)


meta <- data.frame(
    sample = clean_text(
        metadata[[accession_col]]
    ),
    taxon = clean_text(
        metadata[[taxon_col]]
    ),
    country = clean_text(
        metadata[[country_col]]
    ),
    stringsAsFactors = FALSE
)

if (!is.na(within_col)) {

    meta$within_group <- clean_text(
        metadata[[within_col]]
    )

} else {

    meta$within_group <- ""
}

if (!is.na(species_group_col)) {

    meta$species_group <- clean_text(
        metadata[[species_group_col]]
    )

} else {

    meta$species_group <- ""
}

meta <- meta[
    nzchar(meta$sample),
    ,
    drop = FALSE
]

meta <- meta[
    !duplicated(meta$sample),
    ,
    drop = FALSE
]


# ============================================================
# POPULATION DEFINITIONS
#
# Population = taxon + country + existing within-species
# grouping where available.
#
# If taxon/country contains >1 within-species group, the group
# number is added to the visible label; otherwise it is hidden.
# ============================================================

meta$base_population <- ifelse(
    nzchar(meta$country),
    paste(
        meta$taxon,
        meta$country,
        sep = " | "
    ),
    meta$taxon
)

meta$population_key <- paste(
    meta$base_population,
    meta$within_group,
    sep = "||"
)

base_group_counts <- aggregate(
    within_group ~ base_population,
    data = unique(
        meta[
            ,
            c(
                "base_population",
                "within_group"
            )
        ]
    ),
    FUN = function(x) {
        length(
            unique(
                x[nzchar(x)]
            )
        )
    }
)

names(base_group_counts)[2] <- "n_within_groups"

meta$n_within_groups <- (
    base_group_counts$n_within_groups[
        match(
            meta$base_population,
            base_group_counts$base_population
        )
    ]
)

meta$population_label <- meta$base_population

use_suffix <- (
    meta$n_within_groups > 1 &
    nzchar(meta$within_group)
)

meta$population_label[
    use_suffix
] <- paste0(
    meta$population_label[
        use_suffix
    ],
    " | pop ",
    sub(
        "\\.0$",
        "",
        meta$within_group[
            use_suffix
        ]
    )
)



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


# ============================================================
# ONE REGION
# ============================================================

prepare_panel <- function(
    panel_name,
    input_file
) {

    cat(
        "\n========================================\n",
        panel_name,
        "\n========================================\n",
        sep = ""
    )

    if (!file.exists(input_file)) {
        stop(
            paste(
                "Missing established genotype-plot input:",
                input_file
            )
        )
    }

    d <- read.delim(
        gzfile(input_file),
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    required <- c(
        "chrom",
        "pos",
        "ref",
        "alt",
        "sample_label",
        "plot_group",
        "display_GT"
    )

    missing_required <- setdiff(
        required,
        names(d)
    )

    if (length(missing_required) > 0) {
        stop(
            paste(
                panel_name,
                "input is missing:",
                paste(
                    missing_required,
                    collapse = ", "
                )
            )
        )
    }

    d$pos <- suppressWarnings(
        as.numeric(
            d$pos
        )
    )

    d$sample_id <- extract_sample_id(
        d
    )

    d$sample_label <- clean_text(
        d$sample_label
    )

    d$plot_group <- clean_text(
        d$plot_group
    )

    d$display_GT <- clean_text(
        d$display_GT
    )


    # --------------------------------------------------------
    # Genomic site order: identical equal-width SNP-column
    # grammar to the existing genotype plots.
    # --------------------------------------------------------

    sites <- unique(
        d[
            order(
                d$pos,
                d$ref,
                d$alt
            ),
            c(
                "chrom",
                "pos",
                "ref",
                "alt"
            )
        ]
    )

    sites$site_key <- paste(
        sites$chrom,
        sites$pos,
        sites$ref,
        sites$alt,
        sep = ":"
    )

    sites$site_index <- seq_len(
        nrow(sites)
    )

    d$site_key <- paste(
        d$chrom,
        d$pos,
        d$ref,
        d$alt,
        sep = ":"
    )

    d$site_index <- sites$site_index[
        match(
            d$site_key,
            sites$site_key
        )
    ]


    # --------------------------------------------------------
    # Existing special rows.
    # --------------------------------------------------------

    is_bwpk <- grepl(
        "Butterfly World piano keys",
        d$plot_group,
        ignore.case = TRUE,
        fixed = TRUE
    )

    is_focal_group <- grepl(
        "Focal pooled samples",
        d$plot_group,
        ignore.case = TRUE,
        fixed = TRUE
    )

    is_left_recip <- (
        grepl(
            "^TN_Left *\\|",
            d$sample_label
        ) &
        grepl(
            "reciprocal",
            d$sample_label,
            ignore.case = TRUE
        )
    )

    is_right_recip <- (
        grepl(
            "^TN_Right *\\|",
            d$sample_label
        ) &
        grepl(
            "reciprocal",
            d$sample_label,
            ignore.case = TRUE
        )
    )


    # --------------------------------------------------------
    # Comparative butterflies only.
    #
    # The previous input already knows exactly which individuals
    # belonged to this region's original comparative cohort.
    # --------------------------------------------------------

    comparative <- d[
        !is_bwpk &
        !is_focal_group &
        !is_left_recip &
        !is_right_recip,
        ,
        drop = FALSE
    ]

    # --------------------------------------------------------
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
    )


    # --------------------------------------------------------
    # Population sample counts.
    # --------------------------------------------------------

    population_manifest <- unique(
        comparative[
            ,
            c(
                "population_key",
                "population_label",
                "species_group",
                "taxon",
                "country",
                "sample_id"
            )
        ]
    )

    population_n <- aggregate(
        sample_id ~
            population_key +
            population_label +
            species_group +
            taxon +
            country,
        data = population_manifest,
        FUN = function(x) {
            length(
                unique(x)
            )
        }
    )

    names(population_n)[
        names(population_n) == "sample_id"
    ] <- "population_n"


    # --------------------------------------------------------
    # Population allele frequency at each site.
    #
    # Mean of:
    #   A/A = 0
    #   P/A = 0.5
    #   P/P = 1
    #
    # N / Uncertain do not contribute to the mean.
    # --------------------------------------------------------

    callable <- comparative[
        is.finite(
            comparative$GT_frequency
        ),
        ,
        drop = FALSE
    ]

    pop_frequency <- aggregate(
        GT_frequency ~
            population_key +
            population_label +
            site_key +
            site_index +
            pos,
        data = callable,
        FUN = mean
    )

    names(pop_frequency)[
        names(pop_frequency) == "GT_frequency"
    ] <- "plot_value"


    # Number callable at each population/site, retained for audit.
    pop_callable <- aggregate(
        sample_id ~
            population_key +
            site_key,
        data = callable,
        FUN = function(x) {
            length(
                unique(x)
            )
        }
    )

    names(pop_callable)[
        names(pop_callable) == "sample_id"
    ] <- "callable_n"


    # --------------------------------------------------------
    # Explicit full population x site grid.
    #
    # This means sites with zero callable individuals are visibly
    # missing rather than silently absent from the plot.
    # --------------------------------------------------------

    pop_grid <- expand.grid(
        population_key = unique(
            population_n$population_key
        ),
        site_key = sites$site_key,
        stringsAsFactors = FALSE
    )

    pop_grid <- merge(
        pop_grid,
        pop_frequency[
            ,
            c(
                "population_key",
                "site_key",
                "plot_value"
            )
        ],
        by = c(
            "population_key",
            "site_key"
        ),
        all.x = TRUE,
        sort = FALSE
    )

    pop_grid <- merge(
        pop_grid,
        pop_callable,
        by = c(
            "population_key",
            "site_key"
        ),
        all.x = TRUE,
        sort = FALSE
    )

    pop_grid <- merge(
        pop_grid,
        population_n,
        by = "population_key",
        all.x = TRUE,
        sort = FALSE
    )

    pop_grid$site_index <- sites$site_index[
        match(
            pop_grid$site_key,
            sites$site_key
        )
    ]

    pop_grid$pos <- sites$pos[
        match(
            pop_grid$site_key,
            sites$site_key
        )
    ]

    pop_grid$display_row <- paste0(
        pop_grid$population_label,
        "  (n=",
        pop_grid$population_n,
        ")"
    )

    pop_grid$display_state <- NA_character_
    pop_grid$plot_group <- "Population frequencies"


    # --------------------------------------------------------
    # Order population rows by established taxonomic context.
    # --------------------------------------------------------

    population_order <- unique(
        population_n[
            order(
                population_n$species_group,
                population_n$taxon,
                population_n$country,
                population_n$population_label
            ),
            c(
                "population_key",
                "population_label",
                "population_n"
            )
        ]
    )

    population_row_order <- paste0(
        population_order$population_label,
        "  (n=",
        population_order$population_n,
        ")"
    )


    # --------------------------------------------------------
    # Final reciprocal L and R only.
    #
    # We deliberately do NOT retain the old original-call row;
    # this plot has one line per focal side.
    # --------------------------------------------------------

    focal <- d[
        is_left_recip |
        is_right_recip,
        ,
        drop = FALSE
    ]

    if (
        !any(is_left_recip) ||
        !any(is_right_recip)
    ) {
        stop(
            paste(
                panel_name,
                ": could not find both reciprocal TN_Left and TN_Right rows."
            )
        )
    }

    focal$display_row <- ifelse(
        grepl(
            "^TN_Left",
            focal$sample_label
        ),
        "TN_Left",
        "TN_Right"
    )

    focal$plot_group <- "Focal pooled samples"

    focal$plot_value <- genotype_frequency(
        focal$display_GT
    )

    focal$display_state <- focal$display_GT

    focal <- focal[
        ,
        c(
            "site_key",
            "site_index",
            "pos",
            "display_row",
            "plot_group",
            "plot_value",
            "display_state"
        )
    ]


    # --------------------------------------------------------
    # Each BWPK individual gets one separate line.
    # --------------------------------------------------------

    bwpk <- d[
        is_bwpk,
        ,
        drop = FALSE
    ]

    if (nrow(bwpk) == 0) {
        stop(
            paste(
                panel_name,
                ": no Butterfly World piano-key rows found."
            )
        )
    }

    bwpk_manifest <- unique(
        bwpk[
            order(
                bwpk$sample_id
            ),
            c(
                "sample_id",
                "sample_label"
            )
        ]
    )

    bwpk$display_row <- bwpk$sample_id

    bwpk$plot_group <- "Butterfly World piano keys"

    bwpk$plot_value <- genotype_frequency(
        bwpk$display_GT
    )

    bwpk$display_state <- bwpk$display_GT

    bwpk <- bwpk[
        ,
        c(
            "site_key",
            "site_index",
            "pos",
            "display_row",
            "plot_group",
            "plot_value",
            "display_state"
        )
    ]


    # ========================================================
    # REMOVE POSITIONS FIXED FOR RAW REFERENCE
    #
    # Scope:
    #   * exact comparative BA cohort for this pattern
    #   * TN_Left reciprocal-mixture row
    #   * TN_Right reciprocal-mixture row
    #   * all Butterfly World piano-key individuals
    #
    # Remove a SNP iff:
    #   * >=1 displayed genotype is callable
    #   * every callable displayed RAW genotype is 0/0
    #
    # N / missing / Uncertain calls are ignored.
    # An all-missing site is retained.
    # ========================================================

    if (!("raw_GT" %in% names(d))) {
        stop(
            paste0(
                panel_name,
                ": raw_GT is absent from the established ",
                "genotype-plot input."
            )
        )
    }


    # --------------------------------------------------------
    # Select precisely the butterflies represented by this
    # figure.
    #
    # allowed_samples was already frozen by Analysis 152 to
    # the exact Band / Ray / Dennis BA cohort.
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
    # Callable means that the established genotype plot itself
    # displayed a resolved genotype.
    # --------------------------------------------------------

    site_calls$callable <- (
        site_calls$display_GT %in%
            c(
                "A/A",
                "P/A",
                "P/P"
            )
    )


    # --------------------------------------------------------
    # Canonicalize RAW REF/ALT genotype.
    # --------------------------------------------------------

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
        !(site_calls$raw_GT_canonical %in%
            c(
                "0/0",
                "0/1",
                "1/1"
            )
        )
    ] <- NA_character_


    site_calls$usable <- (
        site_calls$callable &
        !is.na(
            site_calls$raw_GT_canonical
        )
    )


    # --------------------------------------------------------
    # One fixed-reference classification per physical SNP.
    # --------------------------------------------------------

    calls_by_site <- split(
        site_calls,
        site_calls$site_key
    )

    fixed_ref_rows <- lapply(
        names(calls_by_site),
        function(site_key_value) {

            z <- calls_by_site[[site_key_value]]

            usable_gt <- z$raw_GT_canonical[
                z$usable
            ]

            callable_n <- length(usable_gt)

            ref_hom_n <- sum(
                usable_gt == "0/0"
            )

            het_n <- sum(
                usable_gt == "0/1"
            )

            alt_hom_n <- sum(
                usable_gt == "1/1"
            )

            nonref_n <- (
                het_n +
                alt_hom_n
            )

            data.frame(
                site_key = site_key_value,
                callable_n = callable_n,
                ref_hom_n = ref_hom_n,
                het_n = het_n,
                alt_hom_n = alt_hom_n,
                nonref_n = nonref_n,

                fixed_reference = (
                    callable_n > 0 &&
                    nonref_n == 0
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
    # Join back to the complete original SNP list.
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

    for (count_column in count_columns) {

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

    original_site_n <- nrow(sites)
    removed_site_n <- length(remove_site_keys)
    retained_site_n <- length(keep_site_keys)

    if (retained_site_n < 1L) {
        stop(
            paste0(
                panel_name,
                ": fixed-reference filtering removed every SNP."
            )
        )
    }


    # --------------------------------------------------------
    # Audit.
    # --------------------------------------------------------

    fixed_ref_audit_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            panel_name,
            "_fixed_reference_site_filter.tsv"
        )
    )

    write.table(
        fixed_ref_audit,
        file = fixed_ref_audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    # --------------------------------------------------------
    # Remove fixed-reference columns and re-index retained SNPs
    # consecutively so there are no blank visual gaps.
    # --------------------------------------------------------

    sites <- sites[
        sites$site_key %in% keep_site_keys,
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
        pop_grid$site_key %in% keep_site_keys,
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
        focal$site_key %in% keep_site_keys,
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
        bwpk$site_key %in% keep_site_keys,
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

    pop_plot <- pop_grid[
        ,
        c(
            "site_key",
            "site_index",
            "pos",
            "display_row",
            "plot_group",
            "plot_value",
            "display_state"
        )
    ]

    plot_data <- rbind(
        pop_plot,
        focal,
        bwpk
    )


    # --------------------------------------------------------
    # Stable display order.
    # --------------------------------------------------------

    focal_order <- c(
        "TN_Left",
        "TN_Right"
    )

    bwpk_order <- sort(
        unique(
            bwpk$display_row
        )
    )

    display_order <- c(
        population_row_order,
        focal_order,
        bwpk_order
    )

    plot_data$display_row <- factor(
        plot_data$display_row,
        levels = rev(
            display_order
        ),
        ordered = TRUE
    )

    group_order <- c(
        "Population frequencies",
        "Focal pooled samples",
        "Butterfly World piano keys"
    )

    plot_data$plot_group <- factor(
        plot_data$plot_group,
        levels = group_order,
        ordered = TRUE
    )


    # --------------------------------------------------------
    # Uncertain individual calls get the old dark-grey state.
    # Missing/N remains pale grey through the NA fill.
    # --------------------------------------------------------

    uncertain <- plot_data[
        plot_data$display_state == "Uncertain",
        ,
        drop = FALSE
    ]

    ordinary <- plot_data[
        is.na(plot_data$display_state) |
        plot_data$display_state != "Uncertain",
        ,
        drop = FALSE
    ]


    # --------------------------------------------------------
    # Position ticks, in the same equal-SNP-width grammar as the
    # existing genotype plots.
    # --------------------------------------------------------

    tick_count <- min(
        18L,
        nrow(sites)
    )

    tick_index <- unique(
        round(
            seq(
                1,
                nrow(sites),
                length.out = tick_count
            )
        )
    )

    tick_labels <- format(
        sites$pos[
            tick_index
        ],
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
    )


    # --------------------------------------------------------
    # PLOT
    # --------------------------------------------------------

    p <- ggplot() +

        geom_tile(
            data = ordinary,
            aes(
                x = site_index,
                y = display_row,
                fill = plot_value
            ),
            width = 0.98,
            height = 0.96,
            colour = NA
        ) +

        geom_tile(
            data = uncertain,
            aes(
                x = site_index,
                y = display_row
            ),
            width = 0.98,
            height = 0.96,
            fill = uncertain_colour,
            colour = NA
        ) +

        facet_grid(
            rows = vars(
                plot_group
            ),
            scales = "free_y",
            space = "free_y",
            switch = "y"
        ) +

        scale_fill_gradientn(
            colours = c(
                AA_colour,
                PA_colour,
                PP_colour
            ),
            values = c(
                0,
                0.5,
                1
            ),
            limits = c(
                0,
                1
            ),
            na.value = missing_colour,
            breaks = c(
                0,
                0.5,
                1
            ),
            labels = c(
                "A/A",
                "P/A",
                "P/P"
            ),
            name = paste(
                "Pattern-associated genotype /",
                "population allele frequency"
            )
        ) +

        scale_x_continuous(
            breaks = tick_index,
            labels = tick_labels,
            limits = c(
                0.5,
                nrow(sites) + 0.5
            ),
            expand = expansion(
                add = c(
                    0,
                    0
                )
            )
        ) +

        scale_y_discrete(
            expand = expansion(
                add = c(
                    0.05,
                    0.05
                )
            )
        ) +

        labs(
            title = paste0(
                panel_name,
                ": population frequency and focal genotype context"
            ),

            subtitle = paste(
                format(
                    nrow(sites),
                    big.mark = ","
                ),
                "biallelic SNPs;",
                nrow(population_order),
                "comparative populations;",
                length(bwpk_order),
                "Butterfly World individuals."
            ),

            x = "Genomic position on Hmr18 (bp)",

            y = NULL,

            caption = paste(
                "Population rows contain only butterflies used in the corresponding final comparative BA cohort",
                "and show mean pattern-associated allele dosage among callable individuals: A/A=0, P/A=0.5, P/P=1.",
                "TN_Left and TN_Right use the reciprocal-mixture call only.",
                "Pale grey = unavailable; dark grey = uncertain.",
                "SNP columns are equal width and ordered by genomic position."
            )
        ) +

        guides(
            fill = guide_colourbar(
                title.position = "top",
                title.hjust = 0,
                barwidth = unit(
                    5.0,
                    "cm"
                ),
                barheight = unit(
                    0.34,
                    "cm"
                )
            )
        ) +

        theme_minimal(
            base_size = 9
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 14
            ),

            plot.subtitle = element_text(
                size = 9.5,
                colour = "#444444",
                margin = margin(
                    b = 8
                )
            ),

            plot.caption = element_text(
                size = 8,
                colour = "#444444",
                hjust = 0,
                margin = margin(
                    t = 8
                )
            ),

            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 5.5
            ),

            axis.text.y = element_text(
                size = 6.0,
                colour = "#222222"
            ),

            axis.title.x = element_text(
                size = 10,
                margin = margin(
                    t = 8
                )
            ),

            panel.grid = element_blank(),

            panel.spacing.y = unit(
                0.12,
                "lines"
            ),

            strip.placement = "outside",

            strip.background = element_rect(
                fill = "#E6E6E6",
                colour = NA
            ),

            strip.text.y.left = element_text(
                angle = 0,
                face = "bold",
                size = 8
            ),

            legend.position = "top",

            legend.justification = "left",

            legend.title = element_text(
                face = "bold",
                size = 8.5
            ),

            legend.text = element_text(
                size = 8
            ),

            plot.margin = margin(
                12,
                16,
                12,
                12
            )
        )


    # --------------------------------------------------------
    # Audits
    # --------------------------------------------------------

    population_audit <- pop_grid[
        order(
            pop_grid$population_label,
            pop_grid$site_index
        ),
        c(
            "population_key",
            "population_label",
            "population_n",
            "site_key",
            "site_index",
            "pos",
            "callable_n",
            "plot_value"
        )
    ]

    population_audit_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            panel_name,
            "_population_frequency.tsv.gz"
        )
    )

    gz <- gzfile(
        population_audit_file,
        "wt"
    )

    write.table(
        population_audit,
        file = gz,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    close(gz)


    row_audit_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            panel_name,
            "_row_order.tsv"
        )
    )

    write.table(
        data.frame(
            plot_group = c(
                rep(
                    "Population frequencies",
                    length(population_row_order)
                ),
                rep(
                    "Focal pooled samples",
                    length(focal_order)
                ),
                rep(
                    "Butterfly World piano keys",
                    length(bwpk_order)
                )
            ),
            display_row = display_order,
            stringsAsFactors = FALSE
        ),
        file = row_audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    cat(
        "SNPs: ",
        nrow(sites),
        "\nComparative populations: ",
        nrow(population_order),
        "\nBWPK individuals: ",
        length(bwpk_order),
        "\n",
        sep = ""
    )

    cat(
        "BWPK: ",
        paste(
            bwpk_order,
            collapse = ", "
        ),
        "\n",
        sep = ""
    )

    list(
        panel = panel_name,
        plot = p,
        SNPs = nrow(sites),
        populations = nrow(population_order),
        BWPK = length(bwpk_order),
        rows = length(display_order)
    )
}


# ============================================================
# GENERATE ALL FIVE
# ============================================================

results <- vector(
    "list",
    nrow(panels)
)

for (i in seq_len(nrow(panels))) {

    results[[i]] <- prepare_panel(
        panels$panel[i],
        panels$input[i]
    )

    result <- results[[i]]

    pdf_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            result$panel,
            "_population_frequency_plus_LR_BWPK.pdf"
        )
    )

    png_file <- file.path(
        output_dir,
        paste0(
            "optix_",
            result$panel,
            "_population_frequency_plus_LR_BWPK.png"
        )
    )

    plot_height <- max(
        7.0,
        3.2 + 0.23 * result$rows
    )

    ggsave(
        filename = pdf_file,
        plot = result$plot,
        device = cairo_pdf,
        width = 16,
        height = plot_height,
        units = "in",
        limitsize = FALSE
    )

    ggsave(
        filename = png_file,
        plot = result$plot,
        width = 16,
        height = plot_height,
        units = "in",
        dpi = 300,
        limitsize = FALSE
    )

    cat(
        "PDF: ",
        pdf_file,
        "\n",
        sep = ""
    )
}


# ============================================================
# COMBINED MULTIPAGE PDF
# ============================================================

combined_pdf <- file.path(
    output_dir,
    "optix_Band_Ray_Dennis_population_frequency_plus_LR_BWPK.combined.pdf"
)

grDevices::cairo_pdf(
    filename = combined_pdf,
    width = 16,
    height = 10,
    onefile = TRUE
)

for (result in results) {
    print(
        result$plot
    )
}

dev.off()


# ============================================================
# SUMMARY
# ============================================================

summary <- do.call(
    rbind,
    lapply(
        results,
        function(x) {
            data.frame(
                panel = x$panel,
                SNPs = x$SNPs,
                populations = x$populations,
                BWPK_individuals = x$BWPK,
                displayed_rows = x$rows,
                stringsAsFactors = FALSE
            )
        }
    )
)

summary_file <- file.path(
    output_dir,
    "optix_Band_Ray_Dennis_population_frequency_plus_LR_BWPK.summary.tsv"
)

write.table(
    summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\n========================================\n",
    "ANALYSIS 152 EXACT-BA POPULATION FIGURE COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    summary,
    row.names = FALSE
)

cat(
    "\nCombined PDF:\n",
    combined_pdf,
    "\n\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
