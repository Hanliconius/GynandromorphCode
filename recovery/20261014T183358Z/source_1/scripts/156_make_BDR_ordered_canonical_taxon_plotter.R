#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_plotter <- file.path(
    project,
    "scripts",
    "155_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical.R"
)

source_manifest <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.taxon_canonical.tsv"
)

metadata_file <- file.path(
    project,
    "data",
    "derived",
    "region_panel_metadata.tsv"
)

ordered_manifest <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.taxon_canonical_BDR_order.tsv"
)

phenotype_profile_file <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_canonical_taxon_BDR_profile.tsv"
)

target_plotter <- file.path(
    project,
    "scripts",
    "156_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_order.R"
)


# ============================================================
# INPUT CHECKS
# ============================================================

for (
    path
    in c(
        source_plotter,
        source_manifest,
        metadata_file
    )
) {
    if (!file.exists(path)) {
        stop(
            paste(
                "Missing required input:",
                path
            )
        )
    }
}


# ============================================================
# HELPERS
# ============================================================

clean_name <- function(x) {

    x <- as.character(x)

    x <- gsub(
        "[\r\n\t]+",
        " ",
        x,
        perl = TRUE
    )

    x <- gsub(
        "\\s+",
        " ",
        x,
        perl = TRUE
    )

    trimws(x)
}


# Metadata headers in this project can contain punctuation,
# whitespace and embedded line breaks. Match columns using the
# established normalize-to-alphanumeric strategy.
normalize_column_name <- function(x) {

    x <- tolower(
        as.character(x)
    )

    gsub(
        "[^a-z0-9]+",
        "",
        x,
        perl = TRUE
    )
}


resolve_column <- function(
    x,
    candidates
) {

    normalized_columns <- normalize_column_name(
        names(x)
    )

    for (
        candidate
        in candidates
    ) {

        target <- normalize_column_name(
            candidate
        )

        hit <- which(
            normalized_columns ==
                target
        )

        if (length(hit) == 1L) {

            cat(
                sprintf(
                    "Resolved metadata column %-10s -> %s\n",
                    candidate,
                    names(x)[hit]
                )
            )

            return(
                names(x)[hit]
            )
        }
    }

    stop(
        paste0(
            "Could not resolve metadata column from: ",
            paste(
                candidates,
                collapse = ", "
            ),
            "\nAvailable normalized columns:\n",
            paste(
                paste0(
                    names(x),
                    " -> ",
                    normalized_columns
                ),
                collapse = "\n"
            )
        )
    )
}


binary_numeric <- function(x) {

    y <- suppressWarnings(
        as.numeric(
            as.character(x)
        )
    )

    result <- rep(
        NA_real_,
        length(y)
    )

    result[
        is.finite(y) &
        abs(y) < 1e-8
    ] <- 0

    result[
        is.finite(y) &
        abs(y - 1) < 1e-8
    ] <- 1

    result
}


trait_fraction <- function(x) {

    x <- x[
        is.finite(x)
    ]

    if (length(x) == 0L) {
        return(
            NA_real_
        )
    }

    mean(x)
}


observed_states <- function(
    B,
    D,
    R
) {

    keep <- (
        is.finite(B) &
        is.finite(D) &
        is.finite(R)
    )

    if (!any(keep)) {
        return(
            "NA"
        )
    }

    states <- unique(
        paste0(
            "B",
            as.integer(B[keep]),
            "D",
            as.integer(D[keep]),
            "R",
            as.integer(R[keep])
        )
    )

    paste(
        sort(states),
        collapse = ","
    )
}


# ============================================================
# READ THE WORKING CANONICAL-TAXON MANIFEST
#
# Sample membership and taxon collapsing are already frozen.
# We are ONLY adding phenotype-order information.
# ============================================================

m <- read.delim(
    source_manifest,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

required_manifest <- c(
    "BA_pattern",
    "sample",
    "taxon",
    "population",
    "data_type"
)

missing_manifest <- setdiff(
    required_manifest,
    names(m)
)

if (length(missing_manifest) > 0L) {
    stop(
        paste(
            "Canonical manifest missing:",
            paste(
                missing_manifest,
                collapse = ", "
            )
        )
    )
}

m$BA_pattern <- clean_name(
    m$BA_pattern
)

m$sample <- clean_name(
    m$sample
)

m$taxon <- clean_name(
    m$taxon
)


# ============================================================
# HARD CHECK: EXACT BA COHORTS HAVE NOT CHANGED
# ============================================================

expected_n <- c(
    Band = 163L,
    Ray = 368L,
    Dennis = 369L
)

observed_n <- vapply(
    names(expected_n),
    function(pattern_name) {

        length(
            unique(
                m$sample[
                    m$BA_pattern ==
                        pattern_name
                ]
            )
        )
    },
    integer(1)
)

if (
    any(
        observed_n !=
            expected_n
    )
) {
    stop(
        paste(
            "Exact BA cohort changed:",
            paste(
                names(observed_n),
                observed_n,
                sep = "=",
                collapse = ", "
            )
        )
    )
}


# ============================================================
# READ B / D / R PHENOTYPE METADATA
# ============================================================

meta <- read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
)

accession_column <- resolve_column(
    meta,
    c(
        "Accession"
    )
)

B_column <- resolve_column(
    meta,
    c(
        "B"
    )
)

D_column <- resolve_column(
    meta,
    c(
        "D"
    )
)

R_column <- resolve_column(
    meta,
    c(
        "R"
    )
)

metadata_traits <- data.frame(
    sample = clean_name(
        meta[[accession_column]]
    ),

    B = binary_numeric(
        meta[[B_column]]
    ),

    D = binary_numeric(
        meta[[D_column]]
    ),

    R = binary_numeric(
        meta[[R_column]]
    ),

    stringsAsFactors = FALSE
)

metadata_traits <- metadata_traits[
    !duplicated(
        metadata_traits$sample
    ),
    ,
    drop = FALSE
]


# ============================================================
# MAP PHENOTYPES TO THE EXACT BA SAMPLES
# ============================================================

trait_index <- match(
    m$sample,
    metadata_traits$sample
)

if (anyNA(trait_index)) {

    missing_samples <- unique(
        m$sample[
            is.na(trait_index)
        ]
    )

    stop(
        paste0(
            "Exact BA samples absent from phenotype metadata: ",
            paste(
                head(
                    missing_samples,
                    30
                ),
                collapse = ", "
            )
        )
    )
}

m$B <- metadata_traits$B[
    trait_index
]

m$D <- metadata_traits$D[
    trait_index
]

m$R <- metadata_traits$R[
    trait_index
]


# ============================================================
# APPLY THE SAME LOCKED PHENOTYPE CORRECTIONS USED IN THE
# FINAL BAND ANALYSIS.
#
# These affect phenotype ordering only; they do not alter the
# already-frozen BA sample membership or genotype calls.
# ============================================================

taxon_lower <- tolower(
    m$taxon
)

malleti_or_ecuadoriensis <- (
    grepl(
        "melpomene ssp\\. malleti$",
        taxon_lower
    ) |
    grepl(
        "melpomene ssp\\. ecuadoriensis$",
        taxon_lower
    )
)

m$B[
    malleti_or_ecuadoriensis
] <- 0


locked_sample <- (
    m$sample ==
        "ERR3652402"
)

m$B[
    locked_sample
] <- 0

m$D[
    locked_sample
] <- 1

m$R[
    locked_sample
] <- 1


# ============================================================
# ONE B/D/R PROFILE PER CANONICAL TAXON WITHIN EACH BA COHORT
#
# Country is NOT involved.
#
# A canonical taxon remains exactly one plotted row.
# ============================================================

profile_keys <- unique(
    m[
        ,
        c(
            "BA_pattern",
            "taxon"
        ),
        drop = FALSE
    ]
)

profile_rows <- lapply(
    seq_len(
        nrow(profile_keys)
    ),
    function(i) {

        pattern_name <- (
            profile_keys$BA_pattern[i]
        )

        taxon_name <- (
            profile_keys$taxon[i]
        )

        z <- m[
            m$BA_pattern ==
                pattern_name &
            m$taxon ==
                taxon_name,
            ,
            drop = FALSE
        ]

        B_fraction <- trait_fraction(
            z$B
        )

        D_fraction <- trait_fraction(
            z$D
        )

        R_fraction <- trait_fraction(
            z$R
        )

        data.frame(
            BA_pattern = pattern_name,
            taxon = taxon_name,

            butterflies = length(
                unique(
                    z$sample
                )
            ),

            B_fraction = B_fraction,
            D_fraction = D_fraction,
            R_fraction = R_fraction,

            observed_BDR_states = observed_states(
                z$B,
                z$D,
                z$R
            ),

            mixed_BDR = (
                grepl(
                    ",",
                    observed_states(
                        z$B,
                        z$D,
                        z$R
                    ),
                    fixed = TRUE
                )
            ),

            stringsAsFactors = FALSE
        )
    }
)

profile <- do.call(
    rbind,
    profile_rows
)


# ============================================================
# PHENOTYPE-FIRST ROW ORDER
#
# User-requested priority:
#
#   1. B
#   2. D
#   3. R
#   4. taxonomic name as tie-breaker
#
# Fractions allow a country-merged taxon to stay as ONE row
# even if its constituent butterflies are not perfectly
# phenotype-identical.
# ============================================================

profile$B_sort <- ifelse(
    is.finite(
        profile$B_fraction
    ),
    profile$B_fraction,
    -1
)

profile$D_sort <- ifelse(
    is.finite(
        profile$D_fraction
    ),
    profile$D_fraction,
    -1
)

profile$R_sort <- ifelse(
    is.finite(
        profile$R_fraction
    ),
    profile$R_fraction,
    -1
)

profile$BDR_sort <- NA_character_

for (
    pattern_name
    in c(
        "Band",
        "Ray",
        "Dennis"
    )
) {

    indices <- which(
        profile$BA_pattern ==
            pattern_name
    )

    ordered_indices <- indices[
        order(
            -profile$B_sort[
                indices
            ],
            -profile$D_sort[
                indices
            ],
            -profile$R_sort[
                indices
            ],
            profile$taxon[
                indices
            ]
        )
    ]

    profile$BDR_sort[
        ordered_indices
    ] <- sprintf(
        "%03d",
        seq_along(
            ordered_indices
        )
    )
}


if (anyNA(profile$BDR_sort)) {
    stop(
        "Some canonical taxa did not receive a B/D/R row order."
    )
}


# Human-readable block label for the audit.
profile$BDR_profile <- ifelse(
    profile$observed_BDR_states == "B1D0R0",
    "B",
    ifelse(
        profile$observed_BDR_states == "B0D1R0",
        "D",
        ifelse(
            profile$observed_BDR_states == "B0D0R1",
            "R",
            ifelse(
                profile$observed_BDR_states == "B0D1R1",
                "D+R",
                ifelse(
                    profile$observed_BDR_states == "B0D0R0",
                    "none",
                    profile$observed_BDR_states
                )
            )
        )
    )
)


profile <- profile[
    order(
        match(
            profile$BA_pattern,
            c(
                "Band",
                "Ray",
                "Dennis"
            )
        ),
        profile$BDR_sort
    ),
    ,
    drop = FALSE
]


write.table(
    profile[
        ,
        c(
            "BA_pattern",
            "BDR_sort",
            "BDR_profile",
            "taxon",
            "butterflies",
            "B_fraction",
            "D_fraction",
            "R_fraction",
            "observed_BDR_states",
            "mixed_BDR"
        ),
        drop = FALSE
    ],
    file = phenotype_profile_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# ADD THE TAXON-LEVEL SORT CODE TO EVERY SAMPLE IN THE
# CANONICAL MANIFEST.
#
# Every sample from the same taxon gets the SAME BDR_sort,
# ensuring we do NOT recreate duplicate taxon rows.
# ============================================================

manifest_key <- paste(
    m$BA_pattern,
    m$taxon,
    sep = "\r"
)

profile_key <- paste(
    profile$BA_pattern,
    profile$taxon,
    sep = "\r"
)

profile_index <- match(
    manifest_key,
    profile_key
)

if (anyNA(profile_index)) {
    stop(
        "Failed to map some exact BA samples to canonical taxon phenotype profiles."
    )
}

m$BDR_sort <- profile$BDR_sort[
    profile_index
]

m$BDR_profile <- profile$BDR_profile[
    profile_index
]

m$B_fraction_taxon <- profile$B_fraction[
    profile_index
]

m$D_fraction_taxon <- profile$D_fraction[
    profile_index
]

m$R_fraction_taxon <- profile$R_fraction[
    profile_index
]


write.table(
    m,
    file = ordered_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# DERIVE ANALYSIS 156 FROM THE WORKING ANALYSIS 155
# ============================================================

text <- paste(
    readLines(
        source_plotter,
        warn = FALSE
    ),
    collapse = "\n"
)


# ------------------------------------------------------------
# New B/D/R-ordered manifest.
# ------------------------------------------------------------

old_manifest_name <- (
    "exact_BA_comparative_sample_manifest.taxon_canonical.tsv"
)

new_manifest_name <- (
    "exact_BA_comparative_sample_manifest.taxon_canonical_BDR_order.tsv"
)

manifest_hits <- gregexpr(
    old_manifest_name,
    text,
    fixed = TRUE
)[[1]]

manifest_n <- if (
    length(manifest_hits) == 1L &&
    manifest_hits[1] == -1L
) {
    0L
} else {
    length(manifest_hits)
}

if (manifest_n != 1L) {
    stop(
        paste(
            "Expected one canonical manifest filename in 155;",
            "found",
            manifest_n
        )
    )
}

text <- sub(
    old_manifest_name,
    new_manifest_name,
    text,
    fixed = TRUE
)


# ------------------------------------------------------------
# New output directory.
# ------------------------------------------------------------

old_output <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical"
)

new_output <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_order"
)

output_hits <- gregexpr(
    old_output,
    text,
    fixed = TRUE
)[[1]]

output_n <- if (
    length(output_hits) == 1L &&
    output_hits[1] == -1L
) {
    0L
} else {
    length(output_hits)
}

if (output_n != 1L) {
    stop(
        paste(
            "Expected one Analysis-155 output directory;",
            "found",
            output_n
        )
    )
}

text <- sub(
    old_output,
    new_output,
    text,
    fixed = TRUE
)


# ------------------------------------------------------------
# THIS is the actual plot-order change.
#
# Analysis 152 deliberately blanked species_group. The original
# plotter subsequently sorts population/taxon rows using this
# field before taxon name.
#
# We repurpose that hidden sort field for the B/D/R taxon rank.
# It never appears in the visible taxon label.
# ------------------------------------------------------------

old_species_group <- (
    '    comparative$species_group <- ""'
)

new_species_group <- paste(
    '    comparative$species_group <- (',
    '        cohort_manifest$BDR_sort[',
    '            cohort_index',
    '        ]',
    '    )',
    sep = "\n"
)

species_hits <- gregexpr(
    old_species_group,
    text,
    fixed = TRUE
)[[1]]

species_n <- if (
    length(species_hits) == 1L &&
    species_hits[1] == -1L
) {
    0L
} else {
    length(species_hits)
}

if (species_n != 1L) {
    stop(
        paste(
            "Expected one blank species_group assignment;",
            "found",
            species_n
        )
    )
}

text <- sub(
    old_species_group,
    new_species_group,
    text,
    fixed = TRUE
)


# ------------------------------------------------------------
# Version label, if present.
# ------------------------------------------------------------

text <- gsub(
    "ANALYSIS 155 CANONICAL-TAXON NON-FIXED-REF EXACT-BA FIGURE COMPLETE",
    "ANALYSIS 156 BDR-ORDERED CANONICAL-TAXON NON-FIXED-REF EXACT-BA FIGURE COMPLETE",
    text,
    fixed = TRUE
)


writeLines(
    text,
    con = target_plotter,
    useBytes = TRUE
)


# ============================================================
# REPORT EXACT ROW ORDER BEFORE PLOTTING
# ============================================================

cat(
    "\nCANONICAL TAXON ROW ORDER BY B -> D -> R\n"
)

for (
    pattern_name
    in c(
        "Band",
        "Ray",
        "Dennis"
    )
) {

    cat(
        "\n",
        pattern_name,
        "\n",
        sep = ""
    )

    print(
        profile[
            profile$BA_pattern ==
                pattern_name,
            c(
                "BDR_sort",
                "BDR_profile",
                "taxon",
                "butterflies",
                "B_fraction",
                "D_fraction",
                "R_fraction",
                "observed_BDR_states"
            ),
            drop = FALSE
        ],
        row.names = FALSE
    )
}


cat(
    "\nExact BA cohort sizes unchanged:\n"
)

for (
    pattern_name
    in names(expected_n)
) {

    cat(
        sprintf(
            "  %-6s %d\n",
            pattern_name,
            observed_n[
                pattern_name
            ]
        )
    )
}


cat(
    "\nPhenotype profile:\n",
    phenotype_profile_file,
    "\n\nOrdered manifest:\n",
    ordered_manifest,
    "\n\nCreated plotter:\n",
    target_plotter,
    "\n",
    sep = ""
)
