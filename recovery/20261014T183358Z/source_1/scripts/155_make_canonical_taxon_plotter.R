#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_plotter <- file.path(
    project,
    "scripts",
    "154_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_only.R"
)

source_manifest <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.tsv"
)

canonical_manifest <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.taxon_canonical.tsv"
)

canonical_summary <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_taxon_canonical_summary.tsv"
)

canonical_audit <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_taxon_canonicalization_audit.tsv"
)

target_plotter <- file.path(
    project,
    "scripts",
    "155_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical.R"
)


# ============================================================
# CHECK INPUTS
# ============================================================

if (!file.exists(source_plotter)) {
    stop(
        paste(
            "Missing working Analysis-154 plotter:",
            source_plotter
        )
    )
}

if (!file.exists(source_manifest)) {
    stop(
        paste(
            "Missing validated exact-BA manifest:",
            source_manifest
        )
    )
}


# ============================================================
# CANONICAL TAXON-NAME FUNCTION
#
# Intentional normalizations:
#
#   ssp. 'meriana'  -> ssp. meriana
#   spp. melpomene  -> ssp. melpomene
#   spp. tristero   -> ssp. tristero
#   nov. 'vicina'   -> nov. vicina
#
# Also collapse hidden newlines / repeated spaces.
# ============================================================

canonicalize_taxon <- function(x) {

    x <- as.character(x)

    x <- gsub(
        "[\r\n\t]+",
        " ",
        x,
        perl = TRUE
    )

    x <- trimws(
        x
    )

    x <- gsub(
        "\\s+",
        " ",
        x,
        perl = TRUE
    )

    # In this metadata, spp. is being used inconsistently
    # with ssp. for the same taxonomic designation.
    x <- gsub(
        "\\bspp\\.",
        "ssp.",
        x,
        perl = TRUE
    )

    # Remove single quotation marks around provisional names.
    x <- gsub(
        "'",
        "",
        x,
        fixed = TRUE
    )

    x <- gsub(
        "\\s+",
        " ",
        x,
        perl = TRUE
    )

    trimws(
        x
    )
}


# ============================================================
# READ ORIGINAL VALIDATED BA MANIFEST
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

required <- c(
    "BA_pattern",
    "sample",
    "taxon",
    "country",
    "data_type"
)

missing_required <- setdiff(
    required,
    names(m)
)

if (length(missing_required) > 0) {
    stop(
        paste(
            "Manifest missing required columns:",
            paste(
                missing_required,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# PRESERVE ORIGINAL VALUES FOR AUDIT
# ============================================================

m$taxon_original <- as.character(
    m$taxon
)

m$country_original <- as.character(
    m$country
)

m$taxon <- canonicalize_taxon(
    m$taxon
)


# ============================================================
# TAXON IS NOW THE ONLY GROUPING DESIGNATION
#
# population is what the existing plotter uses as its main
# population/taxon key.
#
# country is deliberately made constant so it cannot continue
# to split identical taxa internally during aggregate().
# ============================================================

m$population <- m$taxon

m$country <- "taxon_only"


# ============================================================
# HARD VALIDATION: SAMPLE MEMBERSHIP MUST NOT CHANGE
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
            "Exact BA membership changed unexpectedly:",
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
# WRITE CANONICAL MANIFEST
# ============================================================

write.table(
    m,
    file = canonical_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# CANONICALIZATION AUDIT
# ============================================================

audit <- unique(
    m[
        ,
        c(
            "BA_pattern",
            "taxon_original",
            "taxon"
        ),
        drop = FALSE
    ]
)

audit$n_samples <- vapply(
    seq_len(
        nrow(audit)
    ),
    function(i) {

        length(
            unique(
                m$sample[
                    m$BA_pattern ==
                        audit$BA_pattern[i]
                    &
                    m$taxon_original ==
                        audit$taxon_original[i]
                ]
            )
        )
    },
    integer(1)
)

audit <- audit[
    order(
        audit$BA_pattern,
        audit$taxon,
        audit$taxon_original
    ),
    ,
    drop = FALSE
]

write.table(
    audit,
    file = canonical_audit,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# TAXON SUMMARY
# ============================================================

sample_taxon <- unique(
    m[
        ,
        c(
            "BA_pattern",
            "sample",
            "taxon",
            "data_type"
        ),
        drop = FALSE
    ]
)

summary_n <- aggregate(
    sample ~ BA_pattern + taxon,
    data = sample_taxon,
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(summary_n)[
    names(summary_n) == "sample"
] <- "samples"

summary_wgs <- aggregate(
    sample ~ BA_pattern + taxon,
    data = sample_taxon[
        sample_taxon$data_type == "WGS",
        ,
        drop = FALSE
    ],
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(summary_wgs)[
    names(summary_wgs) == "sample"
] <- "WGS"

summary_cs <- aggregate(
    sample ~ BA_pattern + taxon,
    data = sample_taxon[
        sample_taxon$data_type == "CS",
        ,
        drop = FALSE
    ],
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(summary_cs)[
    names(summary_cs) == "sample"
] <- "CS"

taxon_summary <- merge(
    summary_n,
    summary_wgs,
    by = c(
        "BA_pattern",
        "taxon"
    ),
    all.x = TRUE,
    sort = FALSE
)

taxon_summary <- merge(
    taxon_summary,
    summary_cs,
    by = c(
        "BA_pattern",
        "taxon"
    ),
    all.x = TRUE,
    sort = FALSE
)

taxon_summary$WGS[
    is.na(
        taxon_summary$WGS
    )
] <- 0L

taxon_summary$CS[
    is.na(
        taxon_summary$CS
    )
] <- 0L

taxon_summary <- taxon_summary[
    order(
        taxon_summary$BA_pattern,
        taxon_summary$taxon
    ),
    ,
    drop = FALSE
]

write.table(
    taxon_summary,
    file = canonical_summary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# DERIVE ANALYSIS 155 FROM WORKING 154
# ============================================================

text <- paste(
    readLines(
        source_plotter,
        warn = FALSE
    ),
    collapse = "\n"
)


old_manifest_name <- (
    "exact_BA_comparative_sample_manifest.taxon_only.tsv"
)

new_manifest_name <- (
    "exact_BA_comparative_sample_manifest.taxon_canonical.tsv"
)

if (
    !grepl(
        old_manifest_name,
        text,
        fixed = TRUE
    )
) {
    stop(
        "Could not find taxon-only manifest name in Analysis 154."
    )
}

text <- sub(
    old_manifest_name,
    new_manifest_name,
    text,
    fixed = TRUE
)


old_output_dir <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_only"
)

new_output_dir <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical"
)

if (
    !grepl(
        old_output_dir,
        text,
        fixed = TRUE
    )
) {
    stop(
        "Could not find Analysis-154 output directory."
    )
}

text <- sub(
    old_output_dir,
    new_output_dir,
    text,
    fixed = TRUE
)


text <- gsub(
    "ANALYSIS 154 NON-FIXED-REF EXACT-BA TAXON FIGURE COMPLETE",
    "ANALYSIS 155 CANONICAL-TAXON NON-FIXED-REF EXACT-BA FIGURE COMPLETE",
    text,
    fixed = TRUE
)


writeLines(
    text,
    con = target_plotter,
    useBytes = TRUE
)


# ============================================================
# REPORT THE COLLAPSE
# ============================================================

cat(
    "\nCANONICAL TAXON GROUPING\n"
)

for (
    pattern_name
    in names(expected_n)
) {

    before_taxa <- length(
        unique(
            m$taxon_original[
                m$BA_pattern ==
                    pattern_name
            ]
        )
    )

    after_taxa <- length(
        unique(
            m$taxon[
                m$BA_pattern ==
                    pattern_name
            ]
        )
    )

    cat(
        sprintf(
            "  %-6s %3d butterflies | %2d raw labels -> %2d canonical taxa\n",
            pattern_name,
            observed_n[
                pattern_name
            ],
            before_taxa,
            after_taxa
        )
    )
}


cat(
    "\nLABELS THAT WERE NORMALIZED\n"
)

changed <- audit[
    audit$taxon_original !=
        audit$taxon,
    ,
    drop = FALSE
]

if (nrow(changed) == 0) {

    cat(
        "  none\n"
    )

} else {

    print(
        changed[
            ,
            c(
                "BA_pattern",
                "taxon_original",
                "taxon",
                "n_samples"
            ),
            drop = FALSE
        ],
        row.names = FALSE
    )
}


cat(
    "\nCreated canonical manifest:\n",
    canonical_manifest,
    "\n\nCreated plotter:\n",
    target_plotter,
    "\n",
    sep = ""
)
