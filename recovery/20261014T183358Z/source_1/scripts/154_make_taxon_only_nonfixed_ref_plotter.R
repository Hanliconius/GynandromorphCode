#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_plotter <- file.path(
    project,
    "scripts",
    "153_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref.R"
)

target_plotter <- file.path(
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

taxon_manifest <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_comparative_sample_manifest.taxon_only.tsv"
)

taxon_summary_file <- file.path(
    project,
    "analysis",
    "figures",
    "key_regions_population_context_genoplots",
    "exact_BA_taxon_summary.tsv"
)


# ============================================================
# INPUT CHECKS
# ============================================================

if (!file.exists(source_plotter)) {
    stop(
        paste(
            "Missing working non-fixed-reference plotter:",
            source_plotter
        )
    )
}

if (!file.exists(source_manifest)) {
    stop(
        paste(
            "Missing validated BA manifest:",
            source_manifest
        )
    )
}


# ============================================================
# REBUILD ONLY THE GROUPING VARIABLE
#
# Sample membership does NOT change.
# Population becomes taxon only.
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
            "Manifest missing columns:",
            paste(
                missing_required,
                collapse = ", "
            )
        )
    )
}

for (column in required) {
    m[[column]] <- trimws(
        as.character(
            m[[column]]
        )
    )
}


# Hard validation: exact BA cohorts must remain unchanged.
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
            "BA cohort count changed unexpectedly:",
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
# TAXON-ONLY GROUPING
#
# This is the ONLY analytical change.
# ============================================================

m$population <- m$taxon

write.table(
    m,
    file = taxon_manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# TAXON SUMMARY
# ============================================================

taxon_sample_manifest <- unique(
    m[
        ,
        c(
            "BA_pattern",
            "population",
            "taxon",
            "sample",
            "data_type"
        ),
        drop = FALSE
    ]
)

taxon_n <- aggregate(
    sample ~ BA_pattern + population + taxon,
    data = taxon_sample_manifest,
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(taxon_n)[
    names(taxon_n) == "sample"
] <- "samples"


wgs_n <- aggregate(
    sample ~ BA_pattern + population + taxon,
    data = taxon_sample_manifest[
        taxon_sample_manifest$data_type == "WGS",
        ,
        drop = FALSE
    ],
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(wgs_n)[
    names(wgs_n) == "sample"
] <- "WGS"


cs_n <- aggregate(
    sample ~ BA_pattern + population + taxon,
    data = taxon_sample_manifest[
        taxon_sample_manifest$data_type == "CS",
        ,
        drop = FALSE
    ],
    FUN = function(x) {
        length(
            unique(x)
        )
    }
)

names(cs_n)[
    names(cs_n) == "sample"
] <- "CS"


taxon_summary <- merge(
    taxon_n,
    wgs_n,
    by = c(
        "BA_pattern",
        "population",
        "taxon"
    ),
    all.x = TRUE,
    sort = FALSE
)

taxon_summary <- merge(
    taxon_summary,
    cs_n,
    by = c(
        "BA_pattern",
        "population",
        "taxon"
    ),
    all.x = TRUE,
    sort = FALSE
)

taxon_summary$WGS[
    is.na(taxon_summary$WGS)
] <- 0L

taxon_summary$CS[
    is.na(taxon_summary$CS)
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
    file = taxon_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# DERIVE PLOTTER FROM WORKING 153
# ============================================================

text <- paste(
    readLines(
        source_plotter,
        warn = FALSE
    ),
    collapse = "\n"
)


old_manifest <- "exact_BA_comparative_sample_manifest.tsv"
new_manifest <- "exact_BA_comparative_sample_manifest.taxon_only.tsv"

if (!grepl(old_manifest, text, fixed = TRUE)) {
    stop(
        "Could not locate exact BA manifest filename in 153."
    )
}

text <- sub(
    old_manifest,
    new_manifest,
    text,
    fixed = TRUE
)


old_output <- "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref"

new_output <- "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_only"

if (!grepl(old_output, text, fixed = TRUE)) {
    stop(
        "Could not locate 153 output-directory name."
    )
}

text <- sub(
    old_output,
    new_output,
    text,
    fixed = TRUE
)


# ============================================================
# VISIBLE LABEL CHANGES ONLY
# ============================================================

text <- gsub(
    "comparative populations",
    "comparative taxa",
    text,
    fixed = TRUE
)

text <- gsub(
    "Comparative populations",
    "Comparative taxa",
    text,
    fixed = TRUE
)

text <- gsub(
    "Population rows",
    "Taxon rows",
    text,
    fixed = TRUE
)

text <- gsub(
    "population allele frequency",
    "taxon allele frequency",
    text,
    fixed = TRUE
)

text <- gsub(
    "Population frequencies",
    "Taxon frequencies",
    text,
    fixed = TRUE
)

text <- gsub(
    "ANALYSIS 153 NON-FIXED-REF EXACT-BA POPULATION FIGURE COMPLETE",
    "ANALYSIS 154 NON-FIXED-REF EXACT-BA TAXON FIGURE COMPLETE",
    text,
    fixed = TRUE
)


writeLines(
    text,
    con = target_plotter,
    useBytes = TRUE
)


# ============================================================
# REPORT
# ============================================================

cat(
    "\nEXACT BA COHORTS RETAINED\n"
)

for (pattern_name in names(expected_n)) {

    n_taxa <- length(
        unique(
            m$taxon[
                m$BA_pattern ==
                    pattern_name
            ]
        )
    )

    cat(
        sprintf(
            "  %-6s %3d butterflies -> %2d taxonomic rows\n",
            pattern_name,
            observed_n[
                pattern_name
            ],
            n_taxa
        )
    )
}

cat(
    "\nTaxon-only manifest:\n",
    taxon_manifest,
    "\n\nTaxon summary:\n",
    taxon_summary_file,
    "\n\nCreated plotter:\n",
    target_plotter,
    "\n",
    sep = ""
)
