#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_file <- file.path(
    project,
    "scripts",
    "156_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_order.R"
)

target_file <- file.path(
    project,
    "scripts",
    "157_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_labels.R"
)

if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing working Analysis-156 plotter:",
            source_file
        )
    )
}

text <- paste(
    readLines(
        source_file,
        warn = FALSE
    ),
    collapse = "\n"
)


# ============================================================
# NEW OUTPUT DIRECTORY
# ============================================================

old_output <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_order"
)

new_output <- (
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref_taxon_canonical_BDR_labels"
)

hits <- gregexpr(
    old_output,
    text,
    fixed = TRUE
)[[1]]

n_hits <- if (
    length(hits) == 1L &&
    hits[1] == -1L
) {
    0L
} else {
    length(hits)
}

if (n_hits != 1L) {
    stop(
        paste(
            "Expected exactly one Analysis-156 output directory;",
            "found",
            n_hits
        )
    )
}

text <- sub(
    old_output,
    new_output,
    text,
    fixed = TRUE
)


# ============================================================
# REPLACE PLAIN TAXON LABEL WITH B/D/R PHENOTYPE TAG
#
# + = phenotype present throughout canonical taxon
# - = phenotype absent throughout canonical taxon
# ~ = mixed after taxon/country collapsing
# ? = unavailable
#
# This changes LABELS ONLY.
# ============================================================

old_label_block <- paste(
'    comparative$population_label <- (',
'        cohort_manifest$population[',
'            cohort_index',
'        ]',
'    )',
    sep = "\n"
)

new_label_block <- paste(
'    # --------------------------------------------------------',
'    # Compact B/D/R phenotype annotation.',
'    #',
'    # This is the taxon-level analogue of the phenotype strip',
'    # used in the earlier individual genotype plots.',
'    # --------------------------------------------------------',
'',
'    taxon_B_fraction <- suppressWarnings(',
'        as.numeric(',
'            cohort_manifest$B_fraction_taxon[',
'                cohort_index',
'            ]',
'        )',
'    )',
'',
'    taxon_D_fraction <- suppressWarnings(',
'        as.numeric(',
'            cohort_manifest$D_fraction_taxon[',
'                cohort_index',
'            ]',
'        )',
'    )',
'',
'    taxon_R_fraction <- suppressWarnings(',
'        as.numeric(',
'            cohort_manifest$R_fraction_taxon[',
'                cohort_index',
'            ]',
'        )',
'    )',
'',
'    phenotype_tag <- function(x) {',
'',
'        result <- rep(',
'            "?",',
'            length(x)',
'        )',
'',
'        result[',
'            is.finite(x) &',
'            abs(x) < 1e-8',
'        ] <- "-"',
'',
'        result[',
'            is.finite(x) &',
'            abs(x - 1) < 1e-8',
'        ] <- "+"',
'',
'        result[',
'            is.finite(x) &',
'            x > 1e-8 &',
'            x < (1 - 1e-8)',
'        ] <- "~"',
'',
'        result',
'    }',
'',
'    taxon_B_tag <- phenotype_tag(',
'        taxon_B_fraction',
'    )',
'',
'    taxon_D_tag <- phenotype_tag(',
'        taxon_D_fraction',
'    )',
'',
'    taxon_R_tag <- phenotype_tag(',
'        taxon_R_fraction',
'    )',
'',
'    comparative$population_label <- paste0(',
'        "B",',
'        taxon_B_tag,',
'        "  D",',
'        taxon_D_tag,',
'        "  R",',
'        taxon_R_tag,',
'        "    ",',
'        cohort_manifest$population[',
'            cohort_index',
'        ]',
'    )',
    sep = "\n"
)

hits <- gregexpr(
    old_label_block,
    text,
    fixed = TRUE
)[[1]]

n_hits <- if (
    length(hits) == 1L &&
    hits[1] == -1L
) {
    0L
} else {
    length(hits)
}

if (n_hits != 1L) {
    stop(
        paste(
            "Expected exactly one plain population-label block;",
            "found",
            n_hits
        )
    )
}

text <- sub(
    old_label_block,
    new_label_block,
    text,
    fixed = TRUE
)


# ============================================================
# WRITE DERIVATIVE
# ============================================================

writeLines(
    text,
    con = target_file,
    useBytes = TRUE
)

cat(
    "Created:\n",
    target_file,
    "\n\n",
    sep = ""
)

cat(
    "Taxon labels now use:\n",
    "  B+ D- R- = present/absent/absent\n",
    "  ~         = mixed within canonical taxon\n",
    "  ?         = unavailable\n",
    sep = ""
)
