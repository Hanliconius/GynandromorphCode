#!/usr/bin/env Rscript

target <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "157_plot_population_frequency_plus_LR_BWPK_exact_BA_",
    "nonfixed_ref_taxon_canonical_BDR_labels.R"
)

if (!file.exists(target)) {
    stop(
        paste(
            "Missing target:",
            target
        )
    )
}

text <- paste(
    readLines(
        target,
        warn = FALSE
    ),
    collapse = "\n"
)

old <- paste(
'    observed_populations <- sort(',
'        unique(',
'            comparative$population_label',
'        )',
'    )',
    sep = "\n"
)

new <- paste(
'    # Validate canonical taxon identity, not the decorated',
'    # B/D/R display label.',
'    observed_populations <- sort(',
'        unique(',
'            comparative$population_key',
'        )',
'    )',
    sep = "\n"
)

hits <- gregexpr(
    old,
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
            "Expected exactly one population validation block;",
            "found",
            n_hits
        )
    )
}

text <- sub(
    old,
    new,
    text,
    fixed = TRUE
)

writeLines(
    text,
    con = target,
    useBytes = TRUE
)

cat(
    "Patched population identity validation.\n",
    "Canonical identity: population_key\n",
    "Visible B/D/R label: population_label\n",
    sep = ""
)
