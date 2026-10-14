#!/usr/bin/env Rscript

builder <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "156_make_BDR_ordered_canonical_taxon_plotter.R"
)

if (!file.exists(builder)) {
    stop(
        paste(
            "Missing builder:",
            builder
        )
    )
}

text <- paste(
    readLines(
        builder,
        warn = FALSE
    ),
    collapse = "\n"
)

old_block <- paste(
'meta <- read.delim(',
'    metadata_file,',
'    sep = "\\t",',
'    header = TRUE,',
'    stringsAsFactors = FALSE,',
'    check.names = FALSE,',
'    quote = "",',
'    comment.char = ""',
')',
    sep = "\n"
)

new_block <- paste(
'meta <- read.delim(',
'    metadata_file,',
'    sep = "\\t",',
'    header = TRUE,',
'    stringsAsFactors = FALSE,',
'    check.names = FALSE,',
'    quote = "\\"",',
'    comment.char = ""',
')',
    sep = "\n"
)

hits <- gregexpr(
    old_block,
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
            "Expected exactly one metadata read.delim block;",
            "found",
            n_hits
        )
    )
}

text <- sub(
    old_block,
    new_block,
    text,
    fixed = TRUE
)

writeLines(
    text,
    con = builder,
    useBytes = TRUE
)

cat(
    "Fixed metadata TSV quote handling in:\n",
    builder,
    "\n",
    sep = ""
)
