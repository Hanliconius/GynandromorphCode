#!/usr/bin/env Rscript

target <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "158_plot_science_condensation_trials.R"
)

if (!file.exists(target)) {
    stop(
        paste(
            "Missing Analysis-158 target:",
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
'    )',
'    ordered = TRUE',
')',
    sep = "\n"
)

new <- paste(
'    ),',
'    ordered = TRUE',
')',
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

cat(
    "Malformed factor() endings found: ",
    n_hits,
    "\n",
    sep = ""
)

if (n_hits != 2L) {
    stop(
        paste(
            "Expected exactly 2 malformed factor() calls;",
            "found",
            n_hits,
            "- target left unchanged."
        )
    )
}

text <- gsub(
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
    "Fixed both factor() calls in:\n",
    target,
    "\n",
    sep = ""
)
