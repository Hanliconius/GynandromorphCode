#!/usr/bin/env Rscript

builder <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "156_make_BDR_ordered_canonical_taxon_plotter.R"
)

if (!file.exists(builder)) {
    stop(
        paste(
            "Missing 156 builder:",
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
'resolve_column <- function(',
'    x,',
'    candidates',
') {',
'',
'    cleaned <- clean_name(',
'        names(x)',
'    )',
'',
'    for (',
'        candidate',
'        in candidates',
'    ) {',
'',
'        hit <- which(',
'            tolower(cleaned) ==',
'                tolower(candidate)',
'        )',
'',
'        if (length(hit) == 1L) {',
'            return(',
'                names(x)[hit]',
'            )',
'        }',
'    }',
'',
'    stop(',
'        paste(',
'            "Could not resolve column from:",',
'            paste(',
'                candidates,',
'                collapse = ", "',
'            )',
'        )',
'    )',
'}',
    sep = "\n"
)


new_block <- paste(
'# Metadata headers in this project can contain punctuation,',
'# whitespace and embedded line breaks. Match columns using the',
'# established normalize-to-alphanumeric strategy.',
'normalize_column_name <- function(x) {',
'',
'    x <- tolower(',
'        as.character(x)',
'    )',
'',
'    gsub(',
'        "[^a-z0-9]+",',
'        "",',
'        x,',
'        perl = TRUE',
'    )',
'}',
'',
'',
'resolve_column <- function(',
'    x,',
'    candidates',
') {',
'',
'    normalized_columns <- normalize_column_name(',
'        names(x)',
'    )',
'',
'    for (',
'        candidate',
'        in candidates',
'    ) {',
'',
'        target <- normalize_column_name(',
'            candidate',
'        )',
'',
'        hit <- which(',
'            normalized_columns ==',
'                target',
'        )',
'',
'        if (length(hit) == 1L) {',
'',
'            cat(',
'                sprintf(',
'                    "Resolved metadata column %-10s -> %s\\n",',
'                    candidate,',
'                    names(x)[hit]',
'                )',
'            )',
'',
'            return(',
'                names(x)[hit]',
'            )',
'        }',
'    }',
'',
'    stop(',
'        paste0(',
'            "Could not resolve metadata column from: ",',
'            paste(',
'                candidates,',
'                collapse = ", "',
'            ),',
'            "\\nAvailable normalized columns:\\n",',
'            paste(',
'                paste0(',
'                    names(x),',
'                    " -> ",',
'                    normalized_columns',
'                ),',
'                collapse = "\\n"',
'            )',
'        )',
'    )',
'}',
    sep = "\n"
)


if (!grepl(
    old_block,
    text,
    fixed = TRUE
)) {
    stop(
        paste(
            "Could not find the original 156 resolve_column()",
            "function exactly; builder left unchanged."
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
    "Patched metadata-column resolver in:\n",
    builder,
    "\n",
    sep = ""
)
