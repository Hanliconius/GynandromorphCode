#!/usr/bin/env Rscript

target <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "159_plot_science_trials_pubregions_cartoons.R"
)

if (!file.exists(target)) {
    stop(
        paste(
            "Missing Analysis-159 target:",
            target
        )
    )
}

lines <- readLines(
    target,
    warn = FALSE
)


# ============================================================
# REPLACE PLACEHOLDER PUBLISHED-REGION TABLE
# ============================================================

start <- which(
    trimws(lines) ==
        "published_regions <- data.frame("
)

end <- which(
    grepl(
        "^if \\(any\\(!is\\.finite\\(published_regions\\$xmin\\)",
        trimws(lines)
    )
)

if (length(start) != 1L) {
    stop(
        paste(
            "Expected one published_regions start;",
            "found",
            length(start)
        )
    )
}

if (length(end) != 1L) {
    stop(
        paste(
            "Expected one published_regions validation block;",
            "found",
            length(end)
        )
    )
}

if (end <= start) {
    stop(
        "published_regions block boundaries are invalid."
    )
}


new_regions <- c(
'published_regions <- data.frame(',
'    region = c(',
'        "B",',
'        "R1",',
'        "R2",',
'        "D1",',
'        "D2"',
'    ),',
'',
'    module = c(',
'        "B",',
'        "R",',
'        "R",',
'        "D",',
'        "D"',
'    ),',
'',
'    xmin = c(',
'        1112651,',
'        1156746,',
'        1174547,',
'        1215210,',
'        1218557',
'    ) / 1e6,',
'',
'    xmax = c(',
'        1114314,',
'        1160220,',
'        1187833,',
'        1216124,',
'        1221494',
'    ) / 1e6,',
'',
'    colour = c(',
'        trait_colours[["B"]],',
'        trait_colours[["R"]],',
'        trait_colours[["R"]],',
'        trait_colours[["D"]],',
'        trait_colours[["D"]]',
'    ),',
'',
'    stringsAsFactors = FALSE',
')',
''
)

lines <- c(
    lines[
        seq_len(
            start - 1L
        )
    ],
    new_regions,
    lines[
        end:length(lines)
    ]
)


# ============================================================
# REPLACE REGION-BAR VERTICAL POSITIONS
#
# Most intervals share one thin annotation row.
# D1/D2 are very close, so D2 is offset slightly to avoid
# label collisions.
# ============================================================

start2 <- which(
    trimws(lines) ==
        "pub_trial1 <- published_regions"
)

cartoon_marker <- which(
    trimws(lines) ==
        "# BUTTERFLY PATTERN CARTOONS"
)

if (length(start2) != 1L) {
    stop(
        paste(
            "Expected one pub_trial1 position block;",
            "found",
            length(start2)
        )
    )
}

if (length(cartoon_marker) != 1L) {
    stop(
        paste(
            "Expected one butterfly-cartoon marker;",
            "found",
            length(cartoon_marker)
        )
    )
}

# Remove the separator immediately above the cartoon heading.
position_end <- cartoon_marker - 2L

new_positions <- c(
'pub_trial1 <- published_regions',
'',
'pub_trial1$y <- c(',
'    4.18,   # B',
'    4.18,   # R1',
'    4.18,   # R2',
'    4.18,   # D1',
'    4.06    # D2; offset because D1/D2 are adjacent',
')',
'',
'pub_trial1$label_y <- pub_trial1$y + 0.045',
'',
'',
'pub_trial2 <- published_regions',
'',
'pub_trial2$y <- c(',
'    3.27,   # B',
'    3.27,   # R1',
'    3.27,   # R2',
'    3.27,   # D1',
'    3.16    # D2; offset because D1/D2 are adjacent',
')',
'',
'pub_trial2$label_y <- pub_trial2$y + 0.040',
''
)

lines <- c(
    lines[
        seq_len(
            start2 - 1L
        )
    ],
    new_positions,
    lines[
        (position_end + 1L):length(lines)
    ]
)


writeLines(
    lines,
    con = target,
    useBytes = TRUE
)

cat(
    "Restored published intervals in Analysis 159:\n\n",
    "  B   1.112651 - 1.114314 Mb\n",
    "  R1  1.156746 - 1.160220 Mb\n",
    "  R2  1.174547 - 1.187833 Mb\n",
    "  D1  1.215210 - 1.216124 Mb\n",
    "  D2  1.218557 - 1.221494 Mb\n\n",
    "Patched:\n",
    target,
    "\n",
    sep = ""
)
