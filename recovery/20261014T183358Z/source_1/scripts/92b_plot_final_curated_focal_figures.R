#!/usr/bin/env Rscript

project <- "/scratch/martinlab/joe/GynandroCalling"

ck <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
)

global_ck <- file.path(
    project,
    "analysis",
    "checkpoints"
)

wide_file <- file.path(
    ck,
    "final_curated_focal_biological_evidence.wide.tsv"
)

long_file <- file.path(
    ck,
    "final_curated_focal_biological_evidence.long.tsv"
)

optix_stable_file <- file.path(
    global_ck,
    "optix_51sample_stable_candidate_sites.tsv"
)

optix_selected_file <- file.path(
    global_ck,
    "optix_51sample_candidate_pattern_groups.tsv"
)

yb_scores_file <- file.path(
    ck,
    "ivory_original_502_comparative_scores.tsv"
)

yb_selected_file <- file.path(
    ck,
    "ivory_comparative_informative_sites_BA80.tsv"
)

out_main_dir <- file.path(
    project,
    "analysis",
    "figures",
    "main_final_curated_focal_evidence"
)

out_supp_dir <- file.path(
    project,
    "analysis",
    "figures",
    "supplement_final_curated_focal_evidence"
)

dir.create(
    out_main_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    out_supp_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

wide <- read.delim(
    wide_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

long <- read.delim(
    long_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

module_order <- c(
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb"
)

display_names <- c(
    Band = "Band",
    Ray1 = "Ray 1",
    Ray2 = "Ray 2",
    Dennis1 = "Dennis 1",
    Dennis2 = "Dennis 2",
    Yb = "Yb"
)

expected <- list(
    Band = c(
        Left = 0,
        Right = 0.5
    ),
    Ray1 = c(
        Left = 1,
        Right = 0.5
    ),
    Ray2 = c(
        Left = 1,
        Right = 0.5
    ),
    Dennis1 = c(
        Left = 1,
        Right = 0.5
    ),
    Dennis2 = c(
        Left = 1,
        Right = 0.5
    ),
    Yb = c(
        Left = 0.5,
        Right = 1
    )
)

expected_gt <- list(
    Band = c(
        Left = "AA",
        Right = "PA"
    ),
    Ray1 = c(
        Left = "PP",
        Right = "PA"
    ),
    Ray2 = c(
        Left = "PP",
        Right = "PA"
    ),
    Dennis1 = c(
        Left = "PP",
        Right = "PA"
    ),
    Dennis2 = c(
        Left = "PP",
        Right = "PA"
    ),
    Yb = c(
        Left = "PA",
        Right = "PP"
    )
)

left_col <- "#0072B2"
right_col <- "#D55E00"

# Canonical analysed intervals.
module_ranges <- list(
    Band = c(
        1112651,
        1114314
    ),
    Ray1 = c(
        1156746,
        1160220
    ),
    Ray2 = c(
        1174547,
        1187833
    ),
    Dennis1 = c(
        1215210,
        1216124
    ),
    Dennis2 = c(
        1218557,
        1221494
    ),
    Yb = c(
        1627000,
        1637000
    )
)

scale_value <- function(
    x,
    from,
    to
) {
    if (diff(from) == 0) {
        return(
            mean(to)
        )
    }

    to[1] +
        (
            (
                x - from[1]
            ) /
            diff(from)
        ) *
        diff(to)
}


draw_zoom_box <- function(
    module,
    xmin,
    xmax,
    ymin,
    ymax,
    show_y_labels = FALSE
) {

    d <- wide[
        wide$module == module,
        ,
        drop = FALSE
    ]

    box_range <- module_ranges[
        [module]
    ]

    rect(
        xmin,
        ymin,
        xmax,
        ymax,
        border = "grey35",
        lwd = 0.9
    )

    title_y <- ymax + 0.014

    text(
        mean(
            c(
                xmin,
                xmax
            )
        ),
        title_y,
        labels = sprintf(
            "%s  (n=%d)   expected L %s | R %s",
            display_names[
                module
            ],
            nrow(d),
            expected_gt[
                [module]
            ][
                "Left"
            ],
            expected_gt[
                [module]
            ][
                "Right"
            ]
        ),
        cex = 0.64,
        font = 2
    )

    dose_to_y <- function(
        dose
    ) {
        scale_value(
            dose,
            c(
                0,
                1
            ),
            c(
                ymin + 0.10 *
                    (
                        ymax - ymin
                    ),
                ymax - 0.10 *
                    (
                        ymax - ymin
                    )
            )
        )
    }

    y0 <- dose_to_y(
        0
    )
    y05 <- dose_to_y(
        0.5
    )
    y1 <- dose_to_y(
        1
    )

    abline_y <- function(
        y,
        col
    ) {
        segments(
            xmin,
            y,
            xmax,
            y,
            col = col,
            lty = 3,
            lwd = 0.8
        )
    }

    abline_y(
        dose_to_y(
            expected[
                [module]
            ][
                "Left"
            ]
        ),
        left_col
    )

    abline_y(
        dose_to_y(
            expected[
                [module]
            ][
                "Right"
            ]
        ),
        right_col
    )

    if (show_y_labels) {

        axis_x <- xmin - 0.007

        text(
            axis_x,
            y0,
            "AA",
            adj = 1,
            cex = 0.60
        )

        text(
            axis_x,
            y05,
            "PA",
            adj = 1,
            cex = 0.60
        )

        text(
            axis_x,
            y1,
            "PP",
            adj = 1,
            cex = 0.60
        )
    }

    if (nrow(d) > 0) {

        xs <- scale_value(
            d$pos,
            box_range,
            c(
                xmin + 0.02 *
                    (
                        xmax - xmin
                    ),
                xmax - 0.02 *
                    (
                        xmax - xmin
                    )
            )
        )

        for (i in seq_len(
            nrow(d)
        )) {

            left_ok <- isTRUE(
                d$left_QC_pass_DP2_GQ20[
                    i
                ]
            )

            right_ok <- isTRUE(
                d$right_QC_pass_DP2_GQ20[
                    i
                ]
            )

            left_dose <- suppressWarnings(
                as.numeric(
                    d$left_presence_dosage[
                        i
                    ]
                )
            )

            right_dose <- suppressWarnings(
                as.numeric(
                    d$right_presence_dosage[
                        i
                    ]
                )
            )

            x <- xs[
                i
            ]

            if (
                left_ok &&
                right_ok &&
                is.finite(
                    left_dose
                ) &&
                is.finite(
                    right_dose
                )
            ) {

                segments(
                    x,
                    dose_to_y(
                        left_dose
                    ),
                    x,
                    dose_to_y(
                        right_dose
                    ),
                    col = "grey55",
                    lwd = 0.75
                )
            }

            if (
                left_ok &&
                is.finite(
                    left_dose
                )
            ) {

                points(
                    x,
                    dose_to_y(
                        left_dose
                    ),
                    pch = 16,
                    cex = 0.72,
                    col = left_col
                )
            }

            if (
                right_ok &&
                is.finite(
                    right_dose
                )
            ) {

                points(
                    x,
                    dose_to_y(
                        right_dose
                    ),
                    pch = 1,
                    cex = 0.98,
                    lwd = 1.1,
                    col = right_col
                )
            }

            segments(
                x,
                ymin,
                x,
                ymin + 0.018 *
                    (
                        ymax - ymin
                    ),
                col = "grey45",
                lwd = 0.6
            )
        }
    }

    text(
        xmin,
        ymin - 0.012,
        sprintf(
            "%.4f",
            box_range[
                1
            ] / 1e6
        ),
        adj = 0,
        cex = 0.52
    )

    text(
        xmax,
        ymin - 0.012,
        sprintf(
            "%.4f Mb",
            box_range[
                2
            ] / 1e6
        ),
        adj = 1,
        cex = 0.52
    )
}


# ============================================================
# MAIN FIGURE
# ============================================================

main_pdf <- file.path(
    out_main_dir,
    "focal_locus_cartoon_curated_comparative_sites.pdf"
)

pdf(
    main_pdf,
    width = 14,
    height = 8.5,
    useDingbats = FALSE
)

par(
    mar = c(
        0,
        0,
        0,
        0
    ),
    xpd = NA
)

plot.new()

plot.window(
    xlim = c(
        0,
        1
    ),
    ylim = c(
        0,
        1
    )
)

text(
    0.05,
    0.965,
    "Focal left-right genotypes at comparative-selected patterning sites",
    adj = 0,
    cex = 1.25,
    font = 2
)

text(
    0.05,
    0.928,
    "Optix sites: curated 51-individual comparative panel; Yb sites: comparative BA >= 0.80. Focal genotypes were not used for site ascertainment.",
    adj = 0,
    cex = 0.73
)

points(
    0.69,
    0.965,
    pch = 16,
    cex = 0.8,
    col = left_col
)

text(
    0.705,
    0.965,
    "Left",
    adj = 0,
    cex = 0.67
)

points(
    0.77,
    0.965,
    pch = 1,
    cex = 1.0,
    lwd = 1.1,
    col = right_col
)

text(
    0.785,
    0.965,
    "Right",
    adj = 0,
    cex = 0.67
)

text(
    0.86,
    0.965,
    "AA=0   PA=0.5   PP=1",
    adj = 0,
    cex = 0.64
)

# ------------------------------------------------------------
# Hmr18 locus overview
# ------------------------------------------------------------

optix_modules <- c(
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2"
)

optix_min <- min(
    vapply(
        module_ranges[
            optix_modules
        ],
        `[`,
        numeric(
            1
        ),
        1
    )
)

optix_max <- max(
    vapply(
        module_ranges[
            optix_modules
        ],
        `[`,
        numeric(
            1
        ),
        2
    )
)

overview_x <- c(
    0.07,
    0.94
)

overview_y <- 0.845

segments(
    overview_x[
        1
    ],
    overview_y,
    overview_x[
        2
    ],
    overview_y,
    lwd = 2.2,
    col = "grey30"
)

text(
    0.04,
    overview_y,
    "Hmr18",
    adj = 1,
    cex = 0.76,
    font = 2
)

overview_centres <- numeric(
    length(
        optix_modules
    )
)

names(
    overview_centres
) <- optix_modules

for (module in optix_modules) {

    rr <- module_ranges[
        [module]
    ]

    xx <- scale_value(
        rr,
        c(
            optix_min,
            optix_max
        ),
        overview_x
    )

    rect(
        xx[
            1
        ],
        overview_y - 0.018,
        xx[
            2
        ],
        overview_y + 0.018,
        border = "grey25",
        col = "grey90"
    )

    overview_centres[
        module
    ] <- mean(
        xx
    )

    text(
        mean(
            xx
        ),
        overview_y + 0.036,
        display_names[
            module
        ],
        cex = 0.57
    )

    d <- wide[
        wide$module == module,
        ,
        drop = FALSE
    ]

    if (nrow(d) > 0) {

        site_x <- scale_value(
            d$pos,
            c(
                optix_min,
                optix_max
            ),
            overview_x
        )

        segments(
            site_x,
            overview_y - 0.027,
            site_x,
            overview_y + 0.027,
            lwd = 0.55,
            col = "grey20"
        )
    }
}

text(
    overview_x[
        1
    ],
    overview_y - 0.044,
    sprintf(
        "%.3f Mb",
        optix_min / 1e6
    ),
    adj = 0,
    cex = 0.54
)

text(
    overview_x[
        2
    ],
    overview_y - 0.044,
    sprintf(
        "%.3f Mb",
        optix_max / 1e6
    ),
    adj = 1,
    cex = 0.54
)

# ------------------------------------------------------------
# Optix zoom boxes.
#
# Widths scale with genomic span, with a 2-kb minimum
# so very small modules remain visible.
# ------------------------------------------------------------

spans <- vapply(
    module_ranges[
        optix_modules
    ],
    function(x) {
        diff(x)
    },
    numeric(
        1
    )
)

display_spans <- pmax(
    spans,
    2000
)

total_x <- 0.88

gap <- 0.012

usable_x <- (
    total_x
    -
    gap *
    (
        length(
            optix_modules
        )
        - 1
    )
)

widths <- (
    usable_x
    *
    display_spans
    /
    sum(
        display_spans
    )
)

starts <- numeric(
    length(
        widths
    )
)

starts[
    1
] <- 0.07

if (
    length(
        starts
    ) > 1
) {

    for (i in 2:length(
        starts
    )) {

        starts[
            i
        ] <- (
            starts[
                i - 1
            ]
            +
            widths[
                i - 1
            ]
            +
            gap
        )
    }
}

ends <- starts + widths

zoom_y <- c(
    0.47,
    0.72
)

for (i in seq_along(
    optix_modules
)) {

    module <- optix_modules[
        i
    ]

    segments(
        overview_centres[
            module
        ],
        overview_y - 0.035,
        overview_centres[
            module
        ],
        0.765,
        col = "grey75",
        lwd = 0.6
    )

    segments(
        overview_centres[
            module
        ],
        0.765,
        mean(
            c(
                starts[
                    i
                ],
                ends[
                    i
                ]
            )
        ),
        0.765,
        col = "grey75",
        lwd = 0.6
    )

    segments(
        mean(
            c(
                starts[
                    i
                ],
                ends[
                    i
                ]
            )
        ),
        0.765,
        mean(
            c(
                starts[
                    i
                ],
                ends[
                    i
                ]
            )
        ),
        zoom_y[
            2
        ] + 0.025,
        col = "grey75",
        lwd = 0.6
    )

    draw_zoom_box(
        module,
        starts[
            i
        ],
        ends[
            i
        ],
        zoom_y[
            1
        ],
        zoom_y[
            2
        ],
        show_y_labels = (
            i == 1
        )
    )
}

# ------------------------------------------------------------
# Yb separate chromosome/locus
# ------------------------------------------------------------

yb_overview_y <- 0.35

segments(
    0.17,
    yb_overview_y,
    0.83,
    yb_overview_y,
    lwd = 2.2,
    col = "grey30"
)

text(
    0.14,
    yb_overview_y,
    "Hmr15",
    adj = 1,
    cex = 0.76,
    font = 2
)

rr <- module_ranges[
    ["Yb"]
]

yb_d <- wide[
    wide$module == "Yb",
    ,
    drop = FALSE
]

if (nrow(
    yb_d
) > 0) {

    site_x <- scale_value(
        yb_d$pos,
        rr,
        c(
            0.17,
            0.83
        )
    )

    segments(
        site_x,
        yb_overview_y - 0.024,
        site_x,
        yb_overview_y + 0.024,
        lwd = 0.7
    )
}

text(
    0.50,
    yb_overview_y + 0.036,
    "Yb / cortex-ivory interval",
    cex = 0.63
)

draw_zoom_box(
    "Yb",
    0.23,
    0.77,
    0.075,
    0.255,
    show_y_labels = TRUE
)

text(
    0.50,
    0.025,
    "Vertical connectors join left and right observations at the same SNP. Failed DP>=2/GQ>=20 calls remain as gaps.",
    cex = 0.62
)

dev.off()


# ============================================================
# SUPPLEMENT 1:
# COMPARATIVE SITE ASCERTAINMENT
# ============================================================

stable <- read.delim(
    optix_stable_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

selected_optix <- read.delim(
    optix_selected_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

map_module <- function(
    x
) {
    x <- tolower(
        gsub(
            "[^a-z0-9]",
            "",
            x
        )
    )

    if (grepl(
        "band",
        x
    )) {
        return(
            "Band"
        )
    }

    if (grepl(
        "ray1",
        x
    )) {
        return(
            "Ray1"
        )
    }

    if (grepl(
        "ray2",
        x
    )) {
        return(
            "Ray2"
        )
    }

    if (grepl(
        "dennis1",
        x
    )) {
        return(
            "Dennis1"
        )
    }

    if (grepl(
        "dennis2",
        x
    )) {
        return(
            "Dennis2"
        )
    }

    NA_character_
}

stable$module <- vapply(
    stable$panel,
    map_module,
    character(
        1
    )
)

selected_optix$module <- vapply(
    selected_optix$panel,
    map_module,
    character(
        1
    )
)

stable$minimum_BA <- as.numeric(
    stable$minimum_BA
)

selected_optix$minimum_BA <- as.numeric(
    selected_optix$minimum_BA
)

yb_scores <- read.delim(
    yb_scores_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

yb_selected <- read.delim(
    yb_selected_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

yb_scores$balanced_accuracy <- as.numeric(
    yb_scores$balanced_accuracy
)

yb_selected$minimum_BA <- as.numeric(
    yb_selected$minimum_BA
)

ascertainment_pdf <- file.path(
    out_supp_dir,
    "comparative_site_ascertainment.pdf"
)

pdf(
    ascertainment_pdf,
    width = 12,
    height = 8,
    useDingbats = FALSE
)

par(
    mfrow = c(
        2,
        3
    ),
    mar = c(
        4.0,
        4.2,
        2.8,
        1.0
    ),
    oma = c(
        0,
        0,
        2.2,
        0
    )
)

for (module in module_order) {

    if (module != "Yb") {

        all_d <- stable[
            stable$module == module,
            ,
            drop = FALSE
        ]

        sel_d <- selected_optix[
            selected_optix$module == module,
            ,
            drop = FALSE
        ]

        plot(
            all_d$pos / 1e6,
            all_d$minimum_BA,
            pch = 16,
            cex = 0.42,
            col = "grey72",
            xlab = "Hmr18 position (Mb)",
            ylab = "Minimum comparative BA",
            ylim = c(
                0.45,
                1.02
            ),
            main = display_names[
                module
            ],
            bty = "l"
        )

        abline(
            h = 0.80,
            lty = 2,
            lwd = 0.8
        )

        points(
            sel_d$pos / 1e6,
            sel_d$minimum_BA,
            pch = 21,
            bg = left_col,
            col = left_col,
            cex = 0.75
        )

    } else {

        plot(
            yb_scores$pos / 1e6,
            yb_scores$balanced_accuracy,
            pch = 16,
            cex = 0.42,
            col = "grey72",
            xlab = "Hmr15 position (Mb)",
            ylab = "Comparative BA",
            ylim = c(
                0.45,
                1.02
            ),
            main = "Yb",
            bty = "l"
        )

        abline(
            h = 0.80,
            lty = 2,
            lwd = 0.8
        )

        points(
            yb_selected$pos / 1e6,
            yb_selected$minimum_BA,
            pch = 21,
            bg = left_col,
            col = left_col,
            cex = 0.75
        )
    }
}

mtext(
    "Comparative-only site ascertainment; focal GT, AD, DP, GQ and expected focal zygosity were not used to select sites",
    outer = TRUE,
    side = 3,
    cex = 0.82,
    font = 2
)

dev.off()


# ============================================================
# SUPPLEMENT 2:
# FULL FOCAL PRESENCE-ALLELE FRACTIONS
# ============================================================

fraction_pdf <- file.path(
    out_supp_dir,
    "full_focal_presence_allele_fraction.pdf"
)

pdf(
    fraction_pdf,
    width = 12,
    height = 8,
    useDingbats = FALSE
)

par(
    mfrow = c(
        2,
        3
    ),
    mar = c(
        4.2,
        4.1,
        3.0,
        1.0
    ),
    oma = c(
        0,
        0,
        1.8,
        0
    )
)

for (module in module_order) {

    d <- wide[
        wide$module == module,
        ,
        drop = FALSE
    ]

    x <- seq_len(
        nrow(d)
    )

    left_y <- as.numeric(
        d$left_presence_allele_fraction
    )

    right_y <- as.numeric(
        d$right_presence_allele_fraction
    )

    plot(
        NA,
        xlim = c(
            0.5,
            max(
                1.5,
                nrow(d) + 0.5
            )
        ),
        ylim = c(
            -0.03,
            1.03
        ),
        xlab = "Selected SNP",
        ylab = "Presence-associated read fraction",
        main = sprintf(
            "%s  (n=%d)",
            display_names[
                module
            ],
            nrow(d)
        ),
        bty = "l"
    )

    for (i in seq_along(
        x
    )) {

        if (
            is.finite(
                left_y[
                    i
                ]
            ) &&
            is.finite(
                right_y[
                    i
                ]
            )
        ) {

            segments(
                x[
                    i
                ],
                left_y[
                    i
                ],
                x[
                    i
                ],
                right_y[
                    i
                ],
                col = "grey60",
                lwd = 0.8
            )
        }
    }

    points(
        x,
        left_y,
        pch = 16,
        cex = 0.72,
        col = left_col
    )

    points(
        x,
        right_y,
        pch = 1,
        cex = 0.98,
        lwd = 1.1,
        col = right_col
    )

    abline(
        h = c(
            0,
            0.5,
            1
        ),
        col = "grey88",
        lty = 3
    )
}

mtext(
    "Focal read support at the frozen comparative-selected sites; QC failures remain as gaps",
    outer = TRUE,
    side = 3,
    cex = 0.83,
    font = 2
)

dev.off()


cat(
    "\nOutputs:\n",
    main_pdf,
    "\n",
    ascertainment_pdf,
    "\n",
    fraction_pdf,
    "\n",
    sep = ""
)
