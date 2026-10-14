#!/usr/bin/env Rscript

project <- "/scratch/martinlab/joe/GynandroCalling"

ck <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
)

wide_file <- file.path(
    ck,
    "focal_calls_at_frozen_sites_biological.wide.tsv"
)

sites_file <- file.path(
    ck,
    "frozen_comparative_sites_reference_polarized.tsv"
)

out_main <- file.path(
    project,
    "analysis",
    "figures",
    "main_final_merged_lr_focal_evidence"
)

out_supp <- file.path(
    project,
    "analysis",
    "figures",
    "supplement_final_merged_lr_focal_evidence"
)

dir.create(
    out_main,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    out_supp,
    recursive = TRUE,
    showWarnings = FALSE
)

wide <- read.delim(
    wide_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sites <- read.delim(
    sites_file,
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

expected_dose <- list(
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

# Use exact analysed interval boundaries for the locus cartoons.
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
            rep(
                mean(to),
                length(x)
            )
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


get_bool <- function(
    x
) {

    if (is.logical(x)) {
        return(x)
    }

    tolower(
        as.character(x)
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
}


# ============================================================
# Work out actual column names in the 90c output.
# ============================================================

find_col <- function(
    candidates,
    required = TRUE
) {

    hit <- candidates[
        candidates %in% names(wide)
    ]

    if (length(hit) > 0) {
        return(
            hit[1]
        )
    }

    if (required) {
        stop(
            paste(
                "Could not find any of:",
                paste(
                    candidates,
                    collapse = ", "
                )
            )
        )
    }

    NA_character_
}


left_dose_col <- find_col(
    c(
        "left_presence_dosage",
        "TN_Left_presence_dosage"
    )
)

right_dose_col <- find_col(
    c(
        "right_presence_dosage",
        "TN_Right_presence_dosage"
    )
)

left_frac_col <- find_col(
    c(
        "left_presence_allele_fraction",
        "TN_Left_presence_allele_fraction"
    )
)

right_frac_col <- find_col(
    c(
        "right_presence_allele_fraction",
        "TN_Right_presence_allele_fraction"
    )
)

left_qc_col <- find_col(
    c(
        "left_QC_pass_DP2_GQ20",
        "left_passes_DP2_GQ20",
        "TN_Left_QC_pass_DP2_GQ20"
    )
)

right_qc_col <- find_col(
    c(
        "right_QC_pass_DP2_GQ20",
        "right_passes_DP2_GQ20",
        "TN_Right_QC_pass_DP2_GQ20"
    )
)


wide[[left_dose_col]] <- suppressWarnings(
    as.numeric(
        wide[[left_dose_col]]
    )
)

wide[[right_dose_col]] <- suppressWarnings(
    as.numeric(
        wide[[right_dose_col]]
    )
)

wide[[left_frac_col]] <- suppressWarnings(
    as.numeric(
        wide[[left_frac_col]]
    )
)

wide[[right_frac_col]] <- suppressWarnings(
    as.numeric(
        wide[[right_frac_col]]
    )
)

wide[[left_qc_col]] <- get_bool(
    wide[[left_qc_col]]
)

wide[[right_qc_col]] <- get_bool(
    wide[[right_qc_col]]
)


# ============================================================
# MAIN FIGURE HELPERS
# ============================================================

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

    rr <- module_ranges[[module]]

    rect(
        xmin,
        ymin,
        xmax,
        ymax,
        border = "grey35",
        lwd = 0.9
    )

    text(
        mean(
            c(
                xmin,
                xmax
            )
        ),
        ymax + 0.018,
        sprintf(
            "%s  (n=%d)   expected L %s | R %s",
            display_names[[module]],
            nrow(d),
            expected_gt[[module]][["Left"]],
            expected_gt[[module]][["Right"]]
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

    for (
        dose in c(
            0,
            0.5,
            1
        )
    ) {

        segments(
            xmin,
            dose_to_y(
                dose
            ),
            xmax,
            dose_to_y(
                dose
            ),
            col = "grey90",
            lwd = 0.6
        )
    }

    segments(
        xmin,
        dose_to_y(
            expected_dose[[module]][["Left"]]
        ),
        xmax,
        dose_to_y(
            expected_dose[[module]][["Left"]]
        ),
        col = left_col,
        lty = 3,
        lwd = 0.9
    )

    segments(
        xmin,
        dose_to_y(
            expected_dose[[module]][["Right"]]
        ),
        xmax,
        dose_to_y(
            expected_dose[[module]][["Right"]]
        ),
        col = right_col,
        lty = 3,
        lwd = 0.9
    )

    if (show_y_labels) {

        text(
            xmin - 0.007,
            dose_to_y(
                0
            ),
            "AA",
            adj = 1,
            cex = 0.59
        )

        text(
            xmin - 0.007,
            dose_to_y(
                0.5
            ),
            "PA",
            adj = 1,
            cex = 0.59
        )

        text(
            xmin - 0.007,
            dose_to_y(
                1
            ),
            "PP",
            adj = 1,
            cex = 0.59
        )
    }

    if (nrow(d) > 0) {

        xs <- scale_value(
            d$pos,
            rr,
            c(
                xmin + 0.025 *
                    (
                        xmax - xmin
                    ),
                xmax - 0.025 *
                    (
                        xmax - xmin
                    )
            )
        )

        for (
            i in seq_len(
                nrow(d)
            )
        ) {

            left_ok <- (
                isTRUE(
                    d[[left_qc_col]][i]
                ) &&
                is.finite(
                    d[[left_dose_col]][i]
                )
            )

            right_ok <- (
                isTRUE(
                    d[[right_qc_col]][i]
                ) &&
                is.finite(
                    d[[right_dose_col]][i]
                )
            )

            x <- xs[i]

            # Vertical line = same genomic SNP.
            if (
                left_ok &&
                right_ok
            ) {

                segments(
                    x,
                    dose_to_y(
                        d[[left_dose_col]][i]
                    ),
                    x,
                    dose_to_y(
                        d[[right_dose_col]][i]
                    ),
                    col = "grey55",
                    lwd = 0.75
                )
            }

            if (left_ok) {

                points(
                    x,
                    dose_to_y(
                        d[[left_dose_col]][i]
                    ),
                    pch = 16,
                    cex = 0.74,
                    col = left_col
                )
            }

            if (right_ok) {

                points(
                    x,
                    dose_to_y(
                        d[[right_dose_col]][i]
                    ),
                    pch = 1,
                    cex = 1.0,
                    lwd = 1.15,
                    col = right_col
                )
            }

            # Every frozen site gets a genomic tick,
            # including sites with unavailable focal calls.
            segments(
                x,
                ymin,
                x,
                ymin + 0.02 *
                    (
                        ymax - ymin
                    ),
                col = "grey30",
                lwd = 0.65
            )
        }
    }

    text(
        xmin,
        ymin - 0.014,
        sprintf(
            "%.4f",
            rr[1] / 1e6
        ),
        adj = 0,
        cex = 0.52
    )

    text(
        xmax,
        ymin - 0.014,
        sprintf(
            "%.4f Mb",
            rr[2] / 1e6
        ),
        adj = 1,
        cex = 0.52
    )
}


# ============================================================
# MAIN FIGURE
# ============================================================

main_pdf <- file.path(
    out_main,
    "focal_locus_cartoon_merged_left_right.pdf"
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
    "Merged left-right focal genotypes at comparative-selected sites",
    adj = 0,
    cex = 1.22,
    font = 2
)

text(
    0.05,
    0.929,
    "Optix sites were selected from the curated comparative panel before TN_Left and TN_Right were evaluated; unavailable focal measurements remain as gaps.",
    adj = 0,
    cex = 0.70
)

points(
    0.72,
    0.966,
    pch = 16,
    cex = 0.78,
    col = left_col
)

text(
    0.735,
    0.966,
    "TN_Left",
    adj = 0,
    cex = 0.64
)

points(
    0.82,
    0.966,
    pch = 1,
    cex = 1.0,
    lwd = 1.1,
    col = right_col
)

text(
    0.835,
    0.966,
    "TN_Right",
    adj = 0,
    cex = 0.64
)


# ------------------------------------------------------------
# Hmr18 overview
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
        function(x) {
            x[1]
        },
        numeric(
            1
        )
    )
)

optix_max <- max(
    vapply(
        module_ranges[
            optix_modules
        ],
        function(x) {
            x[2]
        },
        numeric(
            1
        )
    )
)

overview_range <- c(
    0.07,
    0.94
)

overview_y <- 0.845

segments(
    overview_range[1],
    overview_y,
    overview_range[2],
    overview_y,
    lwd = 2,
    col = "grey30"
)

text(
    0.04,
    overview_y,
    "Hmr18",
    adj = 1,
    cex = 0.75,
    font = 2
)

overview_centres <- setNames(
    numeric(
        length(
            optix_modules
        )
    ),
    optix_modules
)

for (
    module in optix_modules
) {

    rr <- module_ranges[[module]]

    xx <- scale_value(
        rr,
        c(
            optix_min,
            optix_max
        ),
        overview_range
    )

    rect(
        xx[1],
        overview_y - 0.017,
        xx[2],
        overview_y + 0.017,
        border = "grey30",
        col = "grey91"
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
        display_names[[module]],
        cex = 0.55
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
            overview_range
        )

        segments(
            site_x,
            overview_y - 0.026,
            site_x,
            overview_y + 0.026,
            col = "grey25",
            lwd = 0.5
        )
    }
}


# ------------------------------------------------------------
# Expanded module panels:
# widths proportional to genomic interval, minimum 2 kb.
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

gap <- 0.012

usable_width <- (
    0.88 -
    gap *
    (
        length(
            optix_modules
        ) - 1
    )
)

widths <- (
    usable_width *
    display_spans /
    sum(
        display_spans
    )
)

starts <- numeric(
    length(
        widths
    )
)

starts[1] <- 0.07

if (
    length(
        starts
    ) > 1
) {

    for (
        i in 2:length(
            starts
        )
    ) {

        starts[i] <- (
            starts[i - 1] +
            widths[i - 1] +
            gap
        )
    }
}

ends <- starts + widths

for (
    i in seq_along(
        optix_modules
    )
) {

    module <- optix_modules[i]

    centre <- mean(
        c(
            starts[i],
            ends[i]
        )
    )

    # Connector from physical locus to zoom;
    # not data, so elbowed rather than implying scale.
    segments(
        overview_centres[module],
        overview_y - 0.035,
        overview_centres[module],
        0.765,
        col = "grey78",
        lwd = 0.55
    )

    segments(
        overview_centres[module],
        0.765,
        centre,
        0.765,
        col = "grey78",
        lwd = 0.55
    )

    segments(
        centre,
        0.765,
        centre,
        0.745,
        col = "grey78",
        lwd = 0.55
    )

    draw_zoom_box(
        module,
        starts[i],
        ends[i],
        0.47,
        0.72,
        show_y_labels = (
            i == 1
        )
    )
}


# ------------------------------------------------------------
# Separate Yb locus
# ------------------------------------------------------------

text(
    0.04,
    0.346,
    "Hmr15",
    adj = 1,
    cex = 0.75,
    font = 2
)

segments(
    0.17,
    0.346,
    0.83,
    0.346,
    lwd = 2,
    col = "grey30"
)

yb_d <- wide[
    wide$module == "Yb",
    ,
    drop = FALSE
]

if (nrow(yb_d) > 0) {

    xx <- scale_value(
        yb_d$pos,
        module_ranges[["Yb"]],
        c(
            0.17,
            0.83
        )
    )

    segments(
        xx,
        0.321,
        xx,
        0.371,
        col = "grey25",
        lwd = 0.6
    )
}

text(
    0.50,
    0.382,
    "Yb / cortex-ivory interval",
    cex = 0.60
)

draw_zoom_box(
    "Yb",
    0.23,
    0.77,
    0.075,
    0.245,
    show_y_labels = TRUE
)

text(
    0.50,
    0.025,
    "AA, PA and PP denote 0, 1 and 2 copies of the phenotype-presence-associated allele. Vertical connectors link the same SNP across merged left and right samples.",
    cex = 0.59
)

dev.off()


# ============================================================
# SUPPLEMENT:
# ACTUAL P-ALLELE READ FRACTIONS AT EVERY FROZEN SITE
# ============================================================

supp_pdf <- file.path(
    out_supp,
    "merged_left_right_presence_allele_fraction.pdf"
)

pdf(
    supp_pdf,
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
        4.2,
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

for (
    module in module_order
) {

    d <- wide[
        wide$module == module,
        ,
        drop = FALSE
    ]

    x <- seq_len(
        nrow(d)
    )

    ly <- d[[left_frac_col]]
    ry <- d[[right_frac_col]]

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
        xlab = "Comparative-selected SNP",
        ylab = "P-associated read fraction",
        main = sprintf(
            "%s  (n=%d)",
            display_names[[module]],
            nrow(d)
        ),
        bty = "l"
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

    for (
        i in seq_along(
            x
        )
    ) {

        left_ok <- (
            isTRUE(
                d[[left_qc_col]][i]
            ) &&
            is.finite(
                ly[i]
            )
        )

        right_ok <- (
            isTRUE(
                d[[right_qc_col]][i]
            ) &&
            is.finite(
                ry[i]
            )
        )

        if (
            left_ok &&
            right_ok
        ) {

            segments(
                x[i],
                ly[i],
                x[i],
                ry[i],
                col = "grey60",
                lwd = 0.75
            )
        }

        if (left_ok) {

            points(
                x[i],
                ly[i],
                pch = 16,
                col = left_col,
                cex = 0.72
            )
        }

        if (right_ok) {

            points(
                x[i],
                ry[i],
                pch = 1,
                col = right_col,
                cex = 0.95,
                lwd = 1.1
            )
        }
    }
}

mtext(
    "Merged TN_Left and TN_Right read support at the frozen comparative-selected sites; unavailable or QC-failed focal calls are gaps",
    outer = TRUE,
    side = 3,
    font = 2,
    cex = 0.78
)

dev.off()


cat(
    "\nMain figure:\n  ",
    main_pdf,
    "\n\nSupplement:\n  ",
    supp_pdf,
    "\n",
    sep = ""
)
