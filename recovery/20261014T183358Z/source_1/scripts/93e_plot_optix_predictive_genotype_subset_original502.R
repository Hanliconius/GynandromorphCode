#!/usr/bin/env Rscript

project <- "/scratch/martinlab/joe/GynandroCalling"

ck <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
)

long_file <- file.path(
    ck,
    "optix_original502_selected_sites_genotype_plot.long.tsv"
)

sites_file <- file.path(
    ck,
    "optix_original502_selected_sites_for_plotting.tsv"
)

out_dir <- file.path(
    project,
    "analysis",
    "figures",
    "supplement_optix_predictive_genotype_plots"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

long <- read.delim(
    long_file,
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
    "Ray2",
    "Dennis2"
)

display_names <- c(
    Band = "Band",
    Ray2 = "Ray 2",
    Dennis2 = "Dennis 2"
)

# Requested Whiting-style palette
gt_cols <- c(
    AA = "#2F6BFF",
    PA = "#F39C12",
    PP = "#D81BCE",
    missing = "#E5E7EB"
)

pheno_cols <- c(
    absent = "#FFFFFF",
    present = "#111111"
)

clean_label <- function(x) {

    x <- as.character(x)

    x <- gsub(
        "[[:space:]]+",
        " ",
        x
    )

    trimws(x)
}

thin_axis_labels <- function(labels) {

    n <- length(labels)

    if (n <= 15) {
        step <- 1
    } else if (n <= 30) {
        step <- 2
    } else if (n <= 70) {
        step <- 4
    } else {
        step <- 5
    }

    out <- rep(
        "",
        n
    )

    keep <- seq(
        1,
        n,
        by = step
    )

    out[keep] <- labels[keep]

    out
}

pdf_file <- file.path(
    out_dir,
    paste0(
        "optix_predictive_score_plus_genotype_matrix.",
        "original502.Band_Ray2_Dennis2.pdf"
    )
)

# Large vector pages are intentional:
# 481 individual rows need genuine vertical space.
pdf(
    pdf_file,
    width = 16.5,
    height = 55,
    onefile = TRUE,
    useDingbats = FALSE
)

# Each module uses two panels:
# short BA track + much taller genotype matrix.
layout(
    matrix(
        c(
            1,
            2
        ),
        ncol = 1
    ),
    heights = c(
        1.2,
        17
    )
)

for (module in module_order) {

    d <- long[
        long$module == module,
        ,
        drop = FALSE
    ]

    s <- sites[
        sites$module == module,
        ,
        drop = FALSE
    ]

    if (
        nrow(d) == 0 ||
        nrow(s) == 0
    ) {
        next
    }

    s <- s[
        order(
            s$site_order_within_module,
            s$pos
        ),
        ,
        drop = FALSE
    ]

    sample_df <- unique(
        d[
            ,
            c(
                "sample",
                "sample_order",
                "taxon_name",
                "phenotype_label",
                "phenotype_value"
            )
        ]
    )

    sample_df <- sample_df[
        order(
            sample_df$sample_order
        ),
        ,
        drop = FALSE
    ]

    # Display taxon only.
    # ERR/SRR accession remains an internal join key,
    # but is never printed.
    sample_df$plot_label <- clean_label(
        sample_df$taxon_name
    )

    x_pos <- seq_len(
        nrow(s)
    )


    # ========================================================
    # PANEL A — FULL-PANEL BALANCED ACCURACY
    # ========================================================

    par(
        mar = c(
            0.8,
            18.5,
            3.0,
            4.5
        )
    )

    plot(
        NA,
        xlim = c(
            0.5,
            nrow(s) + 0.5
        ),
        ylim = c(
            0.45,
            1.02
        ),
        xaxt = "n",
        xlab = "",
        ylab = "Balanced accuracy",
        bty = "l"
    )

    abline(
        h = 0.80,
        lty = 2,
        lwd = 1.0,
        col = "grey35"
    )

    abline(
        h = c(
            0.5,
            1.0
        ),
        lty = 3,
        lwd = 0.5,
        col = "grey85"
    )

    valid_ba <- is.finite(
        s$plot_BA
    )

    if (any(valid_ba)) {

        segments(
            x0 = x_pos[
                valid_ba
            ],
            y0 = 0.5,
            x1 = x_pos[
                valid_ba
            ],
            y1 = s$plot_BA[
                valid_ba
            ],
            col = "grey50",
            lwd = 0.8
        )

        points(
            x_pos[
                valid_ba
            ],
            s$plot_BA[
                valid_ba
            ],
            pch = 16,
            cex = 0.85
        )
    }

    title(
        main = sprintf(
            "%s — full comparative panel (%d frozen SNPs; %d individuals)",
            display_names[[module]],
            nrow(s),
            nrow(sample_df)
        ),
        line = 1.0,
        font.main = 2,
        cex.main = 1.1
    )

    text(
        nrow(s) + 0.45,
        0.80,
        "BA = 0.80",
        pos = 2,
        cex = 0.70
    )


    # ========================================================
    # PANEL B — GENOTYPE MATRIX
    # ========================================================

    par(
        mar = c(
            5.0,
            18.5,
            0.5,
            4.5
        )
    )

    plot(
        NA,
        xlim = c(
            0.5,
            nrow(s) + 0.5
        ),
        ylim = c(
            nrow(sample_df) + 0.5,
            0.5
        ),
        xaxt = "n",
        yaxt = "n",
        xlab = "Frozen comparative-selected SNPs in genomic order",
        ylab = "",
        bty = "n"
    )

    par(
        xpd = NA
    )


    # --------------------------------------------------------
    # Sample labels
    #
    # Real y-axis labels in the left margin.
    # Nothing is drawn over genotype cells.
    # --------------------------------------------------------

    axis(
        side = 2,
        at = seq_len(
            nrow(sample_df)
        ),
        labels = sample_df$plot_label,
        las = 1,
        tick = FALSE,
        line = -0.15,
        cex.axis = 0.42
    )


    # --------------------------------------------------------
    # Phenotype strip
    # --------------------------------------------------------

    # Scale strip coordinates to SNP-panel width so it has
    # approximately comparable physical width among modules.
    x_span <- nrow(s)

    strip_right <- (
        0.5 -
        0.011 * x_span
    )

    strip_left <- (
        0.5 -
        0.026 * x_span
    )

    text(
        mean(
            c(
                strip_left,
                strip_right
            )
        ),
        0.12,
        "Phen.",
        cex = 0.65,
        font = 2
    )

    for (
        i in seq_len(
            nrow(sample_df)
        )
    ) {

        phenotype <- as.character(
            sample_df$phenotype_label[i]
        )

        if (
            !phenotype %in%
            names(pheno_cols)
        ) {
            phenotype <- "absent"
        }

        rect(
            strip_left,
            i - 0.5,
            strip_right,
            i + 0.5,
            col = pheno_cols[[phenotype]],
            border = "grey70",
            lwd = 0.20
        )
    }


    # --------------------------------------------------------
    # Divider between phenotype groups
    # --------------------------------------------------------

    phen <- as.character(
        sample_df$phenotype_label
    )

    if (length(phen) > 1) {

        breaks <- which(
            phen[-1]
            !=
            phen[
                -length(phen)
            ]
        )

        if (length(breaks) > 0) {

            for (b in breaks) {

                segments(
                    0.5,
                    b + 0.5,
                    nrow(s) + 0.5,
                    b + 0.5,
                    col = "black",
                    lwd = 1.4
                )
            }
        }
    }


    # --------------------------------------------------------
    # Light column guides
    # --------------------------------------------------------

    for (x in x_pos) {

        segments(
            x + 0.5,
            0.5,
            x + 0.5,
            nrow(sample_df) + 0.5,
            col = "grey94",
            lwd = 0.22
        )
    }


    # --------------------------------------------------------
    # Genotype tiles
    # --------------------------------------------------------

    for (
        j in seq_len(
            nrow(s)
        )
    ) {

        site_key <- s$site_key[j]

        dj <- d[
            d$site_key == site_key,
            ,
            drop = FALSE
        ]

        for (
            i in seq_len(
                nrow(sample_df)
            )
        ) {

            sample_i <- sample_df$sample[i]

            hit <- dj[
                dj$sample == sample_i,
                ,
                drop = FALSE
            ]

            if (nrow(hit) == 1) {

                gt <- as.character(
                    hit$biological_genotype
                )

            } else {

                gt <- "missing"
            }

            if (
                length(gt) != 1 ||
                is.na(gt) ||
                !gt %in% names(gt_cols)
            ) {
                gt <- "missing"
            }

            rect(
                j - 0.5,
                i - 0.5,
                j + 0.5,
                i + 0.5,
                col = gt_cols[[gt]],
                border = "grey80",
                lwd = 0.12
            )
        }
    }


    # --------------------------------------------------------
    # SNP positions
    # --------------------------------------------------------

    axis(
        side = 1,
        at = x_pos,
        labels = thin_axis_labels(
            s$pos
        ),
        las = 2,
        tick = FALSE,
        line = -0.4,
        cex.axis = 0.62
    )


    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    legend(
        "topright",
        inset = c(
            0,
            -0.005
        ),
        horiz = TRUE,
        bty = "n",
        xpd = NA,
        legend = c(
            "AA",
            "PA",
            "PP",
            "missing"
        ),
        fill = gt_cols[
            c(
                "AA",
                "PA",
                "PP",
                "missing"
            )
        ],
        border = "grey55",
        cex = 0.72
    )

    # After the genotype panel is complete the next plot()
    # automatically starts a new PDF page and reuses the layout.
}

dev.off()

cat(
    "\nOutput:\n  ",
    pdf_file,
    "\n",
    sep = ""
)
