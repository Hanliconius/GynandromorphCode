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

cohort_file <- file.path(
    ck,
    "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

scores_file <- file.path(
    ck,
    "original502_frozen_site_scores.WGS_vs_WGSplusCS.tsv"
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


# ============================================================
# Read inputs
# ============================================================

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

cohorts <- read.delim(
    cohort_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

scores <- read.delim(
    scores_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


# ============================================================
# Constants
# ============================================================

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

cohort_order <- c(
    "WGS_only",
    "WGS_plus_CS"
)

cohort_display <- c(
    WGS_only = "WGS only",
    WGS_plus_CS = "WGS + selective sequencing"
)


# Genotype colours in biological P/A space.
gt_cols <- c(
    AA = "#2F6BFF",
    PA = "#F39C12",
    PP = "#D81BCE",
    missing = "#E5E7EB"
)

# Phenotype annotation.
phen_cols <- c(
    absent = "#FFFFFF",
    present = "#111111"
)

# Neutral sequencing-type annotation.
type_cols <- c(
    WGS = "#333333",
    CS = "#BDBDBD"
)


# ============================================================
# Helpers
# ============================================================

clean_label <- function(x) {

    x <- as.character(x)

    x <- gsub(
        "[[:space:]]+",
        " ",
        x
    )

    trimws(x)
}


to_logical <- function(x) {

    tolower(
        trimws(
            as.character(x)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
    )
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


# ============================================================
# Prepare fields
# ============================================================

cohorts$keep_bool <- to_logical(
    cohorts$keep
)

cohorts$callable_fraction <- suppressWarnings(
    as.numeric(
        cohorts$callable_fraction
    )
)

scores$balanced_accuracy <- suppressWarnings(
    as.numeric(
        scores$balanced_accuracy
    )
)

sites$pos <- as.integer(
    sites$pos
)


# ============================================================
# Plot one cohort to one multipage PDF
# ============================================================

plot_cohort <- function(cohort_name) {

    cm <- cohorts[
        cohorts$cohort == cohort_name &
        cohorts$keep_bool,
        ,
        drop = FALSE
    ]

    if (nrow(cm) == 0) {
        stop(
            paste(
                "No retained samples for",
                cohort_name
            )
        )
    }

    module_n <- tapply(
        cm$sample,
        cm$module,
        function(x) {
            length(
                unique(x)
            )
        }
    )

    max_n <- max(
        module_n,
        na.rm = TRUE
    )

    # WGS-only needs a moderately tall page.
    # WGS+CS needs a genuinely tall vector page for ~420 rows.
    page_height <- max(
        20,
        8 + max_n * 0.115
    )

    pdf_file <- file.path(
        out_dir,
        paste0(
            "optix_filtered_",
            cohort_name,
            "_Band_Ray2_Dennis2.pdf"
        )
    )

    pdf(
        pdf_file,
        width = 17,
        height = page_height,
        onefile = TRUE,
        useDingbats = FALSE
    )

    # One page/module:
    # upper BA panel + lower genotype matrix.
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
            16
        )
    )


    for (module in module_order) {

        # ----------------------------------------------------
        # Retained samples for this particular module.
        # ----------------------------------------------------

        keep_rows <- cm[
            cm$module == module,
            ,
            drop = FALSE
        ]

        keep_samples <- unique(
            keep_rows$sample
        )

        d <- long[
            long$module == module &
            long$sample %in% keep_samples,
            ,
            drop = FALSE
        ]

        s <- sites[
            sites$module == module,
            ,
            drop = FALSE
        ]

        sc <- scores[
            scores$cohort == cohort_name &
            scores$module == module,
            ,
            drop = FALSE
        ]

        if (
            nrow(d) == 0 ||
            nrow(s) == 0
        ) {
            next
        }


        # ----------------------------------------------------
        # Frozen-site order remains unchanged.
        # ----------------------------------------------------

        s <- s[
            order(
                s$site_order_within_module,
                s$pos
            ),
            ,
            drop = FALSE
        ]

        # Attach filtered-cohort BA to the frozen site table.
        ba_lookup <- sc[
            ,
            c(
                "site_key",
                "balanced_accuracy",
                "source_status"
            ),
            drop = FALSE
        ]

        names(
            ba_lookup
        )[names(ba_lookup) == "balanced_accuracy"] <- "filtered_BA"

        names(
            ba_lookup
        )[names(ba_lookup) == "source_status"] <- "filtered_source_status"

        s <- merge(
            s,
            ba_lookup,
            by = "site_key",
            all.x = TRUE,
            sort = FALSE
        )

        # Restore exact frozen-site order after merge.
        s <- s[
            order(
                s$site_order_within_module
            ),
            ,
            drop = FALSE
        ]


        # ----------------------------------------------------
        # Build filtered sample annotation table.
        # ----------------------------------------------------

        sample_df <- unique(
            d[
                ,
                c(
                    "sample",
                    "taxon_name",
                    "phenotype_label",
                    "phenotype_value"
                ),
                drop = FALSE
            ]
        )

        sample_df <- merge(
            sample_df,
            keep_rows[
                ,
                c(
                    "sample",
                    "data_type",
                    "callable_fraction"
                ),
                drop = FALSE
            ],
            by = "sample",
            all.x = TRUE,
            sort = FALSE
        )

        sample_df$plot_label <- clean_label(
            sample_df$taxon_name
        )

        sample_df$phenotype_value <- as.numeric(
            sample_df$phenotype_value
        )

        # Sort by phenotype first.
        #
        # For WGS+CS, then sort by data type so we can visually
        # compare WGS and CS blocks within phenotype classes.
        #
        # Finally taxon and accession provide deterministic order,
        # although accession itself is never printed.
        sample_df <- sample_df[
            order(
                sample_df$phenotype_value,
                sample_df$taxon_name,
                sample_df$sample
            ),
            ,
            drop = FALSE
        ]

        rownames(
            sample_df
        ) <- NULL

        sample_df$plot_order <- seq_len(
            nrow(
                sample_df
            )
        )

        x_pos <- seq_len(
            nrow(s)
        )


        # ====================================================
        # TOP: filtered-cohort balanced accuracy
        # ====================================================

        par(
            mar = c(
                0.8,
                19.5,
                3.0,
                4.7
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
            col = "grey86"
        )

        valid <- is.finite(
            s$filtered_BA
        )

        if (any(valid)) {

            segments(
                x0 = x_pos[valid],
                y0 = 0.5,
                x1 = x_pos[valid],
                y1 = s$filtered_BA[valid],
                col = "grey55",
                lwd = 0.8
            )

            points(
                x_pos[valid],
                s$filtered_BA[valid],
                pch = 16,
                cex = 0.82
            )
        }

        title(
            main = sprintf(
                "%s — %s (%d frozen SNPs; %d retained individuals)",
                display_names[module],
                cohort_display[cohort_name],
                nrow(s),
                nrow(sample_df)
            ),
            line = 1.0,
            font.main = 2,
            cex.main = 1.05
        )

        text(
            nrow(s) + 0.45,
            0.80,
            "BA = 0.80",
            pos = 2,
            cex = 0.68
        )


        # ====================================================
        # BOTTOM: genotype matrix
        # ====================================================

        par(
            mar = c(
                5.0,
                19.5,
                0.5,
                4.7
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


        # ----------------------------------------------------
        # Taxon names: genuine y-axis labels.
        #
        # No accessions are plotted and no taxon text is placed
        # inside the genotype matrix.
        # ----------------------------------------------------

        label_cex <- if (
            nrow(sample_df) > 300
        ) {
            0.38
        } else if (
            nrow(sample_df) > 180
        ) {
            0.43
        } else {
            0.54
        }

        axis(
            side = 2,
            at = seq_len(
                nrow(sample_df)
            ),
            labels = sample_df$plot_label,
            las = 1,
            tick = FALSE,
            line = -0.15,
            cex.axis = label_cex
        )


        # ----------------------------------------------------
        # Annotation strips just outside first SNP.
        # ----------------------------------------------------

        x_span <- nrow(s)

        phen_right <- 0.5 - 0.010 * x_span
        phen_left  <- 0.5 - 0.025 * x_span

        type_right <- 0.5 - 0.031 * x_span
        type_left  <- 0.5 - 0.046 * x_span

        text(
            mean(
                c(
                    phen_left,
                    phen_right
                )
            ),
            0.10,
            "Phen.",
            cex = 0.58,
            font = 2
        )

        text(
            mean(
                c(
                    type_left,
                    type_right
                )
            ),
            0.10,
            "Type",
            cex = 0.58,
            font = 2
        )


        # ----------------------------------------------------
        # Annotation strip rows.
        # ----------------------------------------------------

        for (
            i in seq_len(
                nrow(sample_df)
            )
        ) {

            phen <- as.character(
                sample_df$phenotype_label[i]
            )

            if (
                !phen %in% names(phen_cols)
            ) {
                phen <- "absent"
            }

            dtype <- as.character(
                sample_df$data_type[i]
            )

            if (
                !dtype %in% names(type_cols)
            ) {
                dtype <- "CS"
            }

            rect(
                phen_left,
                i - 0.5,
                phen_right,
                i + 0.5,
                col = phen_cols[phen],
                border = "grey65",
                lwd = 0.20
            )

            rect(
                type_left,
                i - 0.5,
                type_right,
                i + 0.5,
                col = type_cols[dtype],
                border = "grey65",
                lwd = 0.20
            )
        }


        # ----------------------------------------------------
        # Phenotype separators.
        # ----------------------------------------------------

        phen_vec <- as.character(
            sample_df$phenotype_label
        )

        if (
            length(
                phen_vec
            ) > 1
        ) {

            breaks <- which(
                phen_vec[-1] !=
                phen_vec[-length(phen_vec)]
            )

            if (
                length(
                    breaks
                ) > 0
            ) {

                for (
                    b in breaks
                ) {

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


        # ----------------------------------------------------
        # Sequencing-type separators within phenotype groups.
        # ----------------------------------------------------

        type_vec <- as.character(
            sample_df$data_type
        )

        if (
            length(
                type_vec
            ) > 1
        ) {

            type_breaks <- which(
                type_vec[-1] !=
                type_vec[-length(type_vec)]
            )

            if (
                length(
                    type_breaks
                ) > 0
            ) {

                for (
                    b in type_breaks
                ) {

                    # Don't overdraw a phenotype boundary.
                    is_phen_break <- (
                        length(phen_vec) > 1 &&
                        phen_vec[b] != phen_vec[b + 1]
                    )

                    if (
                        !is_phen_break
                    ) {

                        segments(
                            0.5,
                            b + 0.5,
                            nrow(s) + 0.5,
                            b + 0.5,
                            col = "grey55",
                            lwd = 0.55
                        )
                    }
                }
            }
        }


        # ----------------------------------------------------
        # Genotype matrix.
        # ----------------------------------------------------

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

                if (
                    nrow(hit) == 1
                ) {

                    gt <- as.character(
                        hit$biological_genotype[1]
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
                    col = gt_cols[gt],
                    border = "grey80",
                    lwd = 0.12
                )
            }
        }


        # ----------------------------------------------------
        # SNP position labels.
        # ----------------------------------------------------

        axis(
            side = 1,
            at = x_pos,
            labels = thin_axis_labels(
                s$pos
            ),
            las = 2,
            tick = FALSE,
            line = -0.4,
            cex.axis = 0.61
        )


        # ----------------------------------------------------
        # Genotype legend.
        # ----------------------------------------------------

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
            cex = 0.70
        )


        # Sequencing-type key.
        legend(
            "bottomright",
            inset = c(
                0,
                -0.005
            ),
            horiz = TRUE,
            bty = "n",
            xpd = NA,
            legend = c(
                "WGS",
                "CS"
            ),
            fill = type_cols[
                c(
                    "WGS",
                    "CS"
                )
            ],
            border = "grey55",
            cex = 0.65
        )
    }

    dev.off()

    cat(
        "\nCreated:\n  ",
        pdf_file,
        "\n",
        sep = ""
    )
}


# ============================================================
# Produce both PDFs.
# ============================================================

for (
    cohort_name in cohort_order
) {

    plot_cohort(
        cohort_name
    )
}

cat(
    "\nDone.\n"
)
