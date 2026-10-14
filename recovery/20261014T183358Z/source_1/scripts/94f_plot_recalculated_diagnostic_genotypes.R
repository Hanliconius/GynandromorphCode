#!/usr/bin/env Rscript

project <- "/scratch/martinlab/joe/GynandroCalling"

ck <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
)

infile <- file.path(
    ck,
    "original502_recalculated_diagnostic_genotypes.WGS_vs_WGSplusCS.long.tsv"
)

outdir <- file.path(
    project,
    "analysis",
    "figures",
    "supplement_optix_predictive_genotype_plots"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

d <- read.delim(
    infile,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

modules <- c(
    "Band",
    "Ray2",
    "Dennis2"
)

cohorts <- c(
    "WGS_only",
    "WGS_plus_CS"
)

module_names <- c(
    Band = "Band",
    Ray2 = "Ray 2",
    Dennis2 = "Dennis 2"
)

cohort_names <- c(
    WGS_only = "WGS only",
    WGS_plus_CS = "WGS + selective sequencing"
)

gt_cols <- c(
    AA = "#2F6BFF",
    PA = "#F39C12",
    PP = "#D81BCE",
    missing = "#E5E7EB"
)

type_cols <- c(
    WGS = "#333333",
    CS = "#BDBDBD"
)


as_bool <- function(x) {

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


for (cohort in cohorts) {

    dc <- d[
        d$cohort == cohort,
        ,
        drop = FALSE
    ]

    max_n <- max(
        sapply(
            split(
                dc$sample,
                dc$module
            ),
            function(x) {
                length(
                    unique(x)
                )
            }
        )
    )

    page_height <- max(
        20,
        8 + max_n * 0.115
    )

    outfile <- file.path(
        outdir,
        paste0(
            "optix_recalculated_diagnostic_",
            cohort,
            "_Band_Ray2_Dennis2.pdf"
        )
    )

    pdf(
        outfile,
        width = 17,
        height = page_height,
        onefile = TRUE,
        useDingbats = FALSE
    )

    layout(
        matrix(
            c(
                1,
                2
            ),
            ncol = 1
        ),
        heights = c(
            1.25,
            16
        )
    )

    for (module in modules) {

        dm <- dc[
            dc$module == module,
            ,
            drop = FALSE
        ]

        if (nrow(dm) == 0) {
            next
        }

        sites <- unique(
            dm[
                ,
                c(
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "site_order",
                    "balanced_accuracy",
                    "passes_BA_0.80",
                    "selection_rule"
                )
            ]
        )

        sites$site_order <- as.integer(
            sites$site_order
        )

        sites$balanced_accuracy <- as.numeric(
            sites$balanced_accuracy
        )

        sites$passes_BA_0.80 <- as_bool(
            sites$passes_BA_0.80
        )

        sites <- sites[
            order(
                sites$site_order
            ),
            ,
            drop = FALSE
        ]

        sample_df <- unique(
            dm[
                ,
                c(
                    "sample",
                    "taxon_name",
                    "data_type",
                    "phenotype"
                )
            ]
        )

        sample_df$phenotype <- as.numeric(
            sample_df$phenotype
        )

        sample_df <- sample_df[
            order(
                sample_df$phenotype,
                sample_df$data_type,
                sample_df$taxon_name,
                sample_df$sample
            ),
            ,
            drop = FALSE
        ]

        rownames(
            sample_df
        ) <- NULL

        nsites <- nrow(
            sites
        )

        nsamples <- nrow(
            sample_df
        )

        xpos <- seq_len(
            nsites
        )


        # ====================================================
        # BA track
        # ====================================================

        par(
            mar = c(
                0.8,
                19.5,
                3.2,
                4.7
            )
        )

        plot(
            NA,
            xlim = c(
                0.5,
                nsites + 0.5
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
            lwd = 1
        )

        segments(
            xpos,
            0.5,
            xpos,
            sites$balanced_accuracy,
            col = "grey60",
            lwd = 0.8
        )

        points(
            xpos[
                sites$passes_BA_0.80
            ],
            sites$balanced_accuracy[
                sites$passes_BA_0.80
            ],
            pch = 16,
            cex = 0.8
        )

        if (
            any(
                !sites$passes_BA_0.80
            )
        ) {

            points(
                xpos[
                    !sites$passes_BA_0.80
                ],
                sites$balanced_accuracy[
                    !sites$passes_BA_0.80
                ],
                pch = 1,
                cex = 0.85
            )
        }

        title_text <- sprintf(
            "%s — %s (%d recalculated SNPs; %d individuals)",
            module_names[module],
            cohort_names[cohort],
            nsites,
            nsamples
        )

        if (
            cohort == "WGS_plus_CS"
            &&
            module == "Band"
            &&
            !any(
                sites$passes_BA_0.80
            )
        ) {

            title_text <- paste0(
                title_text,
                "\nNo eligible Band SNP reached BA = 0.80; ",
                "all eligible expected-orientation SNPs shown"
            )
        }

        title(
            main = title_text,
            line = 0.8,
            cex.main = 1.0
        )


        # ====================================================
        # Genotype matrix
        # ====================================================

        par(
            mar = c(
                5,
                19.5,
                0.5,
                4.7
            )
        )

        plot(
            NA,
            xlim = c(
                0.5,
                nsites + 0.5
            ),
            ylim = c(
                nsamples + 0.5,
                0.5
            ),
            xaxt = "n",
            yaxt = "n",
            xlab = "Recalculated diagnostic SNPs in genomic order",
            ylab = "",
            bty = "n"
        )

        label_cex <- if (
            nsamples > 300
        ) {
            0.38
        } else if (
            nsamples > 180
        ) {
            0.43
        } else {
            0.54
        }

        axis(
            side = 2,
            at = seq_len(
                nsamples
            ),
            labels = sample_df$taxon_name,
            las = 1,
            tick = FALSE,
            line = -0.15,
            cex.axis = label_cex
        )


        # Sequencing-type strip.
        x_span <- nsites

        type_right <- 0.5 - 0.010 * x_span
        type_left <- 0.5 - 0.028 * x_span

        par(
            xpd = NA
        )

        text(
            mean(
                c(
                    type_left,
                    type_right
                )
            ),
            0.1,
            "Type",
            cex = 0.58,
            font = 2
        )

        for (i in seq_len(nsamples)) {

            dtype <- as.character(
                sample_df$data_type[i]
            )

            if (
                !dtype %in% names(type_cols)
            ) {
                dtype <- "CS"
            }

            rect(
                type_left,
                i - 0.5,
                type_right,
                i + 0.5,
                col = type_cols[dtype],
                border = "grey65",
                lwd = 0.2
            )
        }


        # Phenotype-class boundary.
        phenotype_vec <- sample_df$phenotype

        if (
            length(
                phenotype_vec
            ) > 1
        ) {

            phenotype_breaks <- which(
                phenotype_vec[-1]
                !=
                phenotype_vec[
                    -length(
                        phenotype_vec
                    )
                ]
            )

            for (b in phenotype_breaks) {

                segments(
                    0.5,
                    b + 0.5,
                    nsites + 0.5,
                    b + 0.5,
                    lwd = 1.4
                )
            }
        }


        # Genotypes.
        for (j in seq_len(nsites)) {

            site_pos <- sites$pos[j]

            dj <- dm[
                dm$pos == site_pos,
                ,
                drop = FALSE
            ]

            for (i in seq_len(nsamples)) {

                hit <- dj[
                    dj$sample
                    ==
                    sample_df$sample[i],
                    ,
                    drop = FALSE
                ]

                gt <- "missing"

                if (
                    nrow(hit) == 1
                ) {

                    value <- as.character(
                        hit$biological_genotype[1]
                    )

                    if (
                        !is.na(value)
                        &&
                        value %in% names(gt_cols)
                    ) {
                        gt <- value
                    }
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


        # Position labels.
        label_step <- if (
            nsites <= 20
        ) {
            1
        } else if (
            nsites <= 70
        ) {
            4
        } else {
            6
        }

        position_labels <- rep(
            "",
            nsites
        )

        keep_labels <- seq(
            1,
            nsites,
            by = label_step
        )

        position_labels[
            keep_labels
        ] <- sites$pos[
            keep_labels
        ]

        axis(
            side = 1,
            at = xpos,
            labels = position_labels,
            las = 2,
            tick = FALSE,
            line = -0.4,
            cex.axis = 0.60
        )

        legend(
            "topright",
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

        legend(
            "bottomright",
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
        "Created: ",
        outfile,
        "\n",
        sep = ""
    )
}
