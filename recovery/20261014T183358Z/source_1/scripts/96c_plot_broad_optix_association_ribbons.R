#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "original_502",
    "checkpoints"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "broad_optix_association_ribbons"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

inputs <- c(
    Band = file.path(
        checkpoint_dir,
        paste0(
            "longBand_original502_allSNP_scores.filtered.",
            "WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
        )
    ),
    Ray = file.path(
        checkpoint_dir,
        paste0(
            "longRay_original502_allSNP_scores.filtered.",
            "WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
        )
    ),
    Dennis = file.path(
        checkpoint_dir,
        paste0(
            "longDennis_original502_allSNP_scores.filtered.",
            "WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
        )
    )
)

selected_output <- file.path(
    checkpoint_dir,
    "broad_optix_strict_comparative_SNPs.tsv"
)

summary_output <- file.path(
    checkpoint_dir,
    "broad_optix_strict_comparative_SNPs.summary.tsv"
)

pdf_output <- file.path(
    figure_dir,
    "broad_optix_comparative_association_ribbons.pdf"
)

svg_output <- file.path(
    figure_dir,
    "broad_optix_comparative_association_ribbons.svg"
)


as_boolean <- function(x) {

    tolower(
        trimws(
            as.character(x)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
}


read_module <- function(module_name, filename) {

    d <- read.delim(
        filename,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    required <- c(
        "analysis_cohort",
        "chrom",
        "pos",
        "ref",
        "alt",
        "all_schemes_score_eligible",
        "all_schemes_match_expected_orientation",
        "all_schemes_expected_pattern_eligible",
        "minimum_BA_across_schemes"
    )

    missing_columns <- setdiff(
        required,
        names(d)
    )

    if (length(missing_columns) > 0) {
        stop(
            paste0(
                module_name,
                " table is missing columns: ",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }

    d$plot_module <- module_name
    d$pos <- as.numeric(d$pos)

    d$minimum_BA_across_schemes <- as.numeric(
        d$minimum_BA_across_schemes
    )

    d$eligible_all_schemes <- (
        as_boolean(
            d$all_schemes_score_eligible
        )
        &
        as_boolean(
            d$all_schemes_match_expected_orientation
        )
        &
        as_boolean(
            d$all_schemes_expected_pattern_eligible
        )
    )

    d <- d[
        d$eligible_all_schemes
        &
        is.finite(
            d$minimum_BA_across_schemes
        ),
        ,
        drop = FALSE
    ]

    d
}


missing_inputs <- inputs[
    !file.exists(inputs)
]

if (length(missing_inputs) > 0) {

    cat(
        "The ribbon plot was not generated because these inputs ",
        "are still missing:\n",
        sep = ""
    )

    for (filename in missing_inputs) {
        cat(
            filename,
            "\n"
        )
    }

} else {

    tables <- lapply(
        names(inputs),
        function(module_name) {
            read_module(
                module_name,
                inputs[[module_name]]
            )
        }
    )

    combined <- do.call(
        rbind,
        tables
    )

    site_key <- paste(
        combined$plot_module,
        combined$chrom,
        combined$pos,
        combined$ref,
        combined$alt,
        sep = "\t"
    )

    site_groups <- split(
        combined,
        site_key
    )

    collapsed_rows <- lapply(
        site_groups,
        function(site_rows) {

            cohorts <- unique(
                site_rows$analysis_cohort
            )

            required_cohorts <- c(
                "WGS_only",
                "WGS_plus_CS"
            )

            if (!all(required_cohorts %in% cohorts)) {
                return(NULL)
            }

            site_rows <- site_rows[
                site_rows$analysis_cohort
                %in%
                required_cohorts,
                ,
                drop = FALSE
            ]

            data.frame(
                module = site_rows$plot_module[[1]],
                chrom = site_rows$chrom[[1]],
                pos = site_rows$pos[[1]],
                ref = site_rows$ref[[1]],
                alt = site_rows$alt[[1]],

                WGS_only_minimum_BA = min(
                    site_rows$minimum_BA_across_schemes[
                        site_rows$analysis_cohort
                        ==
                        "WGS_only"
                    ],
                    na.rm = TRUE
                ),

                WGS_plus_CS_minimum_BA = min(
                    site_rows$minimum_BA_across_schemes[
                        site_rows$analysis_cohort
                        ==
                        "WGS_plus_CS"
                    ],
                    na.rm = TRUE
                ),

                minimum_BA_all_analyses = min(
                    site_rows$minimum_BA_across_schemes,
                    na.rm = TRUE
                ),

                stringsAsFactors = FALSE
            )
        }
    )

    collapsed_rows <- collapsed_rows[
        !vapply(
            collapsed_rows,
            is.null,
            logical(1)
        )
    ]

    selected <- do.call(
        rbind,
        collapsed_rows
    )

    selected <- selected[
        selected$minimum_BA_all_analyses
        >=
        0.80,
        ,
        drop = FALSE
    ]

    module_order <- c(
        "Band",
        "Ray",
        "Dennis"
    )

    selected$module <- factor(
        selected$module,
        levels = module_order
    )

    selected <- selected[
        order(
            selected$module,
            selected$pos
        ),
        ,
        drop = FALSE
    ]

    selected$module <- as.character(
        selected$module
    )

    write.table(
        selected,
        selected_output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    summary_table <- do.call(
        rbind,
        lapply(
            module_order,
            function(module_name) {

                x <- selected[
                    selected$module
                    ==
                    module_name,
                    ,
                    drop = FALSE
                ]

                data.frame(
                    module = module_name,
                    selected_SNPs = nrow(x),
                    minimum_BA = if (nrow(x) > 0) {
                        min(
                            x$minimum_BA_all_analyses
                        )
                    } else {
                        NA_real_
                    },
                    median_BA = if (nrow(x) > 0) {
                        median(
                            x$minimum_BA_all_analyses
                        )
                    } else {
                        NA_real_
                    },
                    maximum_BA = if (nrow(x) > 0) {
                        max(
                            x$minimum_BA_all_analyses
                        )
                    } else {
                        NA_real_
                    },
                    stringsAsFactors = FALSE
                )
            }
        )
    )

    write.table(
        summary_table,
        summary_output,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    optix_start <- 938197
    optix_end <- 1239419

    original_intervals <- data.frame(
        label = c(
            "B",
            "R1",
            "R2",
            "D1",
            "D2"
        ),
        module = c(
            "Band",
            "Ray",
            "Ray",
            "Dennis",
            "Dennis"
        ),
        start = c(
            1112651,
            1156746,
            1174547,
            1215210,
            1218557
        ),
        end = c(
            1114314,
            1160220,
            1187833,
            1216124,
            1221494
        ),
        stringsAsFactors = FALSE
    )

    mago <- c(
        Band = "#B21875",
        Ray = "#E68613",
        Dennis = "#2B83BA"
    )

    track_y <- c(
        Band = 2.75,
        Ray = 1.75,
        Dennis = 0.75
    )

    interval_label_y <- c(
        B = 4.22,
        R1 = 4.22,
        R2 = 4.08,
        D1 = 4.22,
        D2 = 4.08
    )

    interval_label_x <- c(
        B = mean(c(1112651, 1114314)),
        R1 = mean(c(1156746, 1160220)),
        R2 = mean(c(1174547, 1187833)),
        D1 = mean(c(1215210, 1216124)) - 3500,
        D2 = mean(c(1218557, 1221494)) + 3500
    )


    draw_ribbon <- function() {

        par(
            mar = c(4.3, 7.2, 2.3, 1.5),
            xpd = NA,
            las = 1
        )

        plot(
            NA,
            xlim = c(
                optix_start,
                optix_end
            ),
            ylim = c(
                0.25,
                4.45
            ),
            axes = FALSE,
            xlab = "",
            ylab = "",
            bty = "n"
        )

        # Very light interval projections through all tracks.
        for (row_index in seq_len(nrow(original_intervals))) {

            interval <- original_intervals[
                row_index,
                ,
                drop = FALSE
            ]

            interval_colour <- mago[
                interval$module
            ]

            rect(
                interval$start,
                0.35,
                interval$end,
                4.02,
                col = adjustcolor(
                    interval_colour,
                    alpha.f = 0.055
                ),
                border = NA
            )
        }

        # Comparative association ribbons.
        for (module_name in module_order) {

            y <- track_y[
                module_name
            ]

            rect(
                optix_start,
                y - 0.23,
                optix_end,
                y + 0.23,
                col = "#F2F2F2",
                border = "#D3D3D3",
                lwd = 0.7
            )

            abline_y <- y

            segments(
                optix_start,
                abline_y,
                optix_end,
                abline_y,
                col = "#BEBEBE",
                lwd = 0.6
            )

            module_sites <- selected[
                selected$module
                ==
                module_name,
                ,
                drop = FALSE
            ]

            if (nrow(module_sites) > 0) {

                scaled_height <- (
                    0.07
                    +
                    0.15
                    *
                    pmin(
                        1,
                        pmax(
                            0,
                            (
                                module_sites$
                                    minimum_BA_all_analyses
                                -
                                0.80
                            )
                            /
                            0.20
                        )
                    )
                )

                segments(
                    module_sites$pos,
                    y - scaled_height,
                    module_sites$pos,
                    y + scaled_height,
                    col = adjustcolor(
                        mago[
                            module_name
                        ],
                        alpha.f = 0.78
                    ),
                    lwd = 0.65,
                    lend = "butt"
                )
            }

            text(
                optix_start
                -
                0.018
                *
                (
                    optix_end
                    -
                    optix_start
                ),
                y,
                labels = paste0(
                    module_name,
                    "  (",
                    nrow(module_sites),
                    ")"
                ),
                adj = 1,
                cex = 0.86,
                font = 2,
                col = mago[
                    module_name
                ]
            )
        }

        # Top annotated Optix interval.
        rect(
            optix_start,
            3.72,
            optix_end,
            4.00,
            col = "#D9D9D9",
            border = "#A6A6A6",
            lwd = 0.8
        )

        for (row_index in seq_len(nrow(original_intervals))) {

            interval <- original_intervals[
                row_index,
                ,
                drop = FALSE
            ]

            interval_colour <- mago[
                interval$module
            ]

            interval_midpoint <- mean(
                c(
                    interval$start,
                    interval$end
                )
            )

            rect(
                interval$start,
                3.72,
                interval$end,
                4.00,
                col = interval_colour,
                border = interval_colour
            )

            label_x <- interval_label_x[
                interval$label
            ]

            label_y <- interval_label_y[
                interval$label
            ]

            segments(
                interval_midpoint,
                4.00,
                label_x,
                label_y - 0.05,
                col = interval_colour,
                lwd = 0.7
            )

            text(
                label_x,
                label_y,
                labels = interval$label,
                cex = 0.80,
                font = 2,
                col = interval_colour
            )
        }

        axis_positions <- pretty(
            c(
                optix_start,
                optix_end
            ),
            n = 7
        )

        axis_positions <- axis_positions[
            axis_positions >= optix_start
            &
            axis_positions <= optix_end
        ]

        axis(
            side = 1,
            at = axis_positions,
            labels = sprintf(
                "%.2f",
                axis_positions / 1000000
            ),
            cex.axis = 0.82,
            tck = -0.018,
            mgp = c(2.4, 0.6, 0)
        )

        mtext(
            "Position on Hmr18 (Mb)",
            side = 1,
            line = 2.5,
            cex = 0.88
        )

        title(
            main = "Comparative genotype–phenotype signal across optix",
            cex.main = 1.02,
            line = 0.7
        )

        mtext(
            paste0(
                "Each vertical mark is an independently selected SNP; ",
                "mark height reflects minimum balanced accuracy"
            ),
            side = 3,
            line = -0.55,
            cex = 0.75,
            col = "#4D4D4D"
        )
    }


    pdf(
        pdf_output,
        width = 10.5,
        height = 4.8,
        useDingbats = FALSE
    )

    draw_ribbon()

    dev.off()


    svg(
        svg_output,
        width = 10.5,
        height = 4.8,
        pointsize = 12
    )

    draw_ribbon()

    dev.off()


    cat(
        "BROAD OPTIX ASSOCIATION RIBBON\n\n"
    )

    print(
        summary_table,
        row.names = FALSE
    )

    cat(
        "\nSelected-site table:\n",
        selected_output,
        "\n\nSummary:\n",
        summary_output,
        "\n\nPDF:\n",
        pdf_output,
        "\n\nSVG:\n",
        svg_output,
        "\n",
        sep = ""
    )
}
