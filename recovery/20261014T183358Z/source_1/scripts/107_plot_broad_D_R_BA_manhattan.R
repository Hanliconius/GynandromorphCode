#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "drop_problematic_reanalysis",
    "checkpoints"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "drop_problematic_reanalysis",
    "broad_D_R_BA_manhattan"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

input_files <- c(
    Dennis = file.path(
        checkpoint_dir,
        paste0(
            "longDennis_original502_allSNP_scores.",
            "filtered.WGS_vs_WGSplusCS.tsv"
        )
    ),
    Ray = file.path(
        checkpoint_dir,
        paste0(
            "longRay_original502_allSNP_scores.",
            "filtered.WGS_vs_WGSplusCS.tsv"
        )
    )
)

perfect_file <- file.path(
    checkpoint_dir,
    "perfect_D_R_comparative_SNPs.DP2_GQ20.tsv"
)

missing_inputs <- c(
    input_files[
        !file.exists(
            input_files
        )
    ],
    if (!file.exists(perfect_file)) {
        perfect_file
    } else {
        character(0)
    }
)

if (length(missing_inputs) > 0) {
    stop(
        paste(
            "Missing required files:",
            paste(
                missing_inputs,
                collapse = "\n"
            ),
            sep = "\n"
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
})


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


read_module_scores <- function(
    module_name,
    filename
) {
    d <- read.delim(
        filename,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )

    required <- c(
        "analysis_cohort",
        "scheme",
        "chrom",
        "pos",
        "ref",
        "alt",
        "balanced_accuracy",
        "score_eligible"
    )

    missing_columns <- setdiff(
        required,
        names(d)
    )

    if (length(missing_columns) > 0) {
        stop(
            paste(
                basename(filename),
                "is missing columns:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }

    d <- d[
        d$analysis_cohort == "WGS_plus_CS"
        & d$scheme == "DP2_GQ20",
        ,
        drop = FALSE
    ]

    d$module <- module_name
    d$pos <- as.numeric(d$pos)
    d$balanced_accuracy <- as.numeric(
        d$balanced_accuracy
    )
    d$score_eligible_bool <- as_bool(
        d$score_eligible
    )

    d <- d[
        is.finite(d$pos)
        & is.finite(d$balanced_accuracy),
        ,
        drop = FALSE
    ]

    d$site_key <- paste(
        d$chrom,
        d$pos,
        d$ref,
        d$alt,
        sep = ":"
    )

    d
}


scores <- rbind(
    read_module_scores(
        "Dennis",
        input_files[["Dennis"]]
    ),
    read_module_scores(
        "Ray",
        input_files[["Ray"]]
    )
)

scores$module <- factor(
    scores$module,
    levels = c(
        "Dennis",
        "Ray"
    )
)

scores$position_Mb <- (
    scores$pos / 1000000
)


perfect <- read.delim(
    perfect_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

required_perfect <- c(
    "module",
    "chrom",
    "pos",
    "ref",
    "alt"
)

missing_perfect_columns <- setdiff(
    required_perfect,
    names(perfect)
)

if (length(missing_perfect_columns) > 0) {
    stop(
        paste(
            "Perfect-site table is missing columns:",
            paste(
                missing_perfect_columns,
                collapse = ", "
            )
        )
    )
}

perfect$pos <- as.numeric(
    perfect$pos
)

perfect$site_key <- paste(
    perfect$chrom,
    perfect$pos,
    perfect$ref,
    perfect$alt,
    sep = ":"
)

perfect_keys <- unique(
    paste(
        perfect$module,
        perfect$site_key,
        sep = "|"
    )
)

scores$perfect_selected <- (
    paste(
        as.character(scores$module),
        scores$site_key,
        sep = "|"
    )
    %in% perfect_keys
)


plot_intervals <- data.frame(
    module = factor(
        c(
            "Dennis",
            "Dennis",
            "Ray",
            "Ray"
        ),
        levels = c(
            "Dennis",
            "Ray"
        )
    ),
    interval = c(
        "Dennis1",
        "Dennis2",
        "Ray1",
        "Ray2"
    ),
    start = c(
        1215210,
        1218557,
        1156746,
        1174547
    ),
    end = c(
        1216124,
        1221494,
        1160220,
        1187833
    )
)

plot_intervals$xmin <- (
    plot_intervals$start / 1000000
)

plot_intervals$xmax <- (
    plot_intervals$end / 1000000
)

plot_intervals$xmid <- (
    (
        plot_intervals$start
        + plot_intervals$end
    )
    / 2000000
)


module_colours <- c(
    Dennis = "#E69F00",
    Ray = "#C51B7D"
)

eligible_scores <- scores[
    scores$score_eligible_bool,
    ,
    drop = FALSE
]

ineligible_scores <- scores[
    !scores$score_eligible_bool,
    ,
    drop = FALSE
]

perfect_scores <- scores[
    scores$perfect_selected,
    ,
    drop = FALSE
]


plot <- ggplot() +

    geom_rect(
        data = plot_intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.12
    ) +

    geom_point(
        data = ineligible_scores,
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        colour = "#BDBDBD",
        alpha = 0.28,
        size = 0.38,
        stroke = 0
    ) +

    geom_point(
        data = eligible_scores,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = module
        ),
        alpha = 0.60,
        size = 0.48,
        stroke = 0
    ) +

    geom_point(
        data = perfect_scores,
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        shape = 21,
        fill = "white",
        colour = "#111111",
        size = 1.75,
        stroke = 0.55
    ) +

    geom_hline(
        yintercept = 0.50,
        colour = "#666666",
        linewidth = 0.35
    ) +

    geom_hline(
        yintercept = 0.80,
        colour = "#222222",
        linewidth = 0.45,
        linetype = "dashed"
    ) +

    geom_hline(
        yintercept = 0.90,
        colour = "#777777",
        linewidth = 0.3,
        linetype = "dotted"
    ) +

    geom_text(
        data = plot_intervals,
        aes(
            x = xmid,
            y = 1.018,
            label = interval
        ),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold",
        colour = "#444444"
    ) +

    facet_grid(
        rows = vars(module)
    ) +

    scale_colour_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = "Position across the broad Optix interval (Mb)",
        limits = range(
            scores$position_Mb,
            finite = TRUE
        ),
        expand = expansion(
            mult = c(
                0.005,
                0.005
            )
        )
    ) +

    scale_y_continuous(
        name = "Balanced accuracy",
        breaks = c(
            0.50,
            0.60,
            0.70,
            0.80,
            0.90,
            1.00
        )
    ) +

    coord_cartesian(
        ylim = c(
            0.48,
            1.035
        ),
        clip = "off"
    ) +

    labs(
        title = paste(
            "Comparative genotype-phenotype association",
            "across the broad Optix interval"
        ),
        subtitle = paste(
            "Final cleaned WGS+CS cohort; DP >= 2 and GQ >= 20.",
            "Black rings mark the perfect SNPs used in the genotype plots."
        ),
        caption = paste(
            "Colored points are score-eligible;",
            "light-gray points have a defined BA but failed",
            "score-eligibility requirements.",
            "Dashed line: BA = 0.80.",
            "The gynandromorph was not used in these calculations."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9.5,
            colour = "#444444",
            margin = margin(
                b = 8
            )
        ),
        plot.caption = element_text(
            size = 8,
            colour = "#555555",
            hjust = 0
        ),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#EEEEEE",
            linewidth = 0.3
        ),
        panel.grid.major.y = element_line(
            colour = "#E2E2E2",
            linewidth = 0.35
        ),
        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = NA
        ),
        strip.text.y = element_text(
            face = "bold",
            size = 10
        ),
        panel.spacing.y = grid::unit(
            0.35,
            "lines"
        ),
        plot.margin = margin(
            10,
            18,
            10,
            10
        )
    )


pdf_file <- file.path(
    figure_dir,
    "broad_D_R_DP2_GQ20_BA_manhattan.pdf"
)

png_file <- file.path(
    figure_dir,
    "broad_D_R_DP2_GQ20_BA_manhattan.png"
)

ggsave(
    filename = pdf_file,
    plot = plot,
    width = 11.5,
    height = 7.2,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = png_file,
    plot = plot,
    width = 11.5,
    height = 7.2,
    units = "in",
    dpi = 450,
    type = "cairo",
    limitsize = FALSE
)


summary_rows <- do.call(
    rbind,
    lapply(
        levels(scores$module),
        function(module_name) {
            x <- scores[
                scores$module == module_name,
                ,
                drop = FALSE
            ]

            data.frame(
                module = module_name,
                SNPs_with_defined_BA = nrow(x),
                score_eligible = sum(
                    x$score_eligible_bool
                ),
                BA_ge_0_80 = sum(
                    x$score_eligible_bool
                    & x$balanced_accuracy >= 0.80
                ),
                BA_ge_0_90 = sum(
                    x$score_eligible_bool
                    & x$balanced_accuracy >= 0.90
                ),
                BA_equal_1 = sum(
                    x$score_eligible_bool
                    & abs(
                        x$balanced_accuracy - 1
                    ) < 1e-12
                ),
                highlighted_perfect_SNPs = sum(
                    x$perfect_selected
                )
            )
        }
    )
)

summary_file <- file.path(
    checkpoint_dir,
    "broad_D_R_DP2_GQ20_BA_manhattan.summary.tsv"
)

write.table(
    summary_rows,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "BROAD D/R BALANCED-ACCURACY MANHATTAN PLOT\n\n"
)

print(
    summary_rows,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
