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

required_files <- c(
    input_files,
    perfect_file
)

missing_files <- required_files[
    !file.exists(
        required_files
    )
]

if (length(missing_files) > 0) {
    stop(
        paste(
            "Missing required files:",
            paste(
                missing_files,
                collapse = "\n"
            ),
            sep = "\n"
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
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


canonical_gt <- function(x) {
    x <- gsub(
        "\\|",
        "/",
        trimws(
            as.character(x)
        )
    )

    x[
        x == "1/0"
    ] <- "0/1"

    x[
        !(x %in% c(
            "0/0",
            "0/1",
            "1/1"
        ))
    ] <- NA_character_

    x
}


read_scores <- function(
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
        "score_eligible",
        "TN_Left_GT",
        "TN_Left_AD",
        "TN_Left_DP",
        "TN_Left_GQ",
        "TN_Right_GT",
        "TN_Right_AD",
        "TN_Right_DP",
        "TN_Right_GQ"
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
    d$position_Mb <- d$pos / 1000000

    d$balanced_accuracy <- as.numeric(
        d$balanced_accuracy
    )

    d$score_eligible_bool <- as_bool(
        d$score_eligible
    )

    d$site_key <- paste(
        d$chrom,
        d$pos,
        d$ref,
        d$alt,
        sep = ":"
    )

    d
}


dennis <- read_scores(
    "Dennis",
    input_files[["Dennis"]]
)

ray <- read_scores(
    "Ray",
    input_files[["Ray"]]
)

scores <- rbind(
    dennis,
    ray
)


perfect <- read.delim(
    perfect_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

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

perfect_Dennis <- unique(
    perfect$site_key[
        perfect$module == "Dennis"
    ]
)

perfect_Ray <- unique(
    perfect$site_key[
        perfect$module == "Ray"
    ]
)


scores$perfect_selected <- (
    (
        scores$module == "Dennis"
        & scores$site_key %in% perfect_Dennis
    )
    |
    (
        scores$module == "Ray"
        & scores$site_key %in% perfect_Ray
    )
)


# One physical copy of every SNP for the focal genotype tracks.
focal <- dennis[
    !duplicated(
        dennis$site_key
    ),
    ,
    drop = FALSE
]

focal$left_GT <- canonical_gt(
    focal$TN_Left_GT
)

focal$right_GT <- canonical_gt(
    focal$TN_Right_GT
)

focal$left_DP <- as.numeric(
    focal$TN_Left_DP
)

focal$left_GQ <- as.numeric(
    focal$TN_Left_GQ
)

focal$right_DP <- as.numeric(
    focal$TN_Right_DP
)

focal$right_GQ <- as.numeric(
    focal$TN_Right_GQ
)

focal$left_QC <- (
    !is.na(focal$left_GT)
    & is.finite(focal$left_DP)
    & is.finite(focal$left_GQ)
    & focal$left_DP >= 2
    & focal$left_GQ >= 20
)

focal$right_QC <- (
    !is.na(focal$right_GT)
    & is.finite(focal$right_DP)
    & is.finite(focal$right_GQ)
    & focal$right_DP >= 2
    & focal$right_GQ >= 20
)

focal$both_QC <- (
    focal$left_QC
    & focal$right_QC
)

focal$same_GT <- (
    focal$both_QC
    & focal$left_GT == focal$right_GT
)

focal$different_GT <- (
    focal$both_QC
    & focal$left_GT != focal$right_GT
)

focal$exact_expected_LR <- (
    focal$both_QC
    & focal$left_GT == "1/1"
    & focal$right_GT == "0/1"
)

focal$pair_class <- "One or both calls fail QC"

focal$pair_class[
    focal$same_GT
] <- "Same genotype"

focal$pair_class[
    focal$different_GT
] <- "Other left-right difference"

focal$pair_class[
    focal$exact_expected_LR
] <- "Exact P/P to P/A"


# Add comparative BA values to the focal-site table.
dennis_index <- match(
    focal$site_key,
    dennis$site_key
)

ray_index <- match(
    focal$site_key,
    ray$site_key
)

focal$Dennis_BA <- dennis$balanced_accuracy[
    dennis_index
]

focal$Dennis_score_eligible <- dennis$score_eligible_bool[
    dennis_index
]

focal$Ray_BA <- ray$balanced_accuracy[
    ray_index
]

focal$Ray_score_eligible <- ray$score_eligible_bool[
    ray_index
]

focal$Dennis_perfect_selected <- (
    focal$site_key %in% perfect_Dennis
)

focal$Ray_perfect_selected <- (
    focal$site_key %in% perfect_Ray
)


# Long representation for separate left and right genotype tracks.
focal_long <- rbind(
    data.frame(
        site_key = focal$site_key,
        position_Mb = focal$position_Mb,
        side = "TN_Left",
        y = 1,
        genotype = ifelse(
            focal$left_QC,
            focal$left_GT,
            "No call"
        ),
        pair_class = focal$pair_class,
        exact_expected_LR = focal$exact_expected_LR
    ),
    data.frame(
        site_key = focal$site_key,
        position_Mb = focal$position_Mb,
        side = "TN_Right",
        y = 0,
        genotype = ifelse(
            focal$right_QC,
            focal$right_GT,
            "No call"
        ),
        pair_class = focal$pair_class,
        exact_expected_LR = focal$exact_expected_LR
    )
)

focal_long$genotype <- factor(
    focal_long$genotype,
    levels = c(
        "0/0",
        "0/1",
        "1/1",
        "No call"
    )
)


shared_rows <- focal_long[
    focal_long$pair_class == "Same genotype",
    ,
    drop = FALSE
]

# Display only the exact expected Dennis/Ray focal state:
# TN_Left = P/P and TN_Right = P/A.
different_rows <- focal_long[
    focal_long$exact_expected_LR,
    ,
    drop = FALSE
]

failed_rows <- focal_long[
    focal_long$pair_class
    == "One or both calls fail QC",
    ,
    drop = FALSE
]

# The focal ribbons show only informative, QC-passing
# left-right genotype differences. Identical calls and sites
# where either pooled side fails QC are omitted from the plot.
shared_rows <- shared_rows[
    FALSE,
    ,
    drop = FALSE
]

failed_rows <- failed_rows[
    FALSE,
    ,
    drop = FALSE
]

exact_segments <- focal[
    focal$exact_expected_LR,
    ,
    drop = FALSE
]


plot_intervals <- data.frame(
    module = c(
        "Dennis",
        "Dennis",
        "Ray",
        "Ray"
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


ba_data <- scores[
    is.finite(
        scores$balanced_accuracy
    ),
    ,
    drop = FALSE
]

ba_data$panel <- factor(
    paste(
        ba_data$module,
        "comparative BA"
    ),
    levels = c(
        "Dennis comparative BA",
        "Ray comparative BA"
    )
)

ba_ineligible <- ba_data[
    !ba_data$score_eligible_bool,
    ,
    drop = FALSE
]

ba_eligible <- ba_data[
    ba_data$score_eligible_bool,
    ,
    drop = FALSE
]

ba_perfect <- ba_data[
    ba_data$perfect_selected,
    ,
    drop = FALSE
]


ba_intervals <- rbind(
    data.frame(
        panel = factor(
            "Dennis comparative BA",
            levels = levels(
                ba_data$panel
            )
        ),
        plot_intervals[
            plot_intervals$module == "Dennis",
            c(
                "xmin",
                "xmax"
            )
        ]
    ),
    data.frame(
        panel = factor(
            "Ray comparative BA",
            levels = levels(
                ba_data$panel
            )
        ),
        plot_intervals[
            plot_intervals$module == "Ray",
            c(
                "xmin",
                "xmax"
            )
        ]
    )
)


module_colours <- c(
    Dennis = "#E69F00",
    Ray = "#C51B7D"
)

genotype_colours <- c(
    "0/0" = "#3B6FB6",
    "0/1" = "#E69F00",
    "1/1" = "#C51B7D",
    "No call" = "#BDBDBD"
)

x_limits <- range(
    scores$position_Mb,
    finite = TRUE
)


# OPTIX_ZOOM_1060_1240_AND_GENE_ANNOTATION

# Presentation-only zoom and gene annotation.
# These values are in Mb because the plotting tables use pos_mb.
x_limits <- c(
    1.060,
    1.240
)

optix_start_mb <- 1088602 / 1e6
optix_end_mb <- 1097846 / 1e6
optix_midpoint_mb <- (
    optix_start_mb + optix_end_mb
) / 2

ba_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_mb,
        xmax = optix_end_mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.11
    ) +

    annotate(
        "segment",
        x = optix_start_mb,
        xend = optix_end_mb,
        y = 1.027,
        yend = 1.027,
        colour = "#2166AC",
        linewidth = 2.2,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = optix_midpoint_mb,
        y = 1.038,
        label = "optix",
        colour = "#2166AC",
        fontface = "italic",
        size = 3.5
    ) +

    geom_rect(
        data = ba_intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.11
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

    geom_point(
        data = ba_ineligible,
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
        data = ba_eligible,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = module
        ),
        alpha = 0.62,
        size = 0.48,
        stroke = 0
    ) +

    geom_point(
        data = ba_perfect,
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        shape = 21,
        fill = "white",
        colour = "#111111",
        size = 1.7,
        stroke = 0.55
    ) +

    facet_grid(
        rows = vars(
            panel
        )
    ) +

    scale_colour_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = NULL,
        limits = x_limits,
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
            1.045
        ),
        clip = "off"
    ) +

    labs(
        title = paste(
            "Comparative association and pooled focal genotypes",
            "across the broad Optix interval"
        ),
        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
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
            size = 9
        ),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.y = unit(
            0.35,
            "lines"
        ),
        plot.margin = margin(
            10,
            16,
            0,
            10
        )
    )


focal_intervals <- plot_intervals[
    ,
    c(
        "xmin",
        "xmax"
    )
]


focal_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_mb,
        xmax = optix_end_mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.11
    ) +

    geom_rect(
        data = focal_intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.11
    ) +

    geom_point(
        data = shared_rows,
        aes(
            x = position_Mb,
            y = y,
            colour = genotype
        ),
        shape = 124,
        size = 1.8,
        alpha = 0.16
    ) +

    geom_point(
        data = failed_rows,
        aes(
            x = position_Mb,
            y = y
        ),
        shape = 124,
        size = 1.5,
        colour = "#BDBDBD",
        alpha = 0.20
    ) +

    geom_segment(
        data = exact_segments,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = 0,
            yend = 1
        ),
        colour = "#009E73",
        alpha = 0.40,
        linewidth = 0.32
    ) +

    geom_point(
        data = different_rows,
        aes(
            x = position_Mb,
            y = y,
            colour = genotype
        ),
        shape = 124,
        size = 2.2,
        alpha = 0.88
    ) +

    geom_point(
        data = different_rows[
            different_rows$exact_expected_LR,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            y = y,
            fill = genotype
        ),
        shape = 21,
        colour = "#006B50",
        size = 1.25,
        stroke = 0.35
    ) +

    scale_colour_manual(
        values = genotype_colours,
        drop = FALSE,
        name = "Pooled genotype",
        labels = c(
            "0/0" = "A/A (reference homozygote)",
            "0/1" = "P/A (heterozygote)",
            "1/1" = "P/P (pattern homozygote)",
            "No call" = "No QC-passing call"
        )
    ) +

    scale_fill_manual(
        values = genotype_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = "Position across the broad Optix interval (Mb)",
        limits = x_limits,
        expand = expansion(
            mult = c(
                0.005,
                0.005
            )
        )
    ) +

    scale_y_continuous(
        name = NULL,
        breaks = c(
            0,
            1
        ),
        labels = c(
            "TN_Right",
            "TN_Left"
        ),
        limits = c(
            -0.35,
            1.35
        )
    ) +

    labs(
        title = "Pooled genotypes at exact expected left-right sites",
        subtitle = paste(
            "Only positions with TN_Left=P/P and TN_Right=P/A are shown.",
            "SNP selection and balanced accuracy use comparative samples only."
        ),
        caption = paste(
            "All other genotype combinations, identical calls and",
            "QC-failed positions are omitted from this display.",
            "Genotype color is shown separately for each pooled side."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 11
        ),
        plot.subtitle = element_text(
            size = 8.5,
            colour = "#444444"
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
            colour = "#E0E0E0",
            linewidth = 0.35
        ),
        axis.text.y = element_text(
            face = "bold"
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            2,
            16,
            10,
            10
        )
    )


ba_grob <- ggplotGrob(
    ba_plot
)

focal_grob <- ggplotGrob(
    focal_plot
)

common_widths <- unit.pmax(
    ba_grob$widths,
    focal_grob$widths
)

ba_grob$widths <- common_widths
focal_grob$widths <- common_widths


draw_combined <- function() {
    grid.newpage()

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 2,
                ncol = 1,
                heights = unit(
                    c(
                        2.25,
                        1
                    ),
                    "null"
                )
            )
        )
    )

    pushViewport(
        viewport(
            layout.pos.row = 1,
            layout.pos.col = 1
        )
    )

    grid.draw(
        ba_grob
    )

    popViewport()

    pushViewport(
        viewport(
            layout.pos.row = 2,
            layout.pos.col = 1
        )
    )

    grid.draw(
        focal_grob
    )

    popViewport(2)
}


pdf_file <- file.path(
    figure_dir,
    "broad_D_R_DP2_GQ20_BA_manhattan.pdf"
)

png_file <- file.path(
    figure_dir,
    "broad_D_R_DP2_GQ20_BA_manhattan.png"
)

cairo_pdf(
    filename = pdf_file,
    width = 11.5,
    height = 9.4
)

draw_combined()

dev.off()


png(
    filename = png_file,
    width = 11.5,
    height = 9.4,
    units = "in",
    res = 450
)

draw_combined()

dev.off()


pair_summary <- data.frame(
    focal_pair_state = c(
        "both_A_A",
        "both_P_A",
        "both_P_P",
        "exact_left_P_P_right_P_A",
        "other_left_right_difference",
        "one_or_both_fail_QC"
    ),
    SNPs = c(
        sum(
            focal$both_QC
            & focal$left_GT == "0/0"
            & focal$right_GT == "0/0"
        ),
        sum(
            focal$both_QC
            & focal$left_GT == "0/1"
            & focal$right_GT == "0/1"
        ),
        sum(
            focal$both_QC
            & focal$left_GT == "1/1"
            & focal$right_GT == "1/1"
        ),
        sum(
            focal$exact_expected_LR
        ),
        sum(
            focal$different_GT
            & !focal$exact_expected_LR
        ),
        sum(
            !focal$both_QC
        )
    )
)


overlap_summary <- data.frame(
    module = c(
        "Dennis",
        "Ray"
    ),
    focal_exact_sites = c(
        sum(
            focal$exact_expected_LR
        ),
        sum(
            focal$exact_expected_LR
        )
    ),
    focal_exact_and_BA80 = c(
        sum(
            focal$exact_expected_LR
            & focal$Dennis_score_eligible
            & focal$Dennis_BA >= 0.80,
            na.rm = TRUE
        ),
        sum(
            focal$exact_expected_LR
            & focal$Ray_score_eligible
            & focal$Ray_BA >= 0.80,
            na.rm = TRUE
        )
    ),
    focal_exact_and_BA90 = c(
        sum(
            focal$exact_expected_LR
            & focal$Dennis_score_eligible
            & focal$Dennis_BA >= 0.90,
            na.rm = TRUE
        ),
        sum(
            focal$exact_expected_LR
            & focal$Ray_score_eligible
            & focal$Ray_BA >= 0.90,
            na.rm = TRUE
        )
    ),
    focal_exact_and_BA1 = c(
        sum(
            focal$exact_expected_LR
            & focal$Dennis_score_eligible
            & abs(
                focal$Dennis_BA - 1
            ) < 1e-12,
            na.rm = TRUE
        ),
        sum(
            focal$exact_expected_LR
            & focal$Ray_score_eligible
            & abs(
                focal$Ray_BA - 1
            ) < 1e-12,
            na.rm = TRUE
        )
    ),
    focal_exact_and_selected_perfect = c(
        sum(
            focal$exact_expected_LR
            & focal$Dennis_perfect_selected
        ),
        sum(
            focal$exact_expected_LR
            & focal$Ray_perfect_selected
        )
    )
)


exact_table <- focal[
    focal$exact_expected_LR,
    c(
        "chrom",
        "pos",
        "ref",
        "alt",
        "TN_Left_GT",
        "TN_Left_AD",
        "TN_Left_DP",
        "TN_Left_GQ",
        "TN_Right_GT",
        "TN_Right_AD",
        "TN_Right_DP",
        "TN_Right_GQ",
        "Dennis_BA",
        "Dennis_score_eligible",
        "Dennis_perfect_selected",
        "Ray_BA",
        "Ray_score_eligible",
        "Ray_perfect_selected"
    ),
    drop = FALSE
]

pair_summary_file <- file.path(
    checkpoint_dir,
    "broad_D_R_focal_genotype_pair_summary.tsv"
)

overlap_summary_file <- file.path(
    checkpoint_dir,
    "broad_D_R_focal_exact_comparative_BA_overlap.tsv"
)

exact_table_file <- file.path(
    checkpoint_dir,
    "broad_D_R_focal_exact_sites_with_comparative_BA.tsv"
)

write.table(
    pair_summary,
    file = pair_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    overlap_summary,
    file = overlap_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    exact_table,
    file = exact_table_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "BROAD D/R BA WITH SEPARATE FOCAL GENOTYPE TRACKS\n\n"
)

cat(
    "FOCAL GENOTYPE PAIRS\n"
)

print(
    pair_summary,
    row.names = FALSE
)

cat(
    "\nFOCAL-EXACT / COMPARATIVE-BA OVERLAP\n"
)

print(
    overlap_summary,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nPair summary:\n",
    pair_summary_file,
    "\n\nOverlap summary:\n",
    overlap_summary_file,
    "\n\nExact-site table:\n",
    exact_table_file,
    "\n",
    sep = ""
)
