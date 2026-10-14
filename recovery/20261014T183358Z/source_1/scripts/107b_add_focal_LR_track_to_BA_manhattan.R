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
        !(x %in% c(
            "0/0",
            "0/1",
            "1/0",
            "1/1"
        ))
    ] <- NA_character_

    x[
        x == "1/0"
    ] <- "0/1"

    x
}


gt_dosage <- function(x) {
    result <- rep(
        NA_real_,
        length(x)
    )

    result[
        x == "0/0"
    ] <- 0

    result[
        x == "0/1"
    ] <- 0.5

    result[
        x == "1/1"
    ] <- 1

    result
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

perfect_keys <- unique(
    paste(
        perfect$module,
        perfect$site_key,
        sep = "|"
    )
)

scores$perfect_selected <- (
    paste(
        scores$module,
        scores$site_key,
        sep = "|"
    )
    %in% perfect_keys
)


# Use Dennis rows as one nonduplicated copy of the physical
# SNPs. Dennis and Ray tables carry the same focal calls.
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

focal$both_focal_QC <- (
    !is.na(focal$left_GT)
    & !is.na(focal$right_GT)
    & is.finite(focal$left_DP)
    & is.finite(focal$left_GQ)
    & is.finite(focal$right_DP)
    & is.finite(focal$right_GQ)
    & focal$left_DP >= 2
    & focal$left_GQ >= 20
    & focal$right_DP >= 2
    & focal$right_GQ >= 20
)

focal$left_pattern_dosage <- gt_dosage(
    focal$left_GT
)

focal$right_pattern_dosage <- gt_dosage(
    focal$right_GT
)

focal$left_minus_right_dosage <- (
    focal$left_pattern_dosage
    - focal$right_pattern_dosage
)

focal$exact_expected_LR <- (
    focal$both_focal_QC
    & focal$left_GT == "1/1"
    & focal$right_GT == "0/1"
)

focal_plot <- focal[
    focal$both_focal_QC,
    ,
    drop = FALSE
]


ba_data <- scores[
    is.finite(
        scores$balanced_accuracy
    ),
    ,
    drop = FALSE
]

ba_data$panel <- ifelse(
    ba_data$module == "Dennis",
    "Dennis comparative BA",
    "Ray comparative BA"
)

ba_data$plot_value <- (
    ba_data$balanced_accuracy
)

ba_data$point_class <- ifelse(
    ba_data$score_eligible_bool,
    ba_data$module,
    "Not score-eligible"
)


focal_plot$panel <- (
    "Pooled focal L-R dosage"
)

focal_plot$plot_value <- (
    focal_plot$left_minus_right_dosage
)

focal_plot$point_class <- ifelse(
    focal_plot$exact_expected_LR,
    "Exact P/P to P/A",
    "Other focal genotype"
)


plot_data <- rbind(
    ba_data[
        ,
        c(
            "panel",
            "position_Mb",
            "plot_value",
            "point_class",
            "perfect_selected"
        )
    ],
    data.frame(
        panel = focal_plot$panel,
        position_Mb = focal_plot$position_Mb,
        plot_value = focal_plot$plot_value,
        point_class = focal_plot$point_class,
        perfect_selected = FALSE
    )
)

plot_data$panel <- factor(
    plot_data$panel,
    levels = c(
        "Dennis comparative BA",
        "Ray comparative BA",
        "Pooled focal L-R dosage"
    )
)


plot_intervals <- data.frame(
    panel = factor(
        c(
            "Dennis comparative BA",
            "Dennis comparative BA",
            "Ray comparative BA",
            "Ray comparative BA",
            "Pooled focal L-R dosage",
            "Pooled focal L-R dosage",
            "Pooled focal L-R dosage",
            "Pooled focal L-R dosage"
        ),
        levels = levels(
            plot_data$panel
        )
    ),
    interval = c(
        "Dennis1",
        "Dennis2",
        "Ray1",
        "Ray2",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    ),
    start = c(
        1215210,
        1218557,
        1156746,
        1174547,
        1156746,
        1174547,
        1215210,
        1218557
    ),
    end = c(
        1216124,
        1221494,
        1160220,
        1187833,
        1160220,
        1187833,
        1216124,
        1221494
    )
)

plot_intervals$xmin <- (
    plot_intervals$start / 1000000
)

plot_intervals$xmax <- (
    plot_intervals$end / 1000000
)


ba_thresholds <- data.frame(
    panel = factor(
        rep(
            c(
                "Dennis comparative BA",
                "Ray comparative BA"
            ),
            each = 3
        ),
        levels = levels(
            plot_data$panel
        )
    ),
    intercept = rep(
        c(
            0.50,
            0.80,
            0.90
        ),
        times = 2
    ),
    line_type = rep(
        c(
            "solid",
            "dashed",
            "dotted"
        ),
        times = 2
    )
)


focal_reference_lines <- data.frame(
    panel = factor(
        rep(
            "Pooled focal L-R dosage",
            2
        ),
        levels = levels(
            plot_data$panel
        )
    ),
    intercept = c(
        0,
        0.5
    ),
    line_type = c(
        "solid",
        "dashed"
    )
)


background_points <- plot_data[
    plot_data$point_class
    %in% c(
        "Not score-eligible",
        "Other focal genotype"
    ),
    ,
    drop = FALSE
]

comparative_points <- plot_data[
    plot_data$point_class
    %in% c(
        "Dennis",
        "Ray"
    ),
    ,
    drop = FALSE
]

perfect_points <- plot_data[
    plot_data$perfect_selected,
    ,
    drop = FALSE
]

exact_focal_points <- plot_data[
    plot_data$point_class
    == "Exact P/P to P/A",
    ,
    drop = FALSE
]


point_colours <- c(
    Dennis = "#E69F00",
    Ray = "#C51B7D"
)


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
        alpha = 0.11
    ) +

    geom_hline(
        data = ba_thresholds,
        aes(
            yintercept = intercept,
            linetype = line_type
        ),
        colour = "#555555",
        linewidth = 0.38,
        show.legend = FALSE
    ) +

    geom_hline(
        data = focal_reference_lines,
        aes(
            yintercept = intercept,
            linetype = line_type
        ),
        colour = "#3A7563",
        linewidth = 0.45,
        show.legend = FALSE
    ) +

    geom_point(
        data = background_points,
        aes(
            x = position_Mb,
            y = plot_value
        ),
        colour = "#BDBDBD",
        alpha = 0.30,
        size = 0.38,
        stroke = 0
    ) +

    geom_point(
        data = comparative_points,
        aes(
            x = position_Mb,
            y = plot_value,
            colour = point_class
        ),
        alpha = 0.62,
        size = 0.48,
        stroke = 0
    ) +

    geom_point(
        data = perfect_points,
        aes(
            x = position_Mb,
            y = plot_value
        ),
        shape = 21,
        fill = "white",
        colour = "#111111",
        size = 1.7,
        stroke = 0.55
    ) +

    geom_point(
        data = exact_focal_points,
        aes(
            x = position_Mb,
            y = plot_value
        ),
        shape = 21,
        fill = "#009E73",
        colour = "#123D32",
        size = 1.65,
        stroke = 0.5
    ) +

    facet_grid(
        rows = vars(
            panel
        ),
        scales = "free_y",
        space = "fixed"
    ) +

    scale_colour_manual(
        values = point_colours,
        guide = "none"
    ) +

    scale_linetype_identity() +

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
        name = "Balanced accuracy or left-minus-right pattern dosage"
    ) +

    labs(
        title = paste(
            "Comparative association and pooled focal genotypes",
            "across the broad Optix interval"
        ),
        subtitle = paste(
            "Comparative panels show final WGS+CS DP2/GQ20 balanced accuracy.",
            "The focal panel is an independent display of TN_Left and TN_Right."
        ),
        caption = paste(
            "Black rings: perfect comparative SNPs.",
            "Green points: exact focal TN_Left=P/P and TN_Right=P/A.",
            "In the focal panel, the dashed line at 0.5 is the expected",
            "P/P minus P/A dosage difference; zero indicates equal dosage.",
            "Gray focal points are other QC-passing genotype combinations."
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
            size = 9
        ),
        panel.spacing.y = grid::unit(
            0.4,
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
    height = 9.4,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = png_file,
    plot = plot,
    width = 11.5,
    height = 9.4,
    units = "in",
    dpi = 450,
    type = "cairo",
    limitsize = FALSE
)


focal_exact <- focal[
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
        "left_pattern_dosage",
        "right_pattern_dosage",
        "left_minus_right_dosage"
    ),
    drop = FALSE
]

focal_exact_file <- file.path(
    checkpoint_dir,
    "broad_D_R_focal_exact_expected_LR_sites.tsv"
)

write.table(
    focal_exact,
    file = focal_exact_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


summary <- data.frame(
    statistic = c(
        "physical_SNPs",
        "both_focal_DP2_GQ20",
        "exact_P_P_to_P_A",
        "left_dosage_greater_than_right",
        "equal_dosage",
        "left_dosage_less_than_right"
    ),
    value = c(
        nrow(focal),
        sum(focal$both_focal_QC),
        sum(focal$exact_expected_LR),
        sum(
            focal$both_focal_QC
            & focal$left_minus_right_dosage > 0
        ),
        sum(
            focal$both_focal_QC
            & focal$left_minus_right_dosage == 0
        ),
        sum(
            focal$both_focal_QC
            & focal$left_minus_right_dosage < 0
        )
    )
)

summary_file <- file.path(
    checkpoint_dir,
    "broad_D_R_focal_LR_manhattan.summary.tsv"
)

write.table(
    summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "BROAD D/R BA MANHATTAN WITH FOCAL L/R TRACK\n\n"
)

print(
    summary,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nExact focal sites:\n",
    focal_exact_file,
    "\n\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
