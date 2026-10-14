#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "B_western_plus_heurippa_reanalysis",
    "checkpoints"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "B_western_plus_heurippa_reanalysis",
    "broad_B_western_BA_manhattan"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

score_file <- file.path(checkpoint_dir, "B_western_BA_and_adjusted_LR.tsv")

if (!file.exists(score_file)) {
    stop(
        paste(
            "Missing Band score table:",
            score_file
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


band_biological_gt <- function(x) {
    result <- rep(
        "No call",
        length(x)
    )

    # Band reference allele is pattern-associated.
    result[
        x == "0/0"
    ] <- "P/P"

    result[
        x == "0/1"
    ] <- "P/A"

    result[
        x == "1/1"
    ] <- "A/A"

    result
}


d <- read.delim(
    score_file,
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
            "Band score table is missing columns:",
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

d$pos <- as.numeric(
    d$pos
)

d$position_Mb <- (
    d$pos / 1000000
)

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

d <- d[
    !duplicated(
        d$site_key
    ),
    ,
    drop = FALSE
]


d$left_GT <- canonical_gt(
    d$TN_Left_GT
)

d$right_GT <- canonical_gt(
    d$TN_Right_GT
)

d$left_DP <- as.numeric(
    d$TN_Left_DP
)

d$left_GQ <- as.numeric(
    d$TN_Left_GQ
)

d$right_DP <- as.numeric(
    d$TN_Right_DP
)

d$right_GQ <- as.numeric(
    d$TN_Right_GQ
)


if (!("mixture_pair_posterior" %in% names(d))) {
    stop("Missing adjusted-call posterior")
}

d$adjusted_confident <- (
    is.finite(d$mixture_pair_posterior)
    & d$mixture_pair_posterior >= 0.95
)

d$left_QC <- (
    d$adjusted_confident
    & !is.na(d$left_GT)
)

d$right_QC <- (
    d$adjusted_confident
    & !is.na(d$right_GT)
)

d$both_QC <- (
    d$left_QC
    & d$right_QC
)

d$same_GT <- (
    d$both_QC
    & d$left_GT == d$right_GT
)

d$different_GT <- (
    d$both_QC
    & d$left_GT != d$right_GT
)

# Band expectation:
# left biological A/A = raw 1/1
# right biological P/A = raw 0/1
d$exact_expected_LR <- (
    d$both_QC
    & d$left_GT == "1/1"
    & d$right_GT == "0/1"
)

d$left_biological_GT <- band_biological_gt(
    d$left_GT
)

d$right_biological_GT <- band_biological_gt(
    d$right_GT
)

d$perfect_BA <- (
    d$score_eligible_bool
    & is.finite(d$balanced_accuracy)
    & abs(
        d$balanced_accuracy - 1
    ) < 1e-12
)


ba_data <- d[
    is.finite(
        d$balanced_accuracy
    ),
    ,
    drop = FALSE
]

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
    ba_data$perfect_BA,
    ,
    drop = FALSE
]


# Display only the exact expected Band focal state:
# TN_Left biological A/A and TN_Right biological P/A.
difference_sites <- d[
    d$exact_expected_LR,
    ,
    drop = FALSE
]

focal_long <- rbind(
    data.frame(
        site_key = difference_sites$site_key,
        position_Mb = difference_sites$position_Mb,
        side = rep("TN_Left", nrow(difference_sites)),
        y = rep(1, nrow(difference_sites)),
        biological_GT = difference_sites$left_biological_GT,
        exact_expected_LR = difference_sites$exact_expected_LR
    ),
    data.frame(
        site_key = difference_sites$site_key,
        position_Mb = difference_sites$position_Mb,
        side = rep("TN_Right", nrow(difference_sites)),
        y = rep(0, nrow(difference_sites)),
        biological_GT = difference_sites$right_biological_GT,
        exact_expected_LR = difference_sites$exact_expected_LR
    )
)

focal_long$biological_GT <- factor(
    focal_long$biological_GT,
    levels = c(
        "A/A",
        "P/A",
        "P/P",
        "No call"
    )
)

exact_segments <- d[
    d$exact_expected_LR,
    ,
    drop = FALSE
]


band_interval <- data.frame(
    xmin = 1112651 / 1000000,
    xmax = 1114314 / 1000000
)

# Requested display interval:
# Hmr18:1,060,000-1,240,000
x_limits <- c(
    1060000 / 1000000,
    1240000 / 1000000
)

# optix annotation in the same Hmr18 coordinate system.
optix_start_Mb <- 1088602 / 1000000
optix_end_Mb <- 1097846 / 1000000
optix_midpoint_Mb <- (
    optix_start_Mb
    + optix_end_Mb
) / 2

genotype_colours <- c(
    "A/A" = "#3B6FB6",
    "P/A" = "#E69F00",
    "P/P" = "#C51B7D",
    "No call" = "#BDBDBD"
)


ba_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_Mb,
        xmax = optix_end_Mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    annotate(
        "segment",
        x = optix_start_Mb,
        xend = optix_end_Mb,
        y = 1.017,
        yend = 1.017,
        colour = "#276B8A",
        linewidth = 1.25,
        lineend = "butt"
    ) +

    annotate(
        "segment",
        x = optix_start_Mb,
        xend = optix_start_Mb,
        y = 1.012,
        yend = 1.022,
        colour = "#276B8A",
        linewidth = 0.55
    ) +

    annotate(
        "segment",
        x = optix_end_Mb,
        xend = optix_end_Mb,
        y = 1.012,
        yend = 1.022,
        colour = "#276B8A",
        linewidth = 0.55
    ) +

    annotate(
        "text",
        x = optix_midpoint_Mb,
        y = 1.027,
        label = "optix",
        colour = "#194B63",
        fontface = "italic",
        size = 3.4,
        vjust = 0
    ) +

    geom_rect(
        data = band_interval,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.13
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
            y = balanced_accuracy
        ),
        colour = "#D55E00",
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
        name = "Band balanced accuracy",
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
            "Band association in western samples plus heurippa",
            "across the broad Optix interval"
        ),
        subtitle = paste(
            "Hmr18:1.060-1.240 Mb.",
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
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(
            10,
            16,
            0,
            10
        )
    )


focal_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_Mb,
        xmax = optix_end_Mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    geom_rect(
        data = band_interval,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.13
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
        alpha = 0.42,
        linewidth = 0.34
    ) +

    geom_point(
        data = focal_long,
        aes(
            x = position_Mb,
            y = y,
            colour = biological_GT
        ),
        shape = 124,
        size = 2.3,
        alpha = 0.88
    ) +

    geom_point(
        data = focal_long[
            focal_long$exact_expected_LR,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            y = y,
            fill = biological_GT
        ),
        shape = 21,
        colour = "#006B50",
        size = 1.3,
        stroke = 0.38
    ) +

    scale_colour_manual(
        values = genotype_colours,
        drop = FALSE,
        name = "Band biological genotype"
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
        title = "Mixture-adjusted genotypes at exact expected left-right sites",
        subtitle = paste(
            "Only positions with TN_Left=A/A and TN_Right=P/A are shown.",
            "Pair posterior >=0.95; BA is calculated independently of focal calls."
        ),
        caption = paste(
            "For Band, reference homozygote 0/0 is P/P;",
            "heterozygote 0/1 is P/A; ALT homozygote 1/1 is A/A.",
            "Other pairs and uncertain adjusted calls are omitted; no original-call fallback."
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
                        1.8,
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
    "broad_B_DP2_GQ20_BA_and_focal_genotypes.pdf"
)

png_file <- file.path(
    figure_dir,
    "broad_B_DP2_GQ20_BA_and_focal_genotypes.png"
)

cairo_pdf(
    filename = pdf_file,
    width = 11.5,
    height = 7.4
)

draw_combined()

dev.off()

png(
    filename = png_file,
    width = 11.5,
    height = 7.4,
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
        "exact_left_A_A_right_P_A",
        "other_left_right_difference",
        "one_or_both_fail_QC"
    ),
    SNPs = c(
        sum(
            d$both_QC
            & d$left_biological_GT == "A/A"
            & d$right_biological_GT == "A/A"
        ),
        sum(
            d$both_QC
            & d$left_biological_GT == "P/A"
            & d$right_biological_GT == "P/A"
        ),
        sum(
            d$both_QC
            & d$left_biological_GT == "P/P"
            & d$right_biological_GT == "P/P"
        ),
        sum(
            d$exact_expected_LR
        ),
        sum(
            d$different_GT
            & !d$exact_expected_LR
        ),
        sum(
            !d$both_QC
        )
    )
)


overlap_summary <- data.frame(
    statistic = c(
        "physical_SNPs",
        "score_eligible",
        "BA_ge_0_80",
        "BA_ge_0_90",
        "BA_equal_1",
        "focal_exact_left_A_A_right_P_A",
        "focal_exact_and_BA80",
        "focal_exact_and_BA90",
        "focal_exact_and_BA1"
    ),
    SNPs = c(
        nrow(d),
        sum(
            d$score_eligible_bool
        ),
        sum(
            d$score_eligible_bool
            & d$balanced_accuracy >= 0.80,
            na.rm = TRUE
        ),
        sum(
            d$score_eligible_bool
            & d$balanced_accuracy >= 0.90,
            na.rm = TRUE
        ),
        sum(
            d$perfect_BA
        ),
        sum(
            d$exact_expected_LR
        ),
        sum(
            d$exact_expected_LR
            & d$score_eligible_bool
            & d$balanced_accuracy >= 0.80,
            na.rm = TRUE
        ),
        sum(
            d$exact_expected_LR
            & d$score_eligible_bool
            & d$balanced_accuracy >= 0.90,
            na.rm = TRUE
        ),
        sum(
            d$exact_expected_LR
            & d$perfect_BA
        )
    )
)


exact_table <- d[
    d$exact_expected_LR,
    c(
        "chrom",
        "pos",
        "ref",
        "alt",
        "balanced_accuracy",
        "score_eligible_bool",
        "TN_Left_GT",
        "left_biological_GT",
        "TN_Left_AD",
        "TN_Left_DP",
        "TN_Left_GQ",
        "TN_Right_GT",
        "right_biological_GT",
        "TN_Right_AD",
        "TN_Right_DP",
        "TN_Right_GQ"
    ),
    drop = FALSE
]


pair_summary_file <- file.path(
    checkpoint_dir,
    "broad_B_focal_genotype_pair_summary.tsv"
)

overlap_summary_file <- file.path(
    checkpoint_dir,
    "broad_B_focal_exact_comparative_BA_overlap.tsv"
)

exact_table_file <- file.path(
    checkpoint_dir,
    "broad_B_focal_exact_sites_with_comparative_BA.tsv"
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
    "BAND BA WITH SEPARATE FOCAL GENOTYPE TRACKS\n\n"
)

cat(
    "FOCAL GENOTYPE PAIRS\n"
)

print(
    pair_summary,
    row.names = FALSE
)

cat(
    "\nBAND / FOCAL-EXACT OVERLAP\n"
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
