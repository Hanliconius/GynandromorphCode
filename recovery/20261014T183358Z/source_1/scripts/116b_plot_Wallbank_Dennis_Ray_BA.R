#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

analysis_dir <- file.path(
    project_dir,
    "analysis",
    "original_502",
    "Wallbank_Dennis_Ray_BA"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "Wallbank_Dennis_Ray_BA"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

score_file <- file.path(
    analysis_dir,
    "Wallbank_Dennis_Ray_dominant_BA.tsv"
)

focal_file <- file.path(
    analysis_dir,
    "Wallbank_Dennis_Ray_focal_calls.tsv"
)

exact_file <- file.path(
    analysis_dir,
    "Wallbank_Dennis_Ray_exact_focal_sites.tsv"
)

exact_summary_file <- file.path(
    analysis_dir,
    "Wallbank_Dennis_Ray_exact_focal_sites.summary.tsv"
)

pdf_file <- file.path(
    figure_dir,
    "Wallbank_Dennis_Ray_dominant_BA_and_focal_genotypes.pdf"
)

png_file <- file.path(
    figure_dir,
    "Wallbank_Dennis_Ray_dominant_BA_and_focal_genotypes.png"
)

if (!file.exists(score_file)) {
    stop(
        paste(
            "Missing score table:",
            score_file
        )
    )
}

if (!file.exists(focal_file)) {
    stop(
        paste(
            "Missing focal table:",
            focal_file
        )
    )
}


as_bool <- function(x) {
    tolower(trimws(as.character(x))) %in% c(
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
        as.character(x)
    )

    x[
        x %in% c(
            "",
            ".",
            "./.",
            ".|.",
            "NA"
        )
    ] <- NA_character_

    x[
        x == "1/0"
    ] <- "0/1"

    x
}


scores <- read.delim(
    score_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

focal <- read.delim(
    focal_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_score_columns <- c(
    "module",
    "chrom",
    "pos",
    "ref",
    "alt",
    "balanced_accuracy",
    "minimum_group_call_rate",
    "selected_orientation",
    "associated_allele",
    "eligible_min2",
    "callable_80pct",
    "complete_callability"
)

missing_score_columns <- setdiff(
    required_score_columns,
    names(scores)
)

if (length(missing_score_columns) > 0) {
    stop(
        paste(
            "Score table lacks columns:",
            paste(
                missing_score_columns,
                collapse = ", "
            )
        )
    )
}


required_focal_columns <- c(
    "chrom",
    "pos",
    "ref",
    "alt",
    "TN_Left_GT",
    "TN_Left_DP",
    "TN_Left_GQ",
    "TN_Left_DP2_GQ20",
    "TN_Right_GT",
    "TN_Right_DP",
    "TN_Right_GQ",
    "TN_Right_DP2_GQ20"
)

missing_focal_columns <- setdiff(
    required_focal_columns,
    names(focal)
)

if (length(missing_focal_columns) > 0) {
    stop(
        paste(
            "Focal table lacks columns:",
            paste(
                missing_focal_columns,
                collapse = ", "
            )
        )
    )
}


scores$pos <- as.numeric(
    scores$pos
)

scores$position_Mb <- (
    scores$pos / 1000000
)

scores$balanced_accuracy <- as.numeric(
    scores$balanced_accuracy
)

scores$minimum_group_call_rate <- as.numeric(
    scores$minimum_group_call_rate
)

scores$eligible_min2_bool <- as_bool(
    scores$eligible_min2
)

scores$callable_80pct_bool <- as_bool(
    scores$callable_80pct
)

scores$complete_callability_bool <- as_bool(
    scores$complete_callability
)

scores$module <- factor(
    scores$module,
    levels = c(
        "Ray",
        "Dennis"
    )
)


focal$pos <- as.numeric(
    focal$pos
)

focal$TN_Left_GT <- canonical_gt(
    focal$TN_Left_GT
)

focal$TN_Right_GT <- canonical_gt(
    focal$TN_Right_GT
)

focal$TN_Left_DP2_GQ20_bool <- as_bool(
    focal$TN_Left_DP2_GQ20
)

focal$TN_Right_DP2_GQ20_bool <- as_bool(
    focal$TN_Right_DP2_GQ20
)


scores$site_key <- paste(
    scores$chrom,
    scores$pos,
    scores$ref,
    scores$alt,
    sep = ":"
)

focal$site_key <- paste(
    focal$chrom,
    focal$pos,
    focal$ref,
    focal$alt,
    sep = ":"
)


joined <- merge(
    scores,
    focal,
    by = c(
        "site_key",
        "chrom",
        "pos",
        "ref",
        "alt"
    ),
    all.x = TRUE,
    sort = FALSE
)

joined$module <- factor(
    joined$module,
    levels = c(
        "Ray",
        "Dennis"
    )
)

joined$position_Mb <- (
    joined$pos / 1000000
)


# ------------------------------------------------------------
# Zygosity-aware focal interpretation.
#
# If ALT is the comparative pattern-associated allele:
# P/P = 1/1 and P/A = 0/1.
#
# If REF is pattern-associated:
# P/P = 0/0 and P/A = 0/1.
# ------------------------------------------------------------

joined$expected_left_raw_GT <- ifelse(
    joined$associated_allele == "1",
    "1/1",
    "0/0"
)

joined$expected_right_raw_GT <- "0/1"

joined$focal_both_QC <- (
    joined$TN_Left_DP2_GQ20_bool
    & joined$TN_Right_DP2_GQ20_bool
)

joined$exact_expected_focal <- (
    joined$focal_both_QC
    & joined$TN_Left_GT
        == joined$expected_left_raw_GT
    & joined$TN_Right_GT
        == joined$expected_right_raw_GT
)


# ------------------------------------------------------------
# Comparative BA plotting classes.
# ------------------------------------------------------------

ba_data <- joined[
    joined$eligible_min2_bool
    & is.finite(
        joined$balanced_accuracy
    ),
    ,
    drop = FALSE
]

ba_data$high_BA <- (
    ba_data$balanced_accuracy >= 0.80
)

ba_data$perfect_BA <- (
    abs(
        ba_data$balanced_accuracy - 1
    ) < 1e-12
)

ba_data$perfect_class <- "Not perfect"

ba_data$perfect_class[
    ba_data$perfect_BA
] <- "Perfect among called; <80% callability"

ba_data$perfect_class[
    ba_data$perfect_BA
    & ba_data$callable_80pct_bool
] <- "Perfect; at least 80% callable"

ba_data$perfect_class[
    ba_data$perfect_BA
    & ba_data$complete_callability_bool
] <- "Perfect; completely called"


ba_background <- ba_data[
    !ba_data$high_BA,
    ,
    drop = FALSE
]

ba_high <- ba_data[
    ba_data$high_BA,
    ,
    drop = FALSE
]

perfect_sparse <- ba_data[
    ba_data$perfect_class
        == "Perfect among called; <80% callability",
    ,
    drop = FALSE
]

perfect_80 <- ba_data[
    ba_data$perfect_class
        == "Perfect; at least 80% callable",
    ,
    drop = FALSE
]

perfect_complete <- ba_data[
    ba_data$perfect_class
        == "Perfect; completely called",
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# Independent focal track: show only positions with the exact
# expected left P/P and right P/A state.
# ------------------------------------------------------------

focal_exact <- joined[
    joined$eligible_min2_bool
    & joined$exact_expected_focal,
    ,
    drop = FALSE
]

focal_exact$left_biological_GT <- "P/P"
focal_exact$right_biological_GT <- "P/A"

focal_exact$comparative_class <- cut(
    focal_exact$balanced_accuracy,
    breaks = c(
        -Inf,
        0.80,
        0.90,
        1.0 - 1e-12,
        Inf
    ),
    labels = c(
        "BA < 0.80",
        "0.80 <= BA < 0.90",
        "0.90 <= BA < 1",
        "BA = 1"
    ),
    right = FALSE
)


focal_long <- rbind(
    data.frame(
        module = focal_exact$module,
        site_key = focal_exact$site_key,
        pos = focal_exact$pos,
        position_Mb = focal_exact$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "P/P",
        raw_GT = focal_exact$TN_Left_GT,
        balanced_accuracy = (
            focal_exact$balanced_accuracy
        ),
        minimum_group_call_rate = (
            focal_exact$minimum_group_call_rate
        ),
        comparative_class = (
            focal_exact$comparative_class
        ),
        stringsAsFactors = FALSE
    ),
    data.frame(
        module = focal_exact$module,
        site_key = focal_exact$site_key,
        pos = focal_exact$pos,
        position_Mb = focal_exact$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        raw_GT = focal_exact$TN_Right_GT,
        balanced_accuracy = (
            focal_exact$balanced_accuracy
        ),
        minimum_group_call_rate = (
            focal_exact$minimum_group_call_rate
        ),
        comparative_class = (
            focal_exact$comparative_class
        ),
        stringsAsFactors = FALSE
    )
)

focal_long$module <- factor(
    focal_long$module,
    levels = c(
        "Ray",
        "Dennis"
    )
)

focal_long$biological_GT <- factor(
    focal_long$biological_GT,
    levels = c(
        "A/A",
        "P/A",
        "P/P"
    )
)


write.table(
    focal_exact[
        ,
        c(
            "module",
            "chrom",
            "pos",
            "ref",
            "alt",
            "selected_orientation",
            "associated_allele",
            "balanced_accuracy",
            "minimum_group_call_rate",
            "callable_80pct_bool",
            "complete_callability_bool",
            "TN_Left_GT",
            "TN_Left_DP",
            "TN_Left_GQ",
            "TN_Right_GT",
            "TN_Right_DP",
            "TN_Right_GQ",
            "expected_left_raw_GT",
            "expected_right_raw_GT"
        )
    ],
    file = exact_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


exact_summary <- do.call(
    rbind,
    lapply(
        levels(joined$module),
        function(module_name) {
            module_rows <- focal_exact[
                focal_exact$module == module_name,
                ,
                drop = FALSE
            ]

            data.frame(
                module = module_name,
                focal_exact_sites = nrow(
                    module_rows
                ),
                focal_exact_BA80 = sum(
                    module_rows$balanced_accuracy
                        >= 0.80,
                    na.rm = TRUE
                ),
                focal_exact_BA90 = sum(
                    module_rows$balanced_accuracy
                        >= 0.90,
                    na.rm = TRUE
                ),
                focal_exact_BA1 = sum(
                    abs(
                        module_rows$balanced_accuracy
                            - 1
                    ) < 1e-12,
                    na.rm = TRUE
                ),
                focal_exact_BA1_80pct_callable = sum(
                    abs(
                        module_rows$balanced_accuracy
                            - 1
                    ) < 1e-12
                    & module_rows$callable_80pct_bool,
                    na.rm = TRUE
                ),
                focal_exact_BA1_complete = sum(
                    abs(
                        module_rows$balanced_accuracy
                            - 1
                    ) < 1e-12
                    & module_rows$complete_callability_bool,
                    na.rm = TRUE
                ),
                stringsAsFactors = FALSE
            )
        }
    )
)

write.table(
    exact_summary,
    file = exact_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Plot annotations and colours.
# ------------------------------------------------------------

x_limits <- c(
    1.060,
    1.240
)

optix_start <- 1088602 / 1000000
optix_end <- 1097846 / 1000000

module_colours <- c(
    "Dennis" = "#B51A7D",
    "Ray" = "#E88924"
)

genotype_colours <- c(
    "A/A" = "#B51A7D",
    "P/A" = "#E88924",
    "P/P" = "#3C9DC7"
)


gene_annotation <- data.frame(
    module = factor(
        c(
            "Ray",
            "Dennis"
        ),
        levels = c(
            "Ray",
            "Dennis"
        )
    ),
    xmin = optix_start,
    xmax = optix_end
)


intervals <- data.frame(
    module = factor(
        c(
            "Ray",
            "Ray",
            "Dennis",
            "Dennis"
        ),
        levels = c(
            "Ray",
            "Dennis"
        )
    ),
    interval = c(
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    ),
    xmin = c(
        1156746,
        1174547,
        1215210,
        1218557
    ) / 1000000,
    xmax = c(
        1160220,
        1187833,
        1216124,
        1221494
    ) / 1000000
)


ba_plot <- ggplot() +

    geom_rect(
        data = gene_annotation,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.10
    ) +

    geom_rect(
        data = intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = 0.475,
            ymax = 0.505
        ),
        inherit.aes = FALSE,
        fill = "#777777",
        alpha = 0.23
    ) +

    geom_text(
        data = intervals,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 0.515,
            label = interval
        ),
        inherit.aes = FALSE,
        size = 2.45,
        colour = "#555555"
    ) +

    geom_hline(
        yintercept = 0.50,
        colour = "#888888",
        linewidth = 0.35,
        linetype = "solid"
    ) +

    geom_hline(
        yintercept = 0.80,
        colour = "#999999",
        linewidth = 0.25,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 0.90,
        colour = "#999999",
        linewidth = 0.25,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 1.00,
        colour = "#555555",
        linewidth = 0.35,
        linetype = "solid"
    ) +

    geom_point(
        data = ba_background,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            alpha = minimum_group_call_rate
        ),
        colour = "#B8B8B8",
        size = 0.34,
        stroke = 0
    ) +

    geom_point(
        data = ba_high,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = module,
            alpha = minimum_group_call_rate
        ),
        size = 0.52,
        stroke = 0
    ) +

    geom_point(
        data = perfect_sparse,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = module
        ),
        shape = 25,
        fill = "white",
        size = 2.05,
        stroke = 0.55
    ) +

    geom_point(
        data = perfect_80,
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            fill = module
        ),
        shape = 25,
        colour = "#222222",
        size = 2.20,
        stroke = 0.45
    ) +

    geom_point(
        data = perfect_complete,
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        shape = 25,
        fill = "#111111",
        colour = "#111111",
        size = 2.15,
        stroke = 0.35
    ) +

    geom_text(
        data = gene_annotation,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 1.012,
            label = "optix"
        ),
        inherit.aes = FALSE,
        fontface = "italic",
        size = 2.8,
        colour = "#277DA1"
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_fill_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_alpha_continuous(
        name = "Minimum group call rate",
        range = c(
            0.08,
            0.78
        ),
        limits = c(
            0,
            1
        )
    ) +

    scale_x_continuous(
        name = NULL,
        limits = x_limits,
        breaks = seq(
            1.06,
            1.24,
            by = 0.02
        ),
        labels = function(x) {
            sprintf("%.2f", x)
        },
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_y_continuous(
        name = "Balanced accuracy",
        breaks = c(
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0
        ),
        limits = c(
            0.47,
            1.025
        )
    ) +

    labs(
        title = paste(
            "Dominant genotype-phenotype association",
            "across optix"
        ),
        subtitle = paste(
            "Dennis: 4 no-D versus 2 D-present specimens;",
            "Ray: 5 no-Ray versus 10 D0R1 contigua.",
            "Comparative calls require DP 4-300 and GQ >=30."
        ),
        caption = paste(
            "Downward triangles mark BA=1:",
            "open = <80% callability;",
            "coloured = at least 80% callable;",
            "black = completely called.",
            "TN_Left and TN_Right were excluded from BA calculation."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#ECECEC",
            linewidth = 0.25
        ),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(
            face = "bold",
            size = 10.5
        ),
        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = NA
        ),
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9.2,
            colour = "#444444"
        ),
        plot.caption = element_text(
            size = 7.7,
            colour = "#555555",
            hjust = 0
        ),
        legend.position = "top",
        legend.justification = "left",
        axis.title.y = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            8,
            12,
            3,
            8
        )
    )


focal_plot <- ggplot() +

    geom_rect(
        data = gene_annotation,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.10
    ) +

    geom_segment(
        data = focal_exact,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = 0,
            yend = 1
        ),
        colour = "#009E73",
        linewidth = 0.30,
        alpha = 0.48
    ) +

    geom_point(
        data = focal_long,
        aes(
            x = position_Mb,
            y = y,
            colour = biological_GT
        ),
        shape = 124,
        size = 2.25,
        alpha = 0.90
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        name = "Biological genotype",
        values = genotype_colours,
        drop = FALSE
    ) +

    scale_x_continuous(
        name = "Position on Hmr18 (Mb)",
        limits = x_limits,
        breaks = seq(
            1.06,
            1.24,
            by = 0.02
        ),
        labels = function(x) {
            sprintf("%.2f", x)
        },
        expand = expansion(
            mult = c(
                0,
                0
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
            "TN_Right  P/A",
            "TN_Left  P/P"
        ),
        limits = c(
            -0.25,
            1.25
        )
    ) +

    labs(
        title = paste(
            "Independent pooled focal test:",
            "exact expected P/P to P/A sites"
        ),
        subtitle = paste(
            "Only DP2/GQ20 focal sites with the expected",
            "zygosity-aware left-right transition are displayed."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        panel.grid = element_blank(),
        strip.text.y = element_text(
            face = "bold",
            size = 10
        ),
        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = NA
        ),
        plot.title = element_text(
            face = "bold",
            size = 11
        ),
        plot.subtitle = element_text(
            size = 8.5,
            colour = "#444444"
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(
            face = "bold"
        ),
        axis.text.y = element_text(
            face = "bold",
            size = 8.5
        ),
        axis.title.x = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            3,
            12,
            8,
            8
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
                        3.25,
                        1.75
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


cairo_pdf(
    filename = pdf_file,
    width = 11.5,
    height = 10.2
)

draw_combined()

dev.off()


png(
    filename = png_file,
    width = 3450,
    height = 3060,
    res = 300,
    type = "cairo"
)

draw_combined()

dev.off()


cat(
    "WALLBANK-COHORT D/R BA FIGURE\n\n"
)

cat(
    "COMPARATIVE BA COUNTS\n"
)

comparative_summary <- do.call(
    rbind,
    lapply(
        levels(scores$module),
        function(module_name) {
            d <- ba_data[
                ba_data$module == module_name,
                ,
                drop = FALSE
            ]

            data.frame(
                module = module_name,
                eligible_SNPs = nrow(d),
                BA80 = sum(
                    d$balanced_accuracy >= 0.80
                ),
                BA90 = sum(
                    d$balanced_accuracy >= 0.90
                ),
                BA1_among_called = sum(
                    d$perfect_BA
                ),
                BA1_80pct_callable = sum(
                    d$perfect_BA
                    & d$callable_80pct_bool
                ),
                BA1_complete = sum(
                    d$perfect_BA
                    & d$complete_callability_bool
                )
            )
        }
    )
)

print(
    comparative_summary,
    row.names = FALSE
)

cat(
    "\nEXACT FOCAL OVERLAP\n"
)

print(
    exact_summary,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nExact focal sites:\n",
    exact_file,
    "\n\nExact focal summary:\n",
    exact_summary_file,
    "\n",
    sep = ""
)
