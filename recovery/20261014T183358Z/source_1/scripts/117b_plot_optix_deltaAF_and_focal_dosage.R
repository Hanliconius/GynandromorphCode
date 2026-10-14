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
    "optix_deltaAF_haplotype_analysis"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "optix_deltaAF_haplotype_analysis"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

all_file <- file.path(
    analysis_dir,
    "optix_B_D_R_deltaAF_all_SNPs.tsv"
)

selected_file <- file.path(
    analysis_dir,
    "optix_B_D_R_deltaAF_informative_SNPs.tsv"
)

summary_file <- file.path(
    analysis_dir,
    "optix_B_D_R_deltaAF_summary.tsv"
)

pdf_file <- file.path(
    figure_dir,
    "optix_comparative_deltaAF_and_focal_DR_dosage.pdf"
)

png_file <- file.path(
    figure_dir,
    "optix_comparative_deltaAF_and_focal_DR_dosage.png"
)

if (!file.exists(all_file)) {
    stop(
        paste(
            "Missing all-SNP table:",
            all_file
        )
    )
}

if (!file.exists(selected_file)) {
    stop(
        paste(
            "Missing informative-SNP table:",
            selected_file
        )
    )
}

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

all_scores <- read.delim(
    all_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

selected <- read.delim(
    selected_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

summary_table <- read.delim(
    summary_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

module_levels <- c(
    "Band",
    "Dennis",
    "Ray"
)

all_scores$module <- factor(
    all_scores$module,
    levels = module_levels
)

selected$module <- factor(
    selected$module,
    levels = module_levels
)

all_scores$pos <- as.numeric(
    all_scores$pos
)

selected$pos <- as.numeric(
    selected$pos
)

all_scores$position_Mb <- (
    all_scores$pos / 1000000
)

selected$position_Mb <- (
    selected$pos / 1000000
)

all_scores$abs_delta_AF <- as.numeric(
    all_scores$abs_delta_AF
)

selected$abs_delta_AF <- as.numeric(
    selected$abs_delta_AF
)

selected$TN_Left_QC_bool <- as_bool(
    selected$TN_Left_QC
)

selected$TN_Right_QC_bool <- as_bool(
    selected$TN_Right_QC
)

selected$exact_bool <- as_bool(
    selected$exact_left1_right05
)

selected$TN_Left_DR_GT_dosage <- as.numeric(
    selected$TN_Left_DR_GT_dosage
)

selected$TN_Right_DR_GT_dosage <- as.numeric(
    selected$TN_Right_DR_GT_dosage
)

selected$TN_Left_DR_AD_fraction <- as.numeric(
    selected$TN_Left_DR_AD_fraction
)

selected$TN_Right_DR_AD_fraction <- as.numeric(
    selected$TN_Right_DR_AD_fraction
)


background <- all_scores[
    is.finite(
        all_scores$abs_delta_AF
    ),
    ,
    drop = FALSE
]


focal_pairs <- selected[
    selected$TN_Left_QC_bool
    & selected$TN_Right_QC_bool
    & is.finite(
        selected$TN_Left_DR_AD_fraction
    )
    & is.finite(
        selected$TN_Right_DR_AD_fraction
    ),
    ,
    drop = FALSE
]

exact_pairs <- focal_pairs[
    focal_pairs$exact_bool,
    ,
    drop = FALSE
]


dose_label <- function(x) {
    result <- rep(
        "Other/no call",
        length(x)
    )

    result[
        is.finite(x)
        & abs(x - 0) < 1e-8
    ] <- "non-DR/non-DR"

    result[
        is.finite(x)
        & abs(x - 0.5) < 1e-8
    ] <- "DR/non-DR"

    result[
        is.finite(x)
        & abs(x - 1) < 1e-8
    ] <- "DR/DR"

    result
}


focal_long <- rbind(
    data.frame(
        module = focal_pairs$module,
        position_Mb = focal_pairs$position_Mb,
        side = "TN_Left",
        DR_AD_fraction = (
            focal_pairs$TN_Left_DR_AD_fraction
        ),
        DR_GT_dosage = (
            focal_pairs$TN_Left_DR_GT_dosage
        ),
        genotype_class = dose_label(
            focal_pairs$TN_Left_DR_GT_dosage
        ),
        stringsAsFactors = FALSE
    ),
    data.frame(
        module = focal_pairs$module,
        position_Mb = focal_pairs$position_Mb,
        side = "TN_Right",
        DR_AD_fraction = (
            focal_pairs$TN_Right_DR_AD_fraction
        ),
        DR_GT_dosage = (
            focal_pairs$TN_Right_DR_GT_dosage
        ),
        genotype_class = dose_label(
            focal_pairs$TN_Right_DR_GT_dosage
        ),
        stringsAsFactors = FALSE
    )
)

focal_long$module <- factor(
    focal_long$module,
    levels = module_levels
)

focal_long$side <- factor(
    focal_long$side,
    levels = c(
        "TN_Left",
        "TN_Right"
    )
)

focal_long$genotype_class <- factor(
    focal_long$genotype_class,
    levels = c(
        "non-DR/non-DR",
        "DR/non-DR",
        "DR/DR",
        "Other/no call"
    )
)


module_colours <- c(
    "Band" = "#3C9DC7",
    "Dennis" = "#B51A7D",
    "Ray" = "#E88924"
)

side_colours <- c(
    "TN_Left" = "#0072B2",
    "TN_Right" = "#D55E00"
)

genotype_shapes <- c(
    "non-DR/non-DR" = 15,
    "DR/non-DR" = 18,
    "DR/DR" = 16,
    "Other/no call" = 4
)

x_limits <- c(
    1.060,
    1.240
)

gene_annotation <- data.frame(
    module = factor(
        module_levels,
        levels = module_levels
    ),
    xmin = 1088602 / 1000000,
    xmax = 1097846 / 1000000
)

intervals <- data.frame(
    module = factor(
        c(
            "Band",
            "Dennis",
            "Dennis",
            "Ray",
            "Ray"
        ),
        levels = module_levels
    ),
    interval = c(
        "B",
        "D1",
        "D2",
        "R1",
        "R2"
    ),
    xmin = c(
        1112651,
        1215210,
        1218557,
        1156746,
        1174547
    ) / 1000000,
    xmax = c(
        1114314,
        1216124,
        1221494,
        1160220,
        1187833
    ) / 1000000
)


comparative_plot <- ggplot() +

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
            ymin = 0,
            ymax = 0.045
        ),
        inherit.aes = FALSE,
        fill = "#777777",
        alpha = 0.24
    ) +

    geom_text(
        data = intervals,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 0.065,
            label = interval
        ),
        inherit.aes = FALSE,
        colour = "#555555",
        size = 2.5
    ) +

    geom_hline(
        yintercept = 0.80,
        colour = "#777777",
        linewidth = 0.35,
        linetype = "dotted"
    ) +

    geom_point(
        data = background,
        aes(
            x = position_Mb,
            y = abs_delta_AF
        ),
        colour = "#BDBDBD",
        alpha = 0.16,
        size = 0.34,
        stroke = 0
    ) +

    geom_point(
        data = selected,
        aes(
            x = position_Mb,
            y = abs_delta_AF,
            colour = module
        ),
        shape = 25,
        size = 1.75,
        stroke = 0.20,
        alpha = 0.82
    ) +

    geom_text(
        data = gene_annotation,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 1.015,
            label = "optix"
        ),
        inherit.aes = FALSE,
        colour = "#277DA1",
        fontface = "italic",
        size = 2.8
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        values = module_colours,
        guide = "none"
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
            sprintf(
                "%.2f",
                x
            )
        },
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_y_continuous(
        name = expression(
            abs(Delta * " allele frequency")
        ),
        limits = c(
            0,
            1.03
        ),
        breaks = c(
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0
        )
    ) +

    labs(
        title = paste(
            "Comparative pattern-associated variation",
            "across optix"
        ),
        subtitle = paste(
            "Triangles have |delta AF| >= 0.80 and",
            "at least 80% callability in both comparative groups.",
            "TN_Left and TN_Right were not used for SNP selection."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#ECECEC",
            linewidth = 0.25
        ),
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
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9,
            colour = "#444444"
        ),
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

    geom_hline(
        yintercept = 0.50,
        colour = side_colours[
            "TN_Right"
        ],
        linewidth = 0.30,
        linetype = "dashed",
        alpha = 0.55
    ) +

    geom_hline(
        yintercept = 1.00,
        colour = side_colours[
            "TN_Left"
        ],
        linewidth = 0.30,
        linetype = "dashed",
        alpha = 0.55
    ) +

    geom_segment(
        data = focal_pairs,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = TN_Right_DR_AD_fraction,
            yend = TN_Left_DR_AD_fraction
        ),
        colour = "#8A8A8A",
        linewidth = 0.24,
        alpha = 0.30
    ) +

    geom_segment(
        data = exact_pairs,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = TN_Right_DR_AD_fraction,
            yend = TN_Left_DR_AD_fraction
        ),
        colour = "#009E73",
        linewidth = 0.38,
        alpha = 0.58
    ) +

    geom_point(
        data = focal_long,
        aes(
            x = position_Mb,
            y = DR_AD_fraction,
            colour = side,
            shape = genotype_class
        ),
        size = 1.55,
        alpha = 0.82,
        stroke = 0.25
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        name = "Pooled side",
        values = side_colours
    ) +

    scale_shape_manual(
        name = "Called genotype",
        values = genotype_shapes,
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
            sprintf(
                "%.2f",
                x
            )
        },
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_y_continuous(
        name = "DR-associated allele fraction",
        limits = c(
            -0.03,
            1.03
        ),
        breaks = c(
            0,
            0.5,
            1
        ),
        labels = c(
            "0\nnon-DR/non-DR",
            "0.5\nDR/non-DR",
            "1\nDR/DR"
        )
    ) +

    labs(
        title = paste(
            "Independent pooled left-right haplotype evidence"
        ),
        subtitle = paste(
            "Points show read-level DR-allele fraction;",
            "shape shows the called genotype.",
            "Green connectors mark called left=DR/DR",
            "and right=DR/non-DR."
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
        panel.grid.major.y = element_line(
            colour = "#E2E2E2",
            linewidth = 0.30
        ),
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
        axis.title = element_text(
            face = "bold"
        ),
        axis.text.y = element_text(
            size = 8
        ),
        plot.margin = margin(
            3,
            12,
            8,
            8
        )
    )


comparative_grob <- ggplotGrob(
    comparative_plot
)

focal_grob <- ggplotGrob(
    focal_plot
)

common_widths <- unit.pmax(
    comparative_grob$widths,
    focal_grob$widths
)

comparative_grob$widths <- common_widths
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
                        3,
                        2.3
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
        comparative_grob
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
    height = 11.2
)

draw_combined()

dev.off()


png(
    filename = png_file,
    width = 3450,
    height = 3360,
    res = 300,
    type = "cairo"
)

draw_combined()

dev.off()


cat(
    "OPTIX DELTA-AF AND FOCAL DR-DOSAGE FIGURE\n\n"
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n",
    sep = ""
)
