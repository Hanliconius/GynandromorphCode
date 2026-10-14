#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project,
    "analysis",
    "B_vs_DR_repolarization",
    "checkpoints"
)

figure_dir <- file.path(
    project,
    "analysis",
    "figures",
    "B_vs_DR_repolarization"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

scores_file <- file.path(
    checkpoint_dir,
    "B_vs_DR_broad_optix_site_scores.tsv.gz"
)

focal_file <- file.path(
    checkpoint_dir,
    "B_vs_DR_supporting_sites_adjusted_LR.tsv"
)

if (!file.exists(scores_file)) {
    stop(
        paste(
            "Missing:",
            scores_file
        )
    )
}

if (!file.exists(focal_file)) {
    stop(
        paste(
            "Missing:",
            focal_file
        )
    )
}


# ============================================================
# READ
# ============================================================

scores <- read.delim(
    gzfile(scores_file),
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

focal <- read.delim(
    focal_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)


scores$position_Mb <- as.numeric(
    scores$position_Mb
)

scores$balanced_accuracy <- as.numeric(
    scores$balanced_accuracy
)

scores$score_eligible <- (
    tolower(
        as.character(
            scores$score_eligible
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
)

scores$BDR_support <- (
    tolower(
        as.character(
            scores$BDR_support
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
)

focal$position_Mb <- as.numeric(
    focal$position_Mb
)

focal$focal_pair_confident <- (
    tolower(
        as.character(
            focal$focal_pair_confident
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
)


# ============================================================
# REGION ANNOTATIONS
# ============================================================

x_limits <- c(
    1.060,
    1.240
)

x_breaks <- seq(
    1.06,
    1.24,
    by = 0.02
)

optix_start <- 1088602 / 1e6
optix_end <- 1097846 / 1e6

published_regions <- data.frame(
    label = c(
        "B",
        "R1",
        "R2",
        "D1",
        "D2"
    ),

    xmin = c(
        1112651,
        1156746,
        1174547,
        1215210,
        1218557
    ) / 1e6,

    xmax = c(
        1114314,
        1160220,
        1187833,
        1216124,
        1221494
    ) / 1e6,

    stringsAsFactors = FALSE
)


# ============================================================
# PALETTE
# ============================================================

support_colours <- c(
    "B" = "#0072B2",
    "DR" = "#CC79A7",
    "neither" = "#CFD3D8"
)

genotype_colours <- c(
    "B/B" = "#0072B2",
    "Het" = "#6F6F6F",
    "DR/DR" = "#CC79A7",
    "No call" = "#D9D9D9"
)

side_shapes <- c(
    "TN_Left" = 16,
    "TN_Right" = 17
)


# ============================================================
# COMPARATIVE BA PANEL
# ============================================================

plot_scores <- scores[
    is.finite(
        scores$balanced_accuracy
    ),
    ,
    drop = FALSE
]

plot_scores$ALT_class <- factor(
    plot_scores$ALT_class,
    levels = c(
        "neither",
        "B",
        "DR"
    )
)


ba_plot <- ggplot() +

    geom_segment(
        data = published_regions,
        aes(
            x = xmin,
            xend = xmax,
            y = 1.035,
            yend = 1.035
        ),
        inherit.aes = FALSE,
        linewidth = 2.1,
        colour = "#6C727A",
        lineend = "butt"
    ) +

    geom_text(
        data = published_regions,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 1.052,
            label = label
        ),
        inherit.aes = FALSE,
        size = 2.5,
        colour = "#4C5158"
    ) +

    geom_segment(
        aes(
            x = optix_start,
            xend = optix_end,
            y = 1.085,
            yend = 1.085
        ),
        inherit.aes = FALSE,
        colour = "#3A7CA5",
        linewidth = 3.0,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = mean(
            c(
                optix_start,
                optix_end
            )
        ),
        y = 1.105,
        label = "optix",
        fontface = "italic",
        colour = "#3A7CA5",
        size = 2.6
    ) +

    geom_hline(
        yintercept = 0.80,
        linewidth = 0.30,
        linetype = "dashed",
        colour = "#777777"
    ) +

    geom_point(
        data = plot_scores[
            !plot_scores$BDR_support,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        colour = "#CFD3D8",
        alpha = 0.42,
        size = 0.42,
        stroke = 0
    ) +

    geom_point(
        data = plot_scores[
            plot_scores$BDR_support,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = ALT_class
        ),
        alpha = 0.92,
        size = 0.72,
        stroke = 0
    ) +

    scale_colour_manual(
        values = support_colours,
        breaks = c(
            "B",
            "DR"
        ),
        labels = c(
            "ALT is B-like",
            "ALT is DR-like"
        ),
        name = NULL
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = x_breaks,
        labels = NULL,
        expand = c(
            0,
            0
        )
    ) +

    scale_y_continuous(
        limits = c(
            0.48,
            1.12
        ),
        breaks = c(
            0.5,
            0.8,
            1.0
        ),
        expand = c(
            0,
            0
        )
    ) +

    labs(
        title = "Unified B-like versus Dennis/Ray-like association",
        subtitle = paste(
            "Each informative SNP is polarized to a B-like and DR-like allele;",
            "colored sites have BA >= 0.80 and agree in direction for Dennis and Ray."
        ),
        y = "Balanced accuracy"
    ) +

    theme_classic(
        base_size = 9
    ) +

    theme(
        legend.position = "top",
        legend.justification = "left",

        plot.title = element_text(
            face = "bold",
            size = 10
        ),

        plot.subtitle = element_text(
            size = 7.4,
            colour = "#505050"
        ),

        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),

        panel.grid = element_blank(),

        plot.margin = margin(
            8,
            14,
            1,
            10
        )
    )


# ============================================================
# ONE FOCAL B / HET / DR TRACK
# ============================================================

focal_confident <- focal[
    focal$focal_pair_confident,
    ,
    drop = FALSE
]

state_y <- c(
    "DR/DR" = 0,
    "Het" = 1,
    "B/B" = 2
)

left <- data.frame(
    site_key = focal_confident$site_key,
    position_Mb = focal_confident$position_Mb,
    side = "TN_Left",
    state = focal_confident$TN_Left_BDR_GT,
    stringsAsFactors = FALSE
)

right <- data.frame(
    site_key = focal_confident$site_key,
    position_Mb = focal_confident$position_Mb,
    side = "TN_Right",
    state = focal_confident$TN_Right_BDR_GT,
    stringsAsFactors = FALSE
)

focal_long <- rbind(
    left,
    right
)

focal_long <- focal_long[
    focal_long$state %in% names(
        state_y
    ),
    ,
    drop = FALSE
]

focal_long$y <- unname(
    state_y[
        focal_long$state
    ]
)

# Tiny display-only x offset so identical L/R states remain visible.
focal_long$x_display <- (
    focal_long$position_Mb
    + ifelse(
        focal_long$side ==
            "TN_Left",
        -0.000035,
        0.000035
    )
)

focal_long$state <- factor(
    focal_long$state,
    levels = c(
        "B/B",
        "Het",
        "DR/DR"
    )
)

focal_long$side <- factor(
    focal_long$side,
    levels = c(
        "TN_Left",
        "TN_Right"
    )
)


segments <- focal_confident[
    focal_confident$TN_Left_BDR_GT %in%
        names(state_y)
    &
    focal_confident$TN_Right_BDR_GT %in%
        names(state_y),
    ,
    drop = FALSE
]

segments$left_y <- unname(
    state_y[
        segments$TN_Left_BDR_GT
    ]
)

segments$right_y <- unname(
    state_y[
        segments$TN_Right_BDR_GT
    ]
)


focal_plot <- ggplot() +

    geom_hline(
        yintercept = c(
            0,
            1,
            2
        ),
        colour = "#E0E2E5",
        linewidth = 0.35
    ) +

    geom_segment(
        data = segments,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = left_y,
            yend = right_y
        ),
        inherit.aes = FALSE,
        colour = "#B7BDC4",
        alpha = 0.70,
        linewidth = 0.30
    ) +

    geom_point(
        data = focal_long,
        aes(
            x = x_display,
            y = y,
            colour = state,
            shape = side
        ),
        size = 1.35,
        alpha = 0.90,
        stroke = 0.20
    ) +

    scale_colour_manual(
        values = genotype_colours,
        name = "Focal genotype"
    ) +

    scale_shape_manual(
        values = side_shapes,
        name = NULL
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = x_breaks,
        labels = sprintf(
            "%.2f",
            x_breaks
        ),
        expand = c(
            0,
            0
        )
    ) +

    scale_y_continuous(
        limits = c(
            -0.25,
            2.25
        ),
        breaks = c(
            2,
            1,
            0
        ),
        labels = c(
            "B/B",
            "Het",
            "DR/DR"
        ),
        expand = c(
            0,
            0
        )
    ) +

    labs(
        title = "Mixture-adjusted focal genotypes at B/DR-supporting SNPs",
        subtitle = paste(
            "Circle = TN_Left; triangle = TN_Right.",
            "Only reciprocal-mixture calls with pair posterior >= 0.95 are shown."
        ),
        x = "Position on Hmr18 (Mb)",
        y = NULL
    ) +

    theme_classic(
        base_size = 9
    ) +

    theme(
        legend.position = "top",
        legend.justification = "left",

        plot.title = element_text(
            face = "bold",
            size = 9
        ),

        plot.subtitle = element_text(
            size = 7.3,
            colour = "#505050"
        ),

        axis.text.y = element_text(
            face = "bold",
            size = 8
        ),

        axis.ticks.y = element_blank(),

        panel.grid = element_blank(),

        plot.margin = margin(
            2,
            14,
            8,
            10
        )
    )


# ============================================================
# ALIGN
# ============================================================

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
                        1.75,
                        1.0
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


# ============================================================
# SAVE
# ============================================================

pdf_file <- file.path(
    figure_dir,
    "B_vs_DR_repolarized_BA_and_single_focal_genotype_track.pdf"
)

png_file <- file.path(
    figure_dir,
    "B_vs_DR_repolarized_BA_and_single_focal_genotype_track.png"
)

cairo_pdf(
    filename = pdf_file,
    width = 7.4,
    height = 5.2
)

draw_combined()

dev.off()


png(
    filename = png_file,
    width = 7.4,
    height = 5.2,
    units = "in",
    res = 450
)

draw_combined()

dev.off()


# ============================================================
# FOCAL PAIR SUMMARY
# ============================================================

pair_summary <- as.data.frame(
    table(
        focal_confident$focal_state_pair
    ),
    stringsAsFactors = FALSE
)

names(pair_summary) <- c(
    "focal_state_pair",
    "SNPs"
)

pair_summary <- pair_summary[
    order(
        -pair_summary$SNPs
    ),
    ,
    drop = FALSE
]

pair_summary_file <- file.path(
    checkpoint_dir,
    "B_vs_DR_focal_state_pair_summary.tsv"
)

write.table(
    pair_summary,
    file = pair_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


cat(
    "\n========================================\n",
    "ANALYSIS 161 B-vs-DR REPOLARIZATION COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "Confident focal B/DR-supporting sites: ",
    nrow(
        focal_confident
    ),
    "\n\n",
    sep = ""
)

print(
    pair_summary,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nPair summary:\n",
    pair_summary_file,
    "\n",
    sep = ""
)
