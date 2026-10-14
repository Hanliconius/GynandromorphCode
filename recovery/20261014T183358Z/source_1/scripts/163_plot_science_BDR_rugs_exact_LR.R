#!/usr/bin/env Rscript


# ============================================================
# ANALYSIS 163
#
# Science-style B / D / R summary:
#
#   TOP:
#     * separate Band / Dennis / Ray BA rugs
#     * published B, R1, R2, D1, D2 intervals
#     * optix
#
#   BOTTOM:
#     * ONLY strongly associated SNPs that ALSO match the
#       exact expected mixture-adjusted focal transition
#
#       B:  TN_Left A/A -> TN_Right P/A
#       D:  TN_Left P/P -> TN_Right P/A
#       R:  TN_Left P/P -> TN_Right P/A
#
#     * reciprocal-mixture confidence is inherited from 149
#     * BA > 0.99 is required for the relevant trait
#
# Stable analyses are sourced, not modified.
# ============================================================


project <- paste0(
    "/gpfs/automountdir/gpfs/scratch/",
    "martinlab/joe/GynandroCalling"
)


source_file <- file.path(
    project,
    "scripts",
    "159_plot_science_trials_pubregions_cartoons.R"
)


if (!file.exists(source_file)) {

    stop(
        paste(
            "Missing established Analysis-159 source:",
            source_file
        )
    )
}


suppressPackageStartupMessages({

    library(
        ggplot2
    )

    library(
        grid
    )

    library(
        gtable
    )
})


# ============================================================
# LOAD ESTABLISHED DATA OBJECTS
#
# 159 inherits:
#   149 exact focal-pattern machinery
#   158 B/D/R BA tracks
#   established cartoon grobs
#
# Suppress old derivative plotting warnings.
# ============================================================

suppressWarnings(

    source(
        source_file,
        local = TRUE
    )
)


required_objects <- c(
    "all_ba_science",
    "high_B",
    "high_D",
    "high_R",
    "band_correct",
    "dr_correct",
    "trait_colours",
    "band_cartoon",
    "DR_cartoon",
    "science_x_limits",
    "science_x_breaks",
    "optix_start_science",
    "optix_end_science"
)


missing_objects <- setdiff(
    required_objects,
    ls()
)


if (length(missing_objects) > 0L) {

    stop(
        paste(
            "After sourcing 159, missing:",
            paste(
                missing_objects,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# OUTPUTS
# ============================================================

output_dir <- file.path(
    project,
    "analysis",
    "figures",
    "final_B_D_R_combined",
    "science_BDR_rugs_exact_LR"
)


dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


pdf_file <- file.path(
    output_dir,
    paste0(
        "Science_BDR_BA_rugs_",
        "exact_LR_only.pdf"
    )
)


png_file <- file.path(
    output_dir,
    paste0(
        "Science_BDR_BA_rugs_",
        "exact_LR_only.png"
    )
)


audit_file <- file.path(
    output_dir,
    paste0(
        "Science_BDR_BA_rugs_",
        "exact_LR_only.sites.tsv"
    )
)


summary_file <- file.path(
    output_dir,
    paste0(
        "Science_BDR_BA_rugs_",
        "exact_LR_only.summary.tsv"
    )
)


# ============================================================
# COMMON COLOURS
# ============================================================

neutral_dark <- "#66707A"

neutral_line <- "#C8CDD2"

neutral_light <- "#E3E6E9"

text_colour <- "#242424"

optix_colour <- "#3A7CA5"


# ============================================================
# EXACT HIGH-BA FOCAL SITES
#
# These already derive from the exact expected L/R sets:
#
#   high_B:
#       Band exact A/A -> P/A
#       Band BA > 0.99
#
#   high_D:
#       D/R exact P/P -> P/A
#       Dennis BA > 0.99
#
#   high_R:
#       D/R exact P/P -> P/A
#       Ray BA > 0.99
#
# Thus the lower panel has NO merely directional sites and
# NO alternative focal genotype pairs.
# ============================================================

prepare_exact <- function(
    x,
    trait_name,
    left_GT,
    right_GT
) {

    if (nrow(x) == 0L) {

        return(
            data.frame(
                trait = character(),
                position_Mb = numeric(),
                BA = numeric(),
                left_GT = character(),
                right_GT = character(),
                stringsAsFactors = FALSE
            )
        )
    }


    result <- data.frame(
        trait = trait_name,

        position_Mb = suppressWarnings(
            as.numeric(
                x$position_Mb
            )
        ),

        BA = suppressWarnings(
            as.numeric(
                x$BA
            )
        ),

        left_GT = left_GT,

        right_GT = right_GT,

        stringsAsFactors = FALSE
    )


    result <- result[
        is.finite(
            result$position_Mb
        ) &
        is.finite(
            result$BA
        ) &
        result$BA > 0.99 &
        result$position_Mb >=
            science_x_limits[1] &
        result$position_Mb <=
            science_x_limits[2],
        ,
        drop = FALSE
    ]


    result <- result[
        !duplicated(
            paste(
                result$trait,
                round(
                    result$position_Mb *
                        1e6
                ),
                sep = ":"
            )
        ),
        ,
        drop = FALSE
    ]


    result[
        order(
            result$position_Mb
        ),
        ,
        drop = FALSE
    ]
}


exact_B <- prepare_exact(
    high_B,
    "B",
    "A/A",
    "P/A"
)


exact_D <- prepare_exact(
    high_D,
    "D",
    "P/P",
    "P/A"
)


exact_R <- prepare_exact(
    high_R,
    "R",
    "P/P",
    "P/A"
)


exact_sites <- rbind(
    exact_B,
    exact_D,
    exact_R
)


if (nrow(exact_sites) == 0L) {

    stop(
        "No exact high-BA focal sites survived."
    )
}


exact_sites$pos <- as.integer(
    round(
        exact_sites$position_Mb *
            1e6
    )
)


# ============================================================
# UPPER B / D / R BA RUGS
# ============================================================

ba <- all_ba_science


ba$trait <- factor(
    as.character(
        ba$trait
    ),
    levels = c(
        "B",
        "D",
        "R"
    ),
    ordered = TRUE
)


rug_bases <- c(
    B = 4.15,
    D = 3.40,
    R = 2.65
)


ba$rug_base <- unname(
    rug_bases[
        as.character(
            ba$trait
        )
    ]
)


# Same basic height encoding as our preferred Trial 2,
# slightly compressed to leave more room for the L/R panel.

ba$tick_height <- ifelse(

    ba$BA > 0.99,
    0.50,

    ifelse(

        ba$BA >= 0.90,
        0.37,

        ifelse(

            ba$BA >= 0.80,
            0.24,

            0.09
        )
    )
)


ba$rug_top <- (
    ba$rug_base +
    ba$tick_height
)


# Highlight only sites that are BOTH:
#
#   relevant-trait BA > 0.99
#   exact expected focal L/R transition.

exact_upper <- exact_sites


exact_upper$rug_base <- unname(
    rug_bases[
        exact_upper$trait
    ]
)


exact_upper$rug_top <- (
    exact_upper$rug_base +
    0.53
)


# ============================================================
# LOWER EXACT L/R PAIRS
#
# Three independent trait rows.
#
# We explicitly show Left and Right as separate y positions,
# but ONLY for SNPs having the exact expected pair.
# ============================================================

lane_y <- list(

    B = c(
        Left = 1.55,
        Right = 1.27
    ),

    D = c(
        Left = 0.82,
        Right = 0.54
    ),

    R = c(
        Left = 0.09,
        Right = -0.19
    )
)


make_lr_segments <- function(
    data,
    trait_name
) {

    if (nrow(data) == 0L) {

        return(
            data.frame(
                trait = character(),
                position_Mb = numeric(),
                BA = numeric(),
                left_GT = character(),
                right_GT = character(),
                y_left = numeric(),
                y_right = numeric(),
                stringsAsFactors = FALSE
            )
        )
    }


    data.frame(
        trait = trait_name,

        position_Mb = data$position_Mb,

        BA = data$BA,

        left_GT = data$left_GT,

        right_GT = data$right_GT,

        y_left = unname(
            lane_y[[trait_name]][
                "Left"
            ]
        ),

        y_right = unname(
            lane_y[[trait_name]][
                "Right"
            ]
        ),

        stringsAsFactors = FALSE
    )
}


lr_segments <- rbind(

    make_lr_segments(
        exact_B,
        "B"
    ),

    make_lr_segments(
        exact_D,
        "D"
    ),

    make_lr_segments(
        exact_R,
        "R"
    )
)


lr_left <- data.frame(
    trait = lr_segments$trait,
    position_Mb = lr_segments$position_Mb,
    y = lr_segments$y_left,
    stringsAsFactors = FALSE
)


lr_right <- data.frame(
    trait = lr_segments$trait,
    position_Mb = lr_segments$position_Mb,
    y = lr_segments$y_right,
    stringsAsFactors = FALSE
)


# ============================================================
# LIGHT LOWER-PANEL GUIDE LINES
# ============================================================

lane_guides <- do.call(
    rbind,
    lapply(
        names(
            lane_y
        ),
        function(trait_name) {

            data.frame(
                trait = trait_name,

                side = c(
                    "Left",
                    "Right"
                ),

                y = c(
                    lane_y[[trait_name]][
                        "Left"
                    ],
                    lane_y[[trait_name]][
                        "Right"
                    ]
                ),

                stringsAsFactors = FALSE
            )
        }
    )
)


# ============================================================
# PUBLISHED REGIONS
# ============================================================

published_regions <- data.frame(

    region = c(
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

    trait = c(
        "B",
        "R",
        "R",
        "D",
        "D"
    ),

    stringsAsFactors = FALSE
)


published_regions$y <- c(
    5.32,
    5.32,
    5.32,
    5.32,
    5.32
)


published_regions$label_y <- 5.40


# ============================================================
# SMALL DATA SUMMARY / AUDIT
# ============================================================

audit <- exact_sites[
    ,
    c(
        "trait",
        "pos",
        "position_Mb",
        "BA",
        "left_GT",
        "right_GT"
    ),
    drop = FALSE
]


audit <- audit[
    order(
        match(
            audit$trait,
            c(
                "B",
                "D",
                "R"
            )
        ),
        audit$pos
    ),
    ,
    drop = FALSE
]


write.table(
    audit,
    file = audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


summary_table <- data.frame(

    trait = c(
        "B",
        "D",
        "R"
    ),

    exact_high_BA_LR_sites = c(
        nrow(
            exact_B
        ),
        nrow(
            exact_D
        ),
        nrow(
            exact_R
        )
    ),

    expected_left_GT = c(
        "A/A",
        "P/P",
        "P/P"
    ),

    expected_right_GT = c(
        "P/A",
        "P/A",
        "P/A"
    ),

    stringsAsFactors = FALSE
)


write.table(
    summary_table,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# MAIN PLOT
# ============================================================

p <- ggplot() +


    # --------------------------------------------------------
    # Published region bars
    # --------------------------------------------------------

    geom_segment(
        data = published_regions,

        aes(
            x = xmin,
            xend = xmax,
            y = y,
            yend = y,
            colour = trait
        ),

        linewidth = 2.35,
        lineend = "butt",
        show.legend = FALSE
    ) +


    geom_text(
        data = published_regions,

        aes(
            x = (
                xmin +
                xmax
            ) / 2,
            y = label_y,
            label = region,
            colour = trait
        ),

        size = 2.35,
        fontface = "bold",
        vjust = 0,
        show.legend = FALSE
    ) +


    # --------------------------------------------------------
    # optix
    # --------------------------------------------------------

    geom_segment(

        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 4.94,
            yend = 4.94
        ),

        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 2.9,
        lineend = "butt"
    ) +


    annotate(
        "text",

        x = mean(
            c(
                optix_start_science,
                optix_end_science
            )
        ),

        y = 5.02,

        label = "optix",

        colour = optix_colour,
        size = 2.45,
        fontface = "italic",
        vjust = 0
    ) +


    # --------------------------------------------------------
    # Upper rug baselines
    # --------------------------------------------------------

    geom_segment(
        data = data.frame(
            trait = names(
                rug_bases
            ),
            y = unname(
                rug_bases
            ),
            stringsAsFactors = FALSE
        ),

        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = y,
            yend = y
        ),

        inherit.aes = FALSE,
        colour = neutral_line,
        linewidth = 0.27
    ) +


    # --------------------------------------------------------
    # All comparative BA sites
    # --------------------------------------------------------

    geom_segment(
        data = ba,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = rug_base,
            yend = rug_top
        ),

        inherit.aes = FALSE,

        colour = neutral_dark,
        alpha = 0.31,
        linewidth = 0.23
    ) +


    # --------------------------------------------------------
    # Exact L/R + BA>0.99 sites
    #
    # These are the only coloured upper-rug ticks.
    # --------------------------------------------------------

    geom_segment(
        data = exact_upper,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = rug_base,
            yend = rug_top,
            colour = trait
        ),

        linewidth = 0.75,
        alpha = 1,
        lineend = "butt",
        show.legend = FALSE
    ) +


    # --------------------------------------------------------
    # B/D/R labels inside upper panel
    # --------------------------------------------------------

    geom_text(
        data = data.frame(
            trait = c(
                "B",
                "D",
                "R"
            ),

            x = science_x_limits[1] +
                0.0015,

            y = unname(
                rug_bases
            ) +
                0.31,

            stringsAsFactors = FALSE
        ),

        aes(
            x = x,
            y = y,
            label = trait,
            colour = trait
        ),

        inherit.aes = FALSE,

        hjust = 0,
        fontface = "bold",
        size = 3.0,
        show.legend = FALSE
    ) +


    # --------------------------------------------------------
    # Separator between comparative and focal sections
    # --------------------------------------------------------

    geom_segment(

        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = 2.15,
            yend = 2.15
        ),

        inherit.aes = FALSE,
        colour = "#D4D7DA",
        linewidth = 0.35
    ) +


    annotate(
        "text",

        x = science_x_limits[1] +
            0.0015,

        y = 2.23,

        label = "Exact focal L-R matches",

        hjust = 0,
        vjust = 0,

        size = 2.3,
        fontface = "bold",
        colour = "#505050"
    ) +


    # --------------------------------------------------------
    # Lower lane guides
    # --------------------------------------------------------

    geom_segment(
        data = lane_guides,

        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = y,
            yend = y
        ),

        inherit.aes = FALSE,
        colour = neutral_light,
        linewidth = 0.24
    ) +


    # --------------------------------------------------------
    # Exact left-right connectors
    # --------------------------------------------------------

    geom_segment(
        data = lr_segments,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left,
            colour = trait
        ),

        linewidth = 0.42,
        alpha = 0.52,
        show.legend = FALSE
    ) +


    # --------------------------------------------------------
    # Left endpoints
    # --------------------------------------------------------

    geom_point(
        data = lr_left,

        aes(
            x = position_Mb,
            y = y,
            colour = trait
        ),

        shape = 124,
        size = 2.25,
        stroke = 0.55,
        alpha = 1,
        show.legend = FALSE
    ) +


    # --------------------------------------------------------
    # Right endpoints
    # --------------------------------------------------------

    geom_point(
        data = lr_right,

        aes(
            x = position_Mb,
            y = y,
            colour = trait
        ),

        shape = 124,
        size = 2.25,
        stroke = 0.55,
        alpha = 1,
        show.legend = FALSE
    ) +


    scale_colour_manual(
        values = trait_colours,
        drop = FALSE
    ) +


    scale_x_continuous(

        name = "Position on Hmr18 (Mb)",

        breaks = science_x_breaks,

        labels = sprintf(
            "%.2f",
            science_x_breaks
        ),

        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +


    coord_cartesian(
        xlim = science_x_limits,
        ylim = c(
            -0.48,
            5.56
        ),
        clip = "off"
    ) +


    labs(

        title = paste(
            "B, Dennis and Ray association",
            "with exact focal L-R matches"
        ),

        subtitle = paste(
            "BA is shown as compact rugs;",
            "colour identifies BA>0.99 sites",
            "with the exact expected focal genotype pair."
        ),

        caption = paste(
            "Lower panel:",
            "B = A/A -> P/A;",
            "D and R = P/P -> P/A.",
            "Focal calls use reciprocal-mixture posterior >=0.95;",
            "no original-call fallback."
        )
    ) +


    theme_classic(
        base_size = 8
    ) +


    theme(

        plot.title = element_text(
            face = "bold",
            size = 9.5,
            colour = text_colour,
            margin = margin(
                b = 2
            )
        ),

        plot.subtitle = element_text(
            size = 7.1,
            colour = "#505050",
            margin = margin(
                b = 4
            )
        ),

        plot.caption = element_text(
            size = 6.5,
            colour = "#555555",
            hjust = 0,
            margin = margin(
                t = 4
            )
        ),

        axis.title.x = element_text(
            size = 7.5,
            margin = margin(
                t = 4
            )
        ),

        axis.text.x = element_text(
            size = 6.8,
            colour = "#404040"
        ),

        axis.ticks.x = element_line(
            linewidth = 0.25
        ),

        axis.line.x = element_line(
            linewidth = 0.30
        ),

        axis.line.y = element_blank(),

        axis.text.y = element_blank(),

        axis.ticks.y = element_blank(),

        axis.title.y = element_blank(),

        panel.grid = element_blank(),

        legend.position = "none",

        plot.margin = margin(
            8,
            18,
            7,
            10
        )
    )


# ============================================================
# FIXED-WIDTH CARTOON SIDECAR
#
# We deliberately use page space, rather than putting cartoons
# at impossible genomic x coordinates.
#
# B gets the established Band cartoon.
# D and R get the established Dennis/Ray cartoon.
#
# Genotype text makes the exact transition unambiguous.
# ============================================================

plot_y_min <- -0.48

plot_y_max <- 5.56


to_npc <- function(y) {

    (
        y -
        plot_y_min
    ) /
        (
            plot_y_max -
            plot_y_min
        )
}


make_lr_sidecar <- function() {

    B_center <- mean(
        lane_y$B
    )

    D_center <- mean(
        lane_y$D
    )

    R_center <- mean(
        lane_y$R
    )


    grobTree(


        # ====================================================
        # BAND
        # ====================================================

        grobTree(
            band_cartoon,

            vp = viewport(
                x = 0.18,
                y = to_npc(
                    B_center
                ),
                width = 0.26,
                height = 0.105
            )
        ),


        textGrob(
            paste0(
                "B  (n=",
                nrow(
                    exact_B
                ),
                ")"
            ),

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    B_center
                ) +
                    0.021,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 7.1,
                fontface = "bold",
                col = trait_colours[
                    "B"
                ]
            )
        ),


        textGrob(
            "Left A/A  ->  Right P/A",

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    B_center
                ) -
                    0.018,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 6.2,
                col = "#454545"
            )
        ),


        # ====================================================
        # DENNIS
        # ====================================================

        grobTree(
            DR_cartoon,

            vp = viewport(
                x = 0.18,
                y = to_npc(
                    D_center
                ),
                width = 0.26,
                height = 0.105
            )
        ),


        textGrob(
            paste0(
                "D  (n=",
                nrow(
                    exact_D
                ),
                ")"
            ),

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    D_center
                ) +
                    0.021,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 7.1,
                fontface = "bold",
                col = trait_colours[
                    "D"
                ]
            )
        ),


        textGrob(
            "Left P/P  ->  Right P/A",

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    D_center
                ) -
                    0.018,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 6.2,
                col = "#454545"
            )
        ),


        # ====================================================
        # RAY
        # ====================================================

        grobTree(
            DR_cartoon,

            vp = viewport(
                x = 0.18,
                y = to_npc(
                    R_center
                ),
                width = 0.26,
                height = 0.105
            )
        ),


        textGrob(
            paste0(
                "R  (n=",
                nrow(
                    exact_R
                ),
                ")"
            ),

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    R_center
                ) +
                    0.021,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 7.1,
                fontface = "bold",
                col = trait_colours[
                    "R"
                ]
            )
        ),


        textGrob(
            "Left P/P  ->  Right P/A",

            x = unit(
                0.36,
                "npc"
            ),

            y = unit(
                to_npc(
                    R_center
                ) -
                    0.018,
                "npc"
            ),

            hjust = 0,

            gp = gpar(
                fontsize = 6.2,
                col = "#454545"
            )
        )
    )
}


add_sidecar <- function(
    plot,
    sidecar_width = unit(
        1.58,
        "in"
    )
) {

    g <- ggplotGrob(
        plot
    )


    panel_rows <- which(
        grepl(
            "^panel",
            g$layout$name
        )
    )


    if (length(panel_rows) == 0L) {

        stop(
            "Could not identify ggplot panel."
        )
    }


    panel_l <- min(
        g$layout$l[
            panel_rows
        ]
    )


    panel_t <- min(
        g$layout$t[
            panel_rows
        ]
    )


    panel_b <- max(
        g$layout$b[
            panel_rows
        ]
    )


    g <- gtable_add_cols(
        g,
        widths = sidecar_width,
        pos = panel_l - 1L
    )


    g <- gtable_add_grob(
        g,
        grobs = make_lr_sidecar(),
        t = panel_t,
        b = panel_b,
        l = panel_l,
        r = panel_l,
        clip = "off",
        name = "exact-LR-cartoon-sidecar"
    )


    g
}


final_grob <- add_sidecar(
    p
)


# ============================================================
# SAVE
# ============================================================

draw_pdf <- function(
    grob,
    filename,
    width,
    height
) {

    cairo_pdf(
        filename = filename,
        width = width,
        height = height,
        onefile = FALSE
    )


    grid.newpage()


    grid.draw(
        grob
    )


    dev.off()
}


draw_png <- function(
    grob,
    filename,
    width,
    height,
    dpi = 450
) {

    png(
        filename = filename,
        width = width,
        height = height,
        units = "in",
        res = dpi
    )


    grid.newpage()


    grid.draw(
        grob
    )


    dev.off()
}


draw_pdf(
    final_grob,
    pdf_file,
    width = 8.65,
    height = 5.05
)


draw_png(
    final_grob,
    png_file,
    width = 8.65,
    height = 5.05
)


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n"
)

cat(
    "ANALYSIS 163 COMPLETE\n"
)

cat(
    "========================================\n\n"
)


print(
    summary_table,
    row.names = FALSE
)


cat(
    "\nPDF:\n",
    pdf_file,
    "\n",
    sep = ""
)


cat(
    "\nPNG:\n",
    png_file,
    "\n",
    sep = ""
)


cat(
    "\nExact-site audit:\n",
    audit_file,
    "\n",
    sep = ""
)


cat(
    "\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
