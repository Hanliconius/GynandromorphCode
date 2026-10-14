#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

script_dir <- file.path(
    project,
    "scripts"
)

output_dir <- file.path(
    project,
    "analysis",
    "figures",
    "final_B_D_R_combined"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dr_script <- file.path(
    script_dir,
    "115b_plot_D_R_Wallbank_fixed_dominant_adjusted_LR.R"
)

band_script <- file.path(
    script_dir,
    "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R"
)

output_pdf <- file.path(
    output_dir,
    "Wallbank_Dennis_Ray_plus_restricted_Band_correct_pattern_ribbons.pdf"
)

dr_audit_file <- file.path(
    output_dir,
    "DR_broad_correct_pattern_ribbon_sites.tsv"
)

band_audit_file <- file.path(
    output_dir,
    "Band_broad_correct_pattern_ribbon_sites.tsv"
)

required <- c(
    dr_script,
    band_script
)

missing <- required[
    !file.exists(required)
]

if (length(missing) > 0) {
    stop(
        paste(
            "Missing required script(s):",
            paste(
                missing,
                collapse = "\n"
            )
        )
    )
}


# ============================================================
# LOAD EXISTING ANALYSIS ENVIRONMENT WITHOUT WRITING ITS PDF
# ============================================================

load_plot_environment <- function(
    script
) {

    lines <- readLines(
        script,
        warn = FALSE
    )

    device_start <- which(
        trimws(lines) == "cairo_pdf("
    )

    if (length(device_start) != 1) {
        stop(
            paste(
                basename(script),
                ": expected exactly one cairo_pdf() call; found",
                length(device_start)
            )
        )
    }

    env <- new.env(
        parent = globalenv()
    )

    eval(
        parse(
            text = paste(
                lines[
                    seq_len(
                        device_start[[1]] - 1
                    )
                ],
                collapse = "\n"
            )
        ),
        envir = env
    )

    if (!exists(
        "ba_plot",
        envir = env,
        inherits = FALSE
    )) {
        stop(
            paste(
                basename(script),
                "did not create ba_plot."
            )
        )
    }

    env
}


cat(
    "\nLoading established D/R analysis...\n"
)

dr_env <- load_plot_environment(
    dr_script
)

cat(
    "Loading established Band analysis...\n"
)

band_env <- load_plot_environment(
    band_script
)


# ============================================================
# REMOVE OBSOLETE BA=1 HIGHLIGHT SYMBOLS
#
# In both established plots the special highlight layer contains
# 19 points:
#
#   D/R  = 16 Dennis + 3 Ray selected markers
#   Band = 19 perfect-BA highlighted points
#
# Remove only GeomPoint layers carrying exactly those 19 rows.
# Ordinary BA clouds are left completely untouched.
# ============================================================

remove_19point_highlight_layers <- function(
    plot_object,
    plot_name
) {

    layer_summary <- data.frame(
        layer = seq_along(
            plot_object$layers
        ),
        geom = vapply(
            plot_object$layers,
            function(x) {
                class(
                    x$geom
                )[1]
            },
            character(1)
        ),
        n_rows = vapply(
            plot_object$layers,
            function(x) {

                d <- x$data

                if (
                    is.null(d) ||
                    inherits(
                        d,
                        "waiver"
                    ) ||
                    !is.data.frame(d)
                ) {
                    return(
                        NA_integer_
                    )
                }

                nrow(d)
            },
            integer(1)
        ),
        stringsAsFactors = FALSE
    )

    cat(
        "\n",
        plot_name,
        " layer audit before removing obsolete highlights:\n",
        sep = ""
    )

    print(
        layer_summary
    )

    candidates <- which(
        layer_summary$geom == "GeomPoint" &
        layer_summary$n_rows == 19
    )

    if (length(candidates) == 0) {
        stop(
            paste(
                plot_name,
                ": could not identify the 19-point",
                "BA-highlight layer."
            )
        )
    }

    cat(
        "\nRemoving ",
        plot_name,
        " GeomPoint layer(s): ",
        paste(
            candidates,
            collapse = ", "
        ),
        "\n",
        sep = ""
    )

    plot_object$layers <- plot_object$layers[
        -candidates
    ]

    plot_object
}


dr_ba_plot <- remove_19point_highlight_layers(
    dr_env$ba_plot,
    "D/R"
)

band_ba_plot <- remove_19point_highlight_layers(
    band_env$ba_plot,
    "Band"
)


# ============================================================
# UPDATE ONLY TEXT THAT REFERRED TO THE OBSOLETE HIGHLIGHTS
# ============================================================

dr_ba_plot <- dr_ba_plot +

    labs(
        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Comparative balanced accuracy is shown directly;",
            "TN_Left and TN_Right were excluded from SNP selection."
        )
    )


band_ba_plot <- band_ba_plot +

    labs(
        subtitle = paste(
            "Hmr18:1.060-1.240 Mb.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right;",
            "correct-pattern focal genotypes are displayed below."
        )
    )


# ============================================================
# BROAD REGION
# ============================================================

chromosome <- "Hmr18"

region_start <- 1060000
region_end <- 1240000

x_limits <- c(
    region_start,
    region_end
) / 1e6

x_breaks <- seq(
    1.06,
    1.24,
    by = 0.02
)


# ============================================================
# D/R RIBBON
#
# Use the FULL reciprocal-mixture table, not the 19 selected
# comparative sites and not the frozen module windows.
#
# Correct inheritance pattern:
#
#     TN_Left  = raw 1/1 = biological P/P
#     TN_Right = raw 0/1 = biological P/A
#
# Require pair posterior >= 0.95 exactly as in the established
# adjusted focal analysis.
# ============================================================

if (!exists(
    "mixture",
    envir = dr_env,
    inherits = FALSE
)) {
    stop(
        "D/R environment does not contain the full mixture table."
    )
}

dr_mix <- dr_env$mixture

required_dr_columns <- c(
    "chrom",
    "pos",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior"
)

missing_dr_columns <- setdiff(
    required_dr_columns,
    names(dr_mix)
)

if (length(missing_dr_columns) > 0) {
    stop(
        paste(
            "D/R mixture table is missing:",
            paste(
                missing_dr_columns,
                collapse = ", "
            )
        )
    )
}

dr_mix$pos <- suppressWarnings(
    as.numeric(
        dr_mix$pos
    )
)

dr_mix$pair_posterior <- suppressWarnings(
    as.numeric(
        dr_mix$pair_posterior
    )
)

dr_correct <- dr_mix[
    dr_mix$chrom == chromosome &
    is.finite(
        dr_mix$pos
    ) &
    dr_mix$pos >= region_start &
    dr_mix$pos <= region_end &
    is.finite(
        dr_mix$pair_posterior
    ) &
    dr_mix$pair_posterior >= 0.95 &
    dr_mix$mixture_left_GT == "1/1" &
    dr_mix$mixture_right_GT == "0/1",
    ,
    drop = FALSE
]

# One plotted connector per genomic position.
dr_correct <- dr_correct[
    !duplicated(
        paste(
            dr_correct$chrom,
            dr_correct$pos,
            sep = ":"
        )
    ),
    ,
    drop = FALSE
]

dr_correct <- dr_correct[
    order(
        dr_correct$pos
    ),
    ,
    drop = FALSE
]

if (nrow(dr_correct) == 0) {
    stop(
        "No broad-region D/R P/P -> P/A sites survived."
    )
}

dr_correct$position_Mb <- (
    dr_correct$pos /
    1e6
)

dr_long <- rbind(
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "P/P",
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        stringsAsFactors = FALSE
    )
)


# ============================================================
# BAND RIBBON
#
# The analysis-142 environment already contains the broad
# adjusted focal calls.  Select only the exact Band inheritance
# pattern across the WHOLE region:
#
#     TN_Left  = biological A/A
#     TN_Right = biological P/A
#
# This is intentionally not restricted to the Band window.
# ============================================================

if (!exists(
    "d",
    envir = band_env,
    inherits = FALSE
)) {
    stop(
        "Band environment does not contain d."
    )
}

band_all <- band_env$d

required_band_columns <- c(
    "chrom",
    "pos",
    "adjusted_confident",
    "exact_expected_LR",
    "left_biological_GT",
    "right_biological_GT"
)

missing_band_columns <- setdiff(
    required_band_columns,
    names(band_all)
)

if (length(missing_band_columns) > 0) {
    stop(
        paste(
            "Band table is missing:",
            paste(
                missing_band_columns,
                collapse = ", "
            )
        )
    )
}

band_all$pos <- suppressWarnings(
    as.numeric(
        band_all$pos
    )
)

band_correct <- band_all[
    band_all$chrom == chromosome &
    is.finite(
        band_all$pos
    ) &
    band_all$pos >= region_start &
    band_all$pos <= region_end &
    band_all$adjusted_confident %in% TRUE &
    band_all$exact_expected_LR %in% TRUE &
    band_all$left_biological_GT == "A/A" &
    band_all$right_biological_GT == "P/A",
    ,
    drop = FALSE
]

band_correct <- band_correct[
    !duplicated(
        paste(
            band_correct$chrom,
            band_correct$pos,
            sep = ":"
        )
    ),
    ,
    drop = FALSE
]

band_correct <- band_correct[
    order(
        band_correct$pos
    ),
    ,
    drop = FALSE
]

if (nrow(band_correct) == 0) {
    stop(
        "No broad-region Band A/A -> P/A sites survived."
    )
}

band_correct$position_Mb <- (
    band_correct$pos /
    1e6
)

band_long <- rbind(
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "A/A",
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        stringsAsFactors = FALSE
    )
)


# ============================================================
# GENOTYPE PALETTE
# ============================================================

genotype_colours <- c(
    "A/A" = "#E69F00",
    "P/A" = "#CC79A7",
    "P/P" = "#0072B2"
)


# ============================================================
# SIMPLE CORRECT-PATTERN RIBBON
#
# Deliberately:
#   * no module-window filtering
#   * no window rectangles
#   * no incorrect-pattern SNPs
#   * no no-call SNPs
#
# Every visible connector therefore has the stated inheritance
# pattern.
# ============================================================

make_correct_pattern_ribbon <- function(
    wide_data,
    long_data,
    title,
    subtitle,
    caption
) {

    ggplot() +

        geom_segment(
            data = wide_data,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = "#009E73",
            alpha = 0.46,
            linewidth = 0.34
        ) +

        geom_point(
            data = long_data,
            aes(
                x = position_Mb,
                y = y,
                colour = biological_GT
            ),
            shape = 124,
            size = 2.35,
            alpha = 0.90
        ) +

        scale_colour_manual(
            values = genotype_colours,
            drop = FALSE,
            name = "Biological genotype"
        ) +

        scale_x_continuous(
            name = "Position across the broad Optix interval (Mb)",
            limits = x_limits,
            breaks = x_breaks,
            labels = sprintf(
                "%.2f",
                x_breaks
            ),
            expand = expansion(
                mult = c(
                    0.002,
                    0.002
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
                -0.30,
                1.30
            )
        ) +

        labs(
            title = title,
            subtitle = subtitle,
            caption = caption
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
                colour = "#444444",
                margin = margin(
                    b = 5
                )
            ),
            plot.caption = element_text(
                size = 7.8,
                colour = "#555555",
                hjust = 0,
                margin = margin(
                    t = 5
                )
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
                9,
                10
            )
        )
}


dr_ribbon <- make_correct_pattern_ribbon(
    wide_data = dr_correct,
    long_data = dr_long,
    title = "Mixture-adjusted D/R-compatible left-right genotypes",
    subtitle = paste0(
        "All ",
        nrow(dr_correct),
        " SNPs across Hmr18:1.060-1.240 Mb with ",
        "TN_Left=P/P and TN_Right=P/A; ",
        "reciprocal-mixture pair posterior >=0.95."
    ),
    caption = paste(
        "Sites are selected only by the focal left-right inheritance pattern",
        "across the full plotted region; Dennis/Ray interval membership is not used."
    )
)


band_ribbon <- make_correct_pattern_ribbon(
    wide_data = band_correct,
    long_data = band_long,
    title = "Mixture-adjusted Band-compatible left-right genotypes",
    subtitle = paste0(
        "All ",
        nrow(band_correct),
        " SNPs across Hmr18:1.060-1.240 Mb with ",
        "TN_Left=A/A and TN_Right=P/A."
    ),
    caption = paste(
        "Sites are selected only by the focal left-right inheritance pattern",
        "across the full plotted region; Band-window membership is not used."
    )
)


# ============================================================
# ALIGN GROBS
# ============================================================

align_pair <- function(
    upper_plot,
    lower_plot
) {

    upper <- ggplotGrob(
        upper_plot
    )

    lower <- ggplotGrob(
        lower_plot
    )

    common_widths <- unit.pmax(
        upper$widths,
        lower$widths
    )

    upper$widths <- common_widths
    lower$widths <- common_widths

    list(
        upper = upper,
        lower = lower
    )
}


dr_grobs <- align_pair(
    dr_ba_plot,
    dr_ribbon
)

band_grobs <- align_pair(
    band_ba_plot,
    band_ribbon
)


# ============================================================
# DRAWING HELPER
# ============================================================

draw_pair <- function(
    upper_grob,
    lower_grob,
    panel_label,
    upper_fraction
) {

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 2,
                ncol = 1,
                heights = unit(
                    c(
                        upper_fraction,
                        1 - upper_fraction
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
        upper_grob
    )

    grid.text(
        panel_label,
        x = unit(
            0.006,
            "npc"
        ),
        y = unit(
            0.995,
            "npc"
        ),
        just = c(
            "left",
            "top"
        ),
        gp = gpar(
            fontface = "bold",
            fontsize = 15
        )
    )

    popViewport()

    pushViewport(
        viewport(
            layout.pos.row = 2,
            layout.pos.col = 1
        )
    )

    grid.draw(
        lower_grob
    )

    popViewport()

    popViewport()
}


# ============================================================
# OUTPUT
# ============================================================

grDevices::cairo_pdf(
    filename = output_pdf,
    width = 11.7,
    height = 18.6,
    onefile = TRUE
)

grid.newpage()

pushViewport(
    viewport(
        layout = grid.layout(
            nrow = 2,
            ncol = 1,
            heights = unit(
                c(
                    10.2,
                    8.4
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

draw_pair(
    dr_grobs$upper,
    dr_grobs$lower,
    "A",
    upper_fraction = 0.67
)

popViewport()

pushViewport(
    viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
    )
)

draw_pair(
    band_grobs$upper,
    band_grobs$lower,
    "B",
    upper_fraction = 0.64
)

popViewport()

popViewport()

dev.off()


# ============================================================
# AUDITS
# ============================================================

write.table(
    dr_correct,
    file = dr_audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    band_correct,
    file = band_audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\n========================================\n",
    "ANALYSIS 145 FIGURE COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "D/R correct-pattern broad-region SNPs: ",
    nrow(dr_correct),
    "\n",
    sep = ""
)

cat(
    "Band correct-pattern broad-region SNPs: ",
    nrow(band_correct),
    "\n\n",
    sep = ""
)

cat(
    "D/R position range: ",
    min(dr_correct$pos),
    " - ",
    max(dr_correct$pos),
    "\n",
    sep = ""
)

cat(
    "Band position range: ",
    min(band_correct$pos),
    " - ",
    max(band_correct$pos),
    "\n\n",
    sep = ""
)

cat(
    "PDF:\n",
    output_pdf,
    "\n\n",
    "D/R audit:\n",
    dr_audit_file,
    "\n\n",
    "Band audit:\n",
    band_audit_file,
    "\n",
    sep = ""
)

if (
    !file.exists(output_pdf) ||
    file.info(output_pdf)$size == 0
) {
    stop(
        "Output PDF was not created."
    )
}

cat(
    "\nPDF bytes: ",
    file.info(output_pdf)$size,
    "\n",
    sep = ""
)
