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
    "Wallbank_Dennis_Ray_plus_restricted_Band_balanced_compact_BAgt099_ribbons.pdf"
)

dr_audit_file <- file.path(
    output_dir,
    "DR_broad_correct_pattern_BAgt099_ribbon_sites.tsv"
)

band_audit_file <- file.path(
    output_dir,
    "Band_broad_correct_pattern_BAgt099_ribbon_sites.tsv"
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
# COMPACT FIGURE STYLING
#
# Comparative BA is deliberately neutral so that colour in the
# lower tracks has a single clear meaning.
# ============================================================

neutral_BA <- "#4B5563"
background_BA <- "#D6D9DE"


recolour_fixed_point_layers <- function(
    plot_object
) {

    for (i in seq_along(plot_object$layers)) {

        layer <- plot_object$layers[[i]]

        if (!inherits(
            layer$geom,
            "GeomPoint"
        )) {
            next
        }

        current_colour <- layer$aes_params$colour

        if (
            !is.null(current_colour) &&
            identical(
                current_colour,
                "#BDBDBD"
            )
        ) {
            plot_object$layers[[i]]$aes_params$colour <- (
                background_BA
            )
        }

        if (
            !is.null(current_colour) &&
            identical(
                current_colour,
                "#D55E00"
            )
        ) {
            plot_object$layers[[i]]$aes_params$colour <- (
                neutral_BA
            )
        }
    }

    plot_object
}


dr_ba_plot <- recolour_fixed_point_layers(
    dr_ba_plot
)

band_ba_plot <- recolour_fixed_point_layers(
    band_ba_plot
)


# D/R score-eligible points are mapped to module; make both
# neutral because Dennis versus Ray is already identified by
# the panel/facet structure.

dr_ba_plot <- dr_ba_plot +

    scale_colour_manual(
        values = c(
            Dennis = neutral_BA,
            Ray = neutral_BA
        ),
        guide = "none"
    ) +

    labs(
        title = "Dennis and Ray comparative association",
        subtitle = paste(
            "Hmr18:1.060–1.240 Mb;",
            "TN_Left and TN_Right excluded from comparative scoring."
        )
    ) +

    theme(
        legend.position = "none",
        plot.title = element_text(
            size = 11,
            face = "bold"
        ),
        plot.subtitle = element_text(
            size = 8.2,
            colour = "#555555",
            margin = margin(
                b = 4
            )
        ),
        plot.margin = margin(
            5,
            12,
            0,
            8
        )
    )


band_ba_plot <- band_ba_plot +

    labs(
        title = "Band comparative association",
        subtitle = paste(
            "Hmr18:1.060–1.240 Mb;",
            "TN_Left and TN_Right excluded from comparative scoring."
        )
    ) +

    theme(
        legend.position = "none",
        plot.title = element_text(
            size = 11,
            face = "bold"
        ),
        plot.subtitle = element_text(
            size = 8.2,
            colour = "#555555",
            margin = margin(
                b = 4
            )
        ),
        plot.margin = margin(
            5,
            12,
            0,
            8
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
# First retain the correct focal inheritance direction across
# the whole broad region:
#
#     TN_Left  = 1/1 = biological P/P
#     TN_Right = 0/1 = biological P/A
#
# Then annotate each site with BOTH Dennis and Ray comparative
# BA.  A site is highlighted when either BA is strictly >0.99.
# Correct-direction sites below that threshold remain grey.
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


# ------------------------------------------------------------
# Add Dennis and Ray BA independently.
# ------------------------------------------------------------

if (
    !exists(
        "dennis",
        envir = dr_env,
        inherits = FALSE
    ) ||
    !exists(
        "ray",
        envir = dr_env,
        inherits = FALSE
    )
) {
    stop(
        "D/R environment does not contain Dennis and Ray BA tables."
    )
}

dennis_ba <- dr_env$dennis
ray_ba <- dr_env$ray

required_ba_columns <- c(
    "chrom",
    "pos",
    "balanced_accuracy"
)

for (
    object_name
    in c(
        "dennis_ba",
        "ray_ba"
    )
) {

    object <- get(
        object_name
    )

    missing_columns <- setdiff(
        required_ba_columns,
        names(object)
    )

    if (length(missing_columns) > 0) {
        stop(
            paste(
                object_name,
                "is missing:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }

    object$pos <- suppressWarnings(
        as.numeric(
            object$pos
        )
    )

    object$balanced_accuracy <- suppressWarnings(
        as.numeric(
            object$balanced_accuracy
        )
    )

    object$position_key <- paste(
        object$chrom,
        object$pos,
        sep = ":"
    )

    if (anyDuplicated(
        object$position_key
    )) {
        stop(
            paste(
                object_name,
                "contains duplicated chromosome-position keys;",
                "refusing to guess which BA belongs to a site."
            )
        )
    }

    assign(
        object_name,
        object
    )
}

dr_correct$position_key <- paste(
    dr_correct$chrom,
    dr_correct$pos,
    sep = ":"
)

dennis_index <- match(
    dr_correct$position_key,
    dennis_ba$position_key
)

ray_index <- match(
    dr_correct$position_key,
    ray_ba$position_key
)

dr_correct$Dennis_BA <- (
    dennis_ba$balanced_accuracy[
        dennis_index
    ]
)

dr_correct$Ray_BA <- (
    ray_ba$balanced_accuracy[
        ray_index
    ]
)


# Strict threshold: > 0.99, not >= 0.99.
dr_correct$BA_gt_099 <- (
    (
        is.finite(
            dr_correct$Dennis_BA
        ) &
        dr_correct$Dennis_BA > 0.99
    ) |
    (
        is.finite(
            dr_correct$Ray_BA
        ) &
        dr_correct$Ray_BA > 0.99
    )
)

dennis_for_max <- ifelse(
    is.finite(
        dr_correct$Dennis_BA
    ),
    dr_correct$Dennis_BA,
    -Inf
)

ray_for_max <- ifelse(
    is.finite(
        dr_correct$Ray_BA
    ),
    dr_correct$Ray_BA,
    -Inf
)

dr_correct$maximum_DR_BA <- pmax(
    dennis_for_max,
    ray_for_max
)

dr_correct$maximum_DR_BA[
    !is.finite(
        dr_correct$maximum_DR_BA
    )
] <- NA_real_


dr_long <- rbind(
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "P/P",
        BA_gt_099 = dr_correct$BA_gt_099,
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        BA_gt_099 = dr_correct$BA_gt_099,
        stringsAsFactors = FALSE
    )
)


# BAND RIBBON
#
# Correct focal inheritance direction:
#
#     TN_Left  = biological A/A
#     TN_Right = biological P/A
#
# All such sites remain visible.  Those with Band BA >0.99
# are highlighted; the remainder are grey.
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
    "balanced_accuracy",
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

band_all$balanced_accuracy <- suppressWarnings(
    as.numeric(
        band_all$balanced_accuracy
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

# Strict threshold: >0.99.
band_correct$BA_gt_099 <- (
    is.finite(
        band_correct$balanced_accuracy
    ) &
    band_correct$balanced_accuracy > 0.99
)

band_long <- rbind(
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "A/A",
        BA_gt_099 = band_correct$BA_gt_099,
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        BA_gt_099 = band_correct$BA_gt_099,
        stringsAsFactors = FALSE
    )
)

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
# Compact presentation:
#
#   * only correct-direction focal sites enter the ribbon
#   * BA > 0.99 sites use genotype colour
#   * remaining correct-direction sites are pale grey
#   * expected genotype is written directly into each row label
#   * no separate genotype legend is needed
# ============================================================

genotype_colours <- c(
    "A/A" = "#0072B2",
    "P/A" = "#E69F00",
    "P/P" = "#CC79A7"
)

low_BA_colour <- "#D6D9DE"
connector_colour <- "#667085"


make_correct_pattern_ribbon <- function(
    wide_data,
    long_data,
    left_label,
    right_label,
    summary_text
) {

    grey_wide <- wide_data[
        !wide_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    high_wide <- wide_data[
        wide_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    grey_long <- long_data[
        !long_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    high_long <- long_data[
        long_data$BA_gt_099,
        ,
        drop = FALSE
    ]


    ggplot() +

        # ---------------------------------------------
        # Correct-direction sites with weaker BA
        # ---------------------------------------------

        geom_segment(
            data = grey_wide,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = low_BA_colour,
            alpha = 0.34,
            linewidth = 0.24
        ) +

        geom_point(
            data = grey_long,
            aes(
                x = position_Mb,
                y = y
            ),
            colour = low_BA_colour,
            shape = 124,
            size = 1.65,
            alpha = 0.60
        ) +


        # ---------------------------------------------
        # Correct-direction sites with BA > 0.99
        # ---------------------------------------------

        geom_segment(
            data = high_wide,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = connector_colour,
            alpha = 0.58,
            linewidth = 0.32
        ) +

        geom_point(
            data = high_long,
            aes(
                x = position_Mb,
                y = y,
                colour = biological_GT
            ),
            shape = 124,
            size = 2.15,
            alpha = 0.98
        ) +


        # ---------------------------------------------
        # Scales
        # ---------------------------------------------

        scale_colour_manual(
            values = genotype_colours,
            guide = "none"
        ) +

        scale_x_continuous(
            name = "Position on Hmr18 (Mb)",
            limits = x_limits,
            breaks = seq(
                1.06,
                1.24,
                by = 0.02
            ),
            labels = sprintf(
                "%.2f",
                seq(
                    1.06,
                    1.24,
                    by = 0.02
                )
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
                right_label,
                left_label
            ),
            limits = c(
                -0.20,
                1.20
            ),
            expand = c(
                0,
                0
            )
        ) +


        # ---------------------------------------------
        # One concise explanatory line
        # ---------------------------------------------

        labs(
            title = NULL,
            subtitle = summary_text
        ) +

        theme_minimal(
            base_size = 9
        ) +

        theme(
            plot.subtitle = element_text(
                size = 7.8,
                colour = "#555555",
                margin = margin(
                    b = 2
                )
            ),

            axis.text.y = element_text(
                size = 8.3,
                face = "bold",
                colour = "#333333"
            ),

            axis.text.x = element_text(
                size = 7.6,
                colour = "#555555"
            ),

            axis.title.x = element_text(
                size = 8.2,
                margin = margin(
                    t = 3
                )
            ),

            axis.ticks.y = element_blank(),

            panel.grid.minor = element_blank(),

            panel.grid.major.y = element_blank(),

            panel.grid.major.x = element_line(
                colour = "#EEEEEE",
                linewidth = 0.25
            ),

            legend.position = "none",

            plot.margin = margin(
                1,
                12,
                3,
                8
            )
        )
}


dr_ribbon <- make_correct_pattern_ribbon(
    wide_data = dr_correct,
    long_data = dr_long,

    left_label = "TN_Left   P/P",
    right_label = "TN_Right   P/A",

    summary_text = paste0(
        "Correct-direction sites: ",
        nrow(dr_correct),
        "; coloured = Dennis or Ray BA > 0.99 (n=",
        sum(
            dr_correct$BA_gt_099
        ),
        "); grey = remaining correct-direction sites."
    )
)


band_ribbon <- make_correct_pattern_ribbon(
    wide_data = band_correct,
    long_data = band_long,

    left_label = "TN_Left   A/A",
    right_label = "TN_Right   P/A",

    summary_text = paste0(
        "Correct-direction sites: ",
        nrow(band_correct),
        "; coloured = Band BA > 0.99 (n=",
        sum(
            band_correct$BA_gt_099
        ),
        "); grey = remaining correct-direction sites."
    )
)

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
    width = 10.8,
    height = 11.6,
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
                    7.2,
                    4.4
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
    upper_fraction = 0.80
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
    upper_fraction = 0.66
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
    "ANALYSIS 148 FIGURE COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "D/R correct-direction SNPs: ",
    nrow(dr_correct),
    "\n",
    "D/R correct-direction with Dennis or Ray BA >0.99: ",
    sum(dr_correct$BA_gt_099),
    "\n\n",
    sep = ""
)

cat(
    "Band correct-direction SNPs: ",
    nrow(band_correct),
    "\n",
    "Band correct-direction with BA >0.99: ",
    sum(band_correct$BA_gt_099),
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
