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
    "Wallbank_Dennis_Ray_plus_restricted_Band_unified_balanced_compact_BAgt099_ribbons.pdf"
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
# Shared upper-panel aesthetic language for both A and B.
# Comparative BA is neutral and consistent; colour is reserved
# mainly for the lower genotype ribbons.
# ============================================================

neutral_BA <- "#4C566A"
background_BA <- "#D8DEE9"
shared_text <- "#4A4A4A"
shared_grid <- "#EBEBEB"
shared_strip_fill <- "#F4F4F4"


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


shared_ba_theme <- theme(
    legend.position = "none",

    plot.title = element_text(
        size = 11,
        face = "bold"
    ),

    plot.subtitle = element_text(
        size = 8.2,
        colour = shared_text,
        margin = margin(
            b = 4
        )
    ),

    axis.title.y = element_text(
        size = 8.8,
        margin = margin(
            r = 8
        )
    ),

    axis.text.y = element_text(
        size = 7.8,
        colour = shared_text
    ),

    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    panel.grid.minor = element_blank(),

    panel.grid.major.x = element_line(
        colour = shared_grid,
        linewidth = 0.30
    ),

    panel.grid.major.y = element_line(
        colour = shared_grid,
        linewidth = 0.30
    ),

    strip.text.y.right = element_text(
        size = 8.3,
        face = "bold",
        colour = shared_text,
        angle = 270
    ),

    strip.background = element_rect(
        fill = shared_strip_fill,
        colour = NA
    ),

    plot.margin = margin(
        5,
        12,
        0,
        8
    )
)


dr_ba_plot <- recolour_fixed_point_layers(
    dr_ba_plot
)

band_ba_plot <- recolour_fixed_point_layers(
    band_ba_plot
)


# D/R: module-mapped eligible points should be visually neutral
# because Dennis and Ray are already distinguished by the two
# upper BA tracks.
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
        ),
        y = "Comparative balanced accuracy"
    ) +

    shared_ba_theme


band_ba_plot <- band_ba_plot +

    labs(
        title = "Band comparative association",
        subtitle = paste(
            "Hmr18:1.060–1.240 Mb;",
            "TN_Left and TN_Right excluded from comparative scoring."
        ),
        y = "Comparative balanced accuracy"
    ) +

    shared_ba_theme

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

# ============================================================
# ANALYSIS 158
# SCIENCE-STYLE CONDENSATION TRIALS
#
# Trial 1 = compact quantitative BA tracks
# Trial 2 = BA-height rug tracks
# Trial 4 = high-BA locus schematic
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})


# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------

science_output_dir <- file.path(
    project,
    "analysis",
    "figures",
    "final_B_D_R_combined",
    "science_condensation_trials"
)

dir.create(
    science_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# COMMON LOCUS GEOMETRY
# ------------------------------------------------------------

science_x_limits <- c(
    1.060,
    1.240
)

science_x_breaks <- seq(
    1.06,
    1.24,
    by = 0.02
)

optix_start_science <- 1.088602
optix_end_science <- 1.097846
optix_mid_science <- mean(
    c(
        optix_start_science,
        optix_end_science
    )
)


# ------------------------------------------------------------
# COLOR LANGUAGE
#
# Color is reserved for the three biological hypotheses.
# Everything else remains neutral.
# ------------------------------------------------------------

trait_colours <- c(
    B = "#0072B2",
    D = "#E69F00",
    R = "#CC79A7"
)

neutral_point <- "#C8CED6"
neutral_dark <- "#66707A"
neutral_line <- "#AEB6BF"
neutral_focal <- "#CDD2D8"
strong_focal <- "#303841"
threshold_colour <- "#666666"
optix_colour <- "#3A7CA5"
text_colour <- "#242424"


# ============================================================
# HELPERS
# ============================================================

clamp01 <- function(x) {

    pmin(
        pmax(
            x,
            0
        ),
        1
    )
}


normalise_track <- function(
    x,
    trait,
    ba_column = "balanced_accuracy"
) {

    if (!ba_column %in% names(x)) {
        stop(
            paste(
                trait,
                "BA source is missing column",
                ba_column
            )
        )
    }

    if ("position_Mb" %in% names(x)) {

        position_Mb <- suppressWarnings(
            as.numeric(
                x$position_Mb
            )
        )

    } else if ("pos" %in% names(x)) {

        position_Mb <- suppressWarnings(
            as.numeric(
                x$pos
            )
        ) / 1e6

    } else {

        stop(
            paste(
                trait,
                "BA source has neither position_Mb nor pos."
            )
        )
    }

    BA <- suppressWarnings(
        as.numeric(
            x[[ba_column]]
        )
    )

    result <- data.frame(
        trait = trait,
        position_Mb = position_Mb,
        BA = BA,
        stringsAsFactors = FALSE
    )

    result$pos <- as.integer(
        round(
            result$position_Mb * 1e6
        )
    )

    result <- result[
        is.finite(
            result$position_Mb
        ) &
        is.finite(
            result$BA
        ) &
        result$position_Mb >=
            science_x_limits[1] &
        result$position_Mb <=
            science_x_limits[2],
        ,
        drop = FALSE
    ]

    result
}


science_theme <- theme_classic(
    base_size = 8
) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 9.2,
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
            size = 7.4,
            margin = margin(
                t = 3
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
            linewidth = 0.3
        ),

        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),

        panel.grid = element_blank(),

        plot.margin = margin(
            8,
            33,
            7,
            31
        ),

        legend.position = "none"
    )


# ============================================================
# EXACT COMPARATIVE BA TRACKS USED BY THE ESTABLISHED FIGURE
# ============================================================


# ------------------------------------------------------------
# Dennis and Ray:
# preferentially reuse the final D/R plot_data table because it
# is exactly what fed the established comparative figure.
# ------------------------------------------------------------

if (
    exists(
        "plot_data",
        envir = dr_env,
        inherits = FALSE
    )
) {

    dr_plot_data_science <- dr_env$plot_data

    required_dr_plot_columns <- c(
        "panel",
        "position_Mb",
        "plot_value"
    )

    missing_dr_plot_columns <- setdiff(
        required_dr_plot_columns,
        names(
            dr_plot_data_science
        )
    )

    if (length(missing_dr_plot_columns) > 0L) {
        stop(
            paste(
                "dr_env$plot_data is missing:",
                paste(
                    missing_dr_plot_columns,
                    collapse = ", "
                )
            )
        )
    }

    dennis_track <- data.frame(
        trait = "D",
        position_Mb = suppressWarnings(
            as.numeric(
                dr_plot_data_science$position_Mb[
                    as.character(
                        dr_plot_data_science$panel
                    ) == "Dennis comparative BA"
                ]
            )
        ),
        BA = suppressWarnings(
            as.numeric(
                dr_plot_data_science$plot_value[
                    as.character(
                        dr_plot_data_science$panel
                    ) == "Dennis comparative BA"
                ]
            )
        ),
        stringsAsFactors = FALSE
    )

    ray_track <- data.frame(
        trait = "R",
        position_Mb = suppressWarnings(
            as.numeric(
                dr_plot_data_science$position_Mb[
                    as.character(
                        dr_plot_data_science$panel
                    ) == "Ray comparative BA"
                ]
            )
        ),
        BA = suppressWarnings(
            as.numeric(
                dr_plot_data_science$plot_value[
                    as.character(
                        dr_plot_data_science$panel
                    ) == "Ray comparative BA"
                ]
            )
        ),
        stringsAsFactors = FALSE
    )

    for (nm in c("dennis_track", "ray_track")) {

        z <- get(nm)

        z$pos <- as.integer(
            round(
                z$position_Mb * 1e6
            )
        )

        z <- z[
            is.finite(z$position_Mb) &
            is.finite(z$BA) &
            z$position_Mb >= science_x_limits[1] &
            z$position_Mb <= science_x_limits[2],
            ,
            drop = FALSE
        ]

        assign(
            nm,
            z
        )
    }

} else {

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
            "Could not recover the established Dennis/Ray BA data."
        )
    }

    dennis_track <- normalise_track(
        dr_env$dennis,
        "D"
    )

    ray_track <- normalise_track(
        dr_env$ray,
        "R"
    )
}


# ------------------------------------------------------------
# Band:
# reuse the established finite BA table if available.
# ------------------------------------------------------------

if (
    exists(
        "ba_data",
        envir = band_env,
        inherits = FALSE
    )
) {

    band_source_science <- band_env$ba_data

} else if (
    exists(
        "d",
        envir = band_env,
        inherits = FALSE
    )
) {

    band_source_science <- band_env$d

} else {

    stop(
        "Could not recover the established Band BA table."
    )
}

band_track <- normalise_track(
    band_source_science,
    "B"
)


all_ba_science <- rbind(
    band_track,
    dennis_track,
    ray_track
)

all_ba_science$trait <- factor(
    all_ba_science$trait,
    levels = c(
        "B",
        "D",
        "R"
    ),
    ordered = TRUE
)


# ============================================================
# HIGH-BA SITES WITH THE CORRECT FOCAL L->R DIRECTION
# ============================================================

required_dr_correct <- c(
    "position_Mb",
    "Dennis_BA",
    "Ray_BA",
    "BA_gt_099"
)

if (
    length(
        setdiff(
            required_dr_correct,
            names(dr_correct)
        )
    ) > 0L
) {
    stop(
        "dr_correct is missing the Analysis-149 BA columns."
    )
}

required_band_correct <- c(
    "position_Mb",
    "balanced_accuracy",
    "BA_gt_099"
)

if (
    length(
        setdiff(
            required_band_correct,
            names(band_correct)
        )
    ) > 0L
) {
    stop(
        "band_correct is missing the Analysis-149 BA columns."
    )
}


high_B <- data.frame(
    trait = "B",
    position_Mb = band_correct$position_Mb[
        band_correct$BA_gt_099 %in% TRUE
    ],
    BA = band_correct$balanced_accuracy[
        band_correct$BA_gt_099 %in% TRUE
    ],
    stringsAsFactors = FALSE
)

high_D <- data.frame(
    trait = "D",
    position_Mb = dr_correct$position_Mb[
        is.finite(dr_correct$Dennis_BA) &
        dr_correct$Dennis_BA > 0.99
    ],
    BA = dr_correct$Dennis_BA[
        is.finite(dr_correct$Dennis_BA) &
        dr_correct$Dennis_BA > 0.99
    ],
    stringsAsFactors = FALSE
)

high_R <- data.frame(
    trait = "R",
    position_Mb = dr_correct$position_Mb[
        is.finite(dr_correct$Ray_BA) &
        dr_correct$Ray_BA > 0.99
    ],
    BA = dr_correct$Ray_BA[
        is.finite(dr_correct$Ray_BA) &
        dr_correct$Ray_BA > 0.99
    ],
    stringsAsFactors = FALSE
)

high_correct <- rbind(
    high_B,
    high_D,
    high_R
)

high_correct$pos <- as.integer(
    round(
        high_correct$position_Mb * 1e6
    )
)

high_correct$trait <- factor(
    high_correct$trait,
    levels = c(
        "B",
        "D",
        "R"
    ),
    ordered = TRUE
)


high_keys <- paste(
    as.character(
        high_correct$trait
    ),
    high_correct$pos,
    sep = ":"
)

all_keys <- paste(
    as.character(
        all_ba_science$trait
    ),
    all_ba_science$pos,
    sep = ":"
)

all_ba_science$high_correct_gt099 <- (
    all_keys %in% high_keys
)


# ============================================================
# FOCAL-COMPATIBLE SITE TRACKS
# ============================================================

focal_B <- data.frame(
    track = "B  A/A -> P/A",
    position_Mb = band_correct$position_Mb,
    high = band_correct$BA_gt_099 %in% TRUE,
    stringsAsFactors = FALSE
)

focal_DR <- data.frame(
    track = "D/R  P/P -> P/A",
    position_Mb = dr_correct$position_Mb,
    high = dr_correct$BA_gt_099 %in% TRUE,
    stringsAsFactors = FALSE
)


# ============================================================
# WRITE A SMALL AUDIT FOR THE THREE TRIALS
# ============================================================

write.table(
    high_correct[
        order(
            high_correct$trait,
            high_correct$position_Mb
        ),
        ,
        drop = FALSE
    ],
    file = file.path(
        science_output_dir,
        "Science_trials_high_correct_BDR_BAgt099_sites.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

focal_count_audit <- data.frame(
    category = c(
        "Band focal-compatible",
        "D/R focal-compatible",
        "Band focal-compatible + BA>0.99",
        "D/R focal-compatible + Dennis or Ray BA>0.99",
        "Dennis BA>0.99 within D/R-compatible",
        "Ray BA>0.99 within D/R-compatible"
    ),
    n = c(
        nrow(band_correct),
        nrow(dr_correct),
        sum(band_correct$BA_gt_099 %in% TRUE),
        sum(dr_correct$BA_gt_099 %in% TRUE),
        nrow(high_D),
        nrow(high_R)
    ),
    stringsAsFactors = FALSE
)

write.table(
    focal_count_audit,
    file = file.path(
        science_output_dir,
        "Science_trials_site_count_audit.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# TRIAL 1
# THREE COMPACT QUANTITATIVE BA TRACKS + SHARED FOCAL TRACK
# ============================================================

ba_bases_1 <- c(
    B = 3.25,
    D = 2.05,
    R = 0.85
)

ba_span_1 <- 0.78

all_ba_science$y_trial1 <- vapply(
    seq_len(
        nrow(all_ba_science)
    ),
    function(i) {

        trait_name <- as.character(
            all_ba_science$trait[i]
        )

        ba_bases_1[[trait_name]] +
            clamp01(
                (all_ba_science$BA[i] - 0.5) / 0.5
            ) * ba_span_1
    },
    numeric(1)
)

high_correct$y_trial1 <- vapply(
    seq_len(
        nrow(high_correct)
    ),
    function(i) {

        trait_name <- as.character(
            high_correct$trait[i]
        )

        ba_bases_1[[trait_name]] +
            clamp01(
                (high_correct$BA[i] - 0.5) / 0.5
            ) * ba_span_1
    },
    numeric(1)
)

track_guides_1 <- do.call(
    rbind,
    lapply(
        names(ba_bases_1),
        function(trait_name) {

            data.frame(
                trait = trait_name,
                baseline = ba_bases_1[[trait_name]],
                threshold = ba_bases_1[[trait_name]] +
                    ((0.8 - 0.5) / 0.5) * ba_span_1,
                top = ba_bases_1[[trait_name]] + ba_span_1,
                stringsAsFactors = FALSE
            )
        }
    )
)

trait_labels_1 <- data.frame(
    x = science_x_limits[1],
    y = unname(ba_bases_1) + ba_span_1 / 2,
    label = names(ba_bases_1),
    stringsAsFactors = FALSE
)

focal_B_1 <- focal_B
focal_B_1$y <- -0.03

focal_DR_1 <- focal_DR
focal_DR_1$y <- -0.38

trial1 <- ggplot() +

    geom_segment(
        data = track_guides_1,
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = baseline,
            yend = baseline
        ),
        inherit.aes = FALSE,
        colour = neutral_line,
        linewidth = 0.25
    ) +

    geom_segment(
        data = track_guides_1,
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = threshold,
            yend = threshold
        ),
        inherit.aes = FALSE,
        colour = threshold_colour,
        linewidth = 0.30,
        linetype = "22"
    ) +

    geom_point(
        data = all_ba_science,
        aes(
            x = position_Mb,
            y = y_trial1
        ),
        colour = neutral_point,
        alpha = 0.52,
        size = 0.42,
        stroke = 0
    ) +

    geom_point(
        data = high_correct,
        aes(
            x = position_Mb,
            y = y_trial1,
            colour = trait
        ),
        alpha = 0.98,
        size = 0.92,
        stroke = 0
    ) +

    geom_segment(
        data = focal_B_1,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.075,
            yend = y + 0.075
        ),
        colour = neutral_focal,
        linewidth = 0.30
    ) +

    geom_segment(
        data = focal_B_1[
            focal_B_1$high,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.085,
            yend = y + 0.085
        ),
        colour = strong_focal,
        linewidth = 0.55
    ) +

    geom_segment(
        data = focal_DR_1,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.075,
            yend = y + 0.075
        ),
        colour = neutral_focal,
        linewidth = 0.30
    ) +

    geom_segment(
        data = focal_DR_1[
            focal_DR_1$high,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.085,
            yend = y + 0.085
        ),
        colour = strong_focal,
        linewidth = 0.55
    ) +

    geom_segment(
        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 4.36,
            yend = 4.36
        ),
        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 3.0,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = optix_mid_science,
        y = 4.48,
        label = "optix",
        colour = optix_colour,
        size = 2.5,
        fontface = "italic"
    ) +

    geom_text(
        data = trait_labels_1,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 1.65,
        fontface = "bold",
        size = 3.0,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.03,
        label = "B   A/A -> P/A",
        hjust = 1.12,
        size = 2.35,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.38,
        label = "D/R   P/P -> P/A",
        hjust = 1.12,
        size = 2.35,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1] + 0.001,
        y = 4.10,
        label = "BA 0.5-1.0; dashed = 0.8",
        hjust = 0,
        size = 2.15,
        colour = "#555555"
    ) +

    scale_colour_manual(
        values = trait_colours,
        drop = FALSE
    ) +

    scale_x_continuous(
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
            -0.62,
            4.55
        ),
        clip = "off"
    ) +

    labs(
        title = "Comparative association and focal L-R genotype context",
        subtitle = paste(
            "B, Dennis (D), and Ray (R) share one Hmr18 axis;",
            "colored points are focal-compatible sites with BA > 0.99."
        ),
        x = "Position on Hmr18 (Mb)",
        caption = paste0(
            "Focal-compatible sites: B n=",
            nrow(band_correct),
            "; D/R n=",
            nrow(dr_correct),
            ".  BA>0.99 focal-compatible: B n=",
            sum(band_correct$BA_gt_099 %in% TRUE),
            "; D/R union n=",
            sum(dr_correct$BA_gt_099 %in% TRUE),
            "."
        )
    ) +

    science_theme


# ============================================================
# TRIAL 2
# BA-HEIGHT RUG TRACKS + SHARED FOCAL TRACK
# ============================================================

rug_bases <- c(
    B = 2.55,
    D = 1.55,
    R = 0.55
)

all_ba_science$tick_height <- ifelse(
    all_ba_science$BA > 0.99,
    0.55,
    ifelse(
        all_ba_science$BA >= 0.90,
        0.42,
        ifelse(
            all_ba_science$BA >= 0.80,
            0.27,
            0.11
        )
    )
)

all_ba_science$rug_base <- unname(
    rug_bases[
        as.character(
            all_ba_science$trait
        )
    ]
)

all_ba_science$rug_top <- (
    all_ba_science$rug_base +
    all_ba_science$tick_height
)

high_rug <- high_correct
high_rug$rug_base <- unname(
    rug_bases[
        as.character(
            high_rug$trait
        )
    ]
)
high_rug$rug_top <- high_rug$rug_base + 0.58

focal_B_2 <- focal_B
focal_B_2$y <- -0.12

focal_DR_2 <- focal_DR
focal_DR_2$y <- -0.42

rug_labels <- data.frame(
    x = science_x_limits[1],
    y = unname(rug_bases) + 0.27,
    label = names(rug_bases),
    stringsAsFactors = FALSE
)

trial2 <- ggplot() +

    geom_segment(
        data = data.frame(
            baseline = unname(rug_bases)
        ),
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = baseline,
            yend = baseline
        ),
        inherit.aes = FALSE,
        colour = neutral_line,
        linewidth = 0.24
    ) +

    geom_segment(
        data = all_ba_science,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = rug_base,
            yend = rug_top
        ),
        colour = neutral_dark,
        alpha = 0.34,
        linewidth = 0.24
    ) +

    geom_segment(
        data = high_rug,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = rug_base,
            yend = rug_top,
            colour = trait
        ),
        alpha = 1,
        linewidth = 0.66
    ) +

    geom_segment(
        data = focal_B_2,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.065,
            yend = y + 0.065
        ),
        colour = neutral_focal,
        linewidth = 0.27
    ) +

    geom_segment(
        data = focal_B_2[
            focal_B_2$high,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.075,
            yend = y + 0.075
        ),
        colour = strong_focal,
        linewidth = 0.52
    ) +

    geom_segment(
        data = focal_DR_2,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.065,
            yend = y + 0.065
        ),
        colour = neutral_focal,
        linewidth = 0.27
    ) +

    geom_segment(
        data = focal_DR_2[
            focal_DR_2$high,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y - 0.075,
            yend = y + 0.075
        ),
        colour = strong_focal,
        linewidth = 0.52
    ) +

    geom_segment(
        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 3.42,
            yend = 3.42
        ),
        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 3.0,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = optix_mid_science,
        y = 3.54,
        label = "optix",
        colour = optix_colour,
        size = 2.5,
        fontface = "italic"
    ) +

    geom_text(
        data = rug_labels,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 1.65,
        fontface = "bold",
        size = 3.0,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.12,
        label = "B   A/A -> P/A",
        hjust = 1.12,
        size = 2.35,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.42,
        label = "D/R   P/P -> P/A",
        hjust = 1.12,
        size = 2.35,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[2],
        y = 3.20,
        label = "tick height: <0.8 | 0.8-0.9 | 0.9-0.99 | >0.99",
        hjust = 1,
        size = 2.05,
        colour = "#555555"
    ) +

    scale_colour_manual(
        values = trait_colours,
        drop = FALSE
    ) +

    scale_x_continuous(
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
            -0.62,
            3.61
        ),
        clip = "off"
    ) +

    labs(
        title = "Comparative association as compact BA rugs",
        subtitle = paste(
            "Tick height encodes BA class;",
            "color marks focal-compatible sites with BA > 0.99."
        ),
        x = "Position on Hmr18 (Mb)"
    ) +

    science_theme


# ============================================================
# TRIAL 4
# AGGRESSIVE HIGH-BA LOCUS SCHEMATIC
#
# No BA cloud and no arbitrary window/bin definition.
# We retain exact positions of focal-compatible BA>0.99 sites.
# ============================================================

schematic_bases <- c(
    B = 2.30,
    D = 1.40,
    R = 0.50
)

high_schematic <- high_correct
high_schematic$base <- unname(
    schematic_bases[
        as.character(
            high_schematic$trait
        )
    ]
)

schematic_labels <- data.frame(
    x = science_x_limits[1],
    y = unname(schematic_bases),
    label = names(schematic_bases),
    stringsAsFactors = FALSE
)

schematic_counts <- data.frame(
    trait = c(
        "B",
        "D",
        "R"
    ),
    x = science_x_limits[2],
    y = unname(schematic_bases),
    label = paste0(
        "n=",
        c(
            nrow(high_B),
            nrow(high_D),
            nrow(high_R)
        )
    ),
    stringsAsFactors = FALSE
)

trial4 <- ggplot() +

    geom_segment(
        data = data.frame(
            baseline = unname(schematic_bases)
        ),
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = baseline,
            yend = baseline
        ),
        inherit.aes = FALSE,
        colour = "#D5D9DE",
        linewidth = 0.55
    ) +

    geom_segment(
        data = high_schematic,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = base - 0.19,
            yend = base + 0.19,
            colour = trait
        ),
        linewidth = 0.82,
        lineend = "butt"
    ) +

    geom_segment(
        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 3.02,
            yend = 3.02
        ),
        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 4.0,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = optix_mid_science,
        y = 3.18,
        label = "optix",
        colour = optix_colour,
        size = 2.8,
        fontface = "italic"
    ) +

    geom_text(
        data = schematic_labels,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 1.65,
        fontface = "bold",
        size = 3.25,
        colour = text_colour
    ) +

    geom_text(
        data = schematic_counts,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = -0.18,
        size = 2.35,
        colour = "#444444"
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.02,
        label = paste0(
            "Focal L->R:   B  A/A->P/A     D/R  P/P->P/A"
        ),
        hjust = 0,
        size = 2.45,
        colour = text_colour
    ) +

    annotate(
        "text",
        x = science_x_limits[1],
        y = -0.34,
        label = paste0(
            nrow(band_correct),
            " B-compatible and ",
            nrow(dr_correct),
            " D/R-compatible SNPs; ",
            sum(band_correct$BA_gt_099 %in% TRUE),
            " B and ",
            sum(dr_correct$BA_gt_099 %in% TRUE),
            " unique D/R sites also exceed BA 0.99."
        ),
        hjust = 0,
        size = 2.15,
        colour = "#555555"
    ) +

    scale_colour_manual(
        values = trait_colours,
        drop = FALSE
    ) +

    scale_x_continuous(
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
            -0.50,
            3.24
        ),
        clip = "off"
    ) +

    labs(
        title = "Focal-compatible high-accuracy association at optix",
        subtitle = "Ticks show exact positions of B, Dennis, or Ray BA > 0.99 sites.",
        x = "Position on Hmr18 (Mb)"
    ) +

    science_theme


# ============================================================
# SAVE ALL THREE OPTIONS
# ============================================================

trial1_pdf <- file.path(
    science_output_dir,
    "Science_trial_1_three_compact_BDR_BA_tracks.pdf"
)

trial1_png <- file.path(
    science_output_dir,
    "Science_trial_1_three_compact_BDR_BA_tracks.png"
)

trial2_pdf <- file.path(
    science_output_dir,
    "Science_trial_2_BDR_BA_rug_tracks.pdf"
)

trial2_png <- file.path(
    science_output_dir,
    "Science_trial_2_BDR_BA_rug_tracks.png"
)

trial4_pdf <- file.path(
    science_output_dir,
    "Science_trial_4_high_BA_locus_schematic.pdf"
)

trial4_png <- file.path(
    science_output_dir,
    "Science_trial_4_high_BA_locus_schematic.png"
)

ggsave(
    filename = trial1_pdf,
    plot = trial1,
    width = 7.2,
    height = 4.35,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = trial1_png,
    plot = trial1,
    width = 7.2,
    height = 4.35,
    units = "in",
    dpi = 450,
    limitsize = FALSE
)

ggsave(
    filename = trial2_pdf,
    plot = trial2,
    width = 7.2,
    height = 3.55,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = trial2_png,
    plot = trial2,
    width = 7.2,
    height = 3.55,
    units = "in",
    dpi = 450,
    limitsize = FALSE
)

ggsave(
    filename = trial4_pdf,
    plot = trial4,
    width = 7.2,
    height = 2.85,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = trial4_png,
    plot = trial4,
    width = 7.2,
    height = 2.85,
    units = "in",
    dpi = 450,
    limitsize = FALSE
)


# ============================================================
# MULTIPAGE COMPARISON PDF
# ============================================================

comparison_pdf <- file.path(
    science_output_dir,
    "Science_trials_1_2_4_comparison.pdf"
)

grDevices::cairo_pdf(
    filename = comparison_pdf,
    width = 7.2,
    height = 4.5,
    onefile = TRUE
)

print(trial1)
print(trial2)
print(trial4)

grDevices::dev.off()


# ============================================================
# TERMINAL SUMMARY
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 158 SCIENCE CONDENSATION TRIALS COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    focal_count_audit,
    row.names = FALSE
)

cat(
    "\nB high sites: ",
    nrow(high_B),
    "\nD high sites: ",
    nrow(high_D),
    "\nR high sites: ",
    nrow(high_R),
    "\nD/R high-site union: ",
    sum(dr_correct$BA_gt_099 %in% TRUE),
    "\n\n",
    sep = ""
)

cat(
    "Trial 1:\n",
    trial1_pdf,
    "\n\nTrial 2:\n",
    trial2_pdf,
    "\n\nTrial 4:\n",
    trial4_pdf,
    "\n\nComparison:\n",
    comparison_pdf,
    "\n",
    sep = ""
)
