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

checkpoint_dir <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
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


# ============================================================
# INPUTS
# ============================================================

dr_script <- file.path(
    script_dir,
    "115b_plot_D_R_Wallbank_fixed_dominant_adjusted_LR.R"
)

band_script <- file.path(
    script_dir,
    "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R"
)

wide_file <- file.path(
    checkpoint_dir,
    "focal_calls_at_frozen_sites_biological.wide.tsv"
)

sites_file <- file.path(
    checkpoint_dir,
    "frozen_comparative_sites_reference_polarized.tsv"
)

output_pdf <- file.path(
    output_dir,
    "Wallbank_Dennis_Ray_plus_restricted_Band_FROZEN_focal_sites.pdf"
)

audit_file <- file.path(
    output_dir,
    "Wallbank_B_D_R_FROZEN_focal_sites_plotted.tsv"
)

required_files <- c(
    dr_script,
    band_script,
    wide_file,
    sites_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {
    stop(
        paste(
            "Missing required file(s):",
            paste(
                missing_files,
                collapse = "\n"
            )
        )
    )
}


# ============================================================
# LOAD EXISTING PLOT ENVIRONMENT WITHOUT WRITING ITS PDF
#
# This deliberately preserves the already-tested comparative
# BA construction in 115b and 142.  We only replace the focal
# genotype panels below them.
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
                "did not create ba_plot"
            )
        )
    }

    env
}


cat(
    "\nLoading already-established comparative BA plots...\n"
)

dr_env <- load_plot_environment(
    dr_script
)

band_env <- load_plot_environment(
    band_script
)


# ============================================================
# FROZEN SITE TABLES
# ============================================================

wide <- read.delim(
    wide_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sites <- read.delim(
    sites_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_wide_columns <- c(
    "site_id",
    "module",
    "chrom",
    "pos",
    "left_QC_pass_DP2_GQ20",
    "left_biological_genotype",
    "left_expected_genotype",
    "right_QC_pass_DP2_GQ20",
    "right_biological_genotype",
    "right_expected_genotype",
    "both_callable",
    "pair_matches_expectation"
)

missing_wide_columns <- setdiff(
    required_wide_columns,
    names(wide)
)

if (length(missing_wide_columns) > 0) {
    stop(
        paste(
            "Frozen focal table is missing columns:",
            paste(
                missing_wide_columns,
                collapse = ", "
            )
        )
    )
}

required_site_columns <- c(
    "site_id",
    "module",
    "chrom",
    "pos",
    "left_expected_genotype",
    "right_expected_genotype"
)

missing_site_columns <- setdiff(
    required_site_columns,
    names(sites)
)

if (length(missing_site_columns) > 0) {
    stop(
        paste(
            "Frozen site table is missing columns:",
            paste(
                missing_site_columns,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# HARD GUARDS ON THE FROZEN SET
#
# These stop this figure silently drifting back into any of
# the later ad-hoc restricted site sets.
# ============================================================

expected_counts <- c(
    Band = 8L,
    Ray1 = 15L,
    Ray2 = 65L,
    Dennis1 = 3L,
    Dennis2 = 21L,
    Yb = 5L
)

observed_counts <- table(
    sites$module
)

for (module_name in names(expected_counts)) {

    observed_n <- unname(
        observed_counts[
            module_name
        ]
    )

    if (
        length(observed_n) == 0 ||
        is.na(observed_n) ||
        observed_n != expected_counts[
            module_name
        ]
    ) {
        stop(
            paste(
                "Frozen-site count changed for",
                module_name,
                "- expected",
                expected_counts[
                    module_name
                ],
                "but found",
                ifelse(
                    length(observed_n) == 0,
                    0,
                    observed_n
                )
            )
        )
    }
}

if (nrow(sites) != 117) {
    stop(
        paste(
            "Expected exactly 117 frozen comparative sites;",
            "found",
            nrow(sites)
        )
    )
}

if (nrow(wide) != 117) {
    stop(
        paste(
            "Expected exactly 117 focal rows at frozen sites;",
            "found",
            nrow(wide)
        )
    )
}

if (
    anyDuplicated(
        sites$site_id
    )
) {
    stop(
        "Duplicate site_id values in frozen site table."
    )
}

if (
    anyDuplicated(
        wide$site_id
    )
) {
    stop(
        "Duplicate site_id values in frozen focal table."
    )
}

if (
    !setequal(
        sites$site_id,
        wide$site_id
    )
) {
    stop(
        paste(
            "Frozen comparative sites and frozen focal calls",
            "do not contain the same site_id set."
        )
    )
}


# ============================================================
# BOOLEAN NORMALISATION
# ============================================================

as_bool <- function(
    x
) {

    if (is.logical(x)) {
        return(x)
    }

    y <- tolower(
        trimws(
            as.character(x)
        )
    )

    out <- rep(
        NA,
        length(y)
    )

    out[
        y %in% c(
            "true",
            "t",
            "1",
            "yes"
        )
    ] <- TRUE

    out[
        y %in% c(
            "false",
            "f",
            "0",
            "no"
        )
    ] <- FALSE

    out
}


wide$left_QC <- as_bool(
    wide$left_QC_pass_DP2_GQ20
)

wide$right_QC <- as_bool(
    wide$right_QC_pass_DP2_GQ20
)

wide$pair_expected <- as_bool(
    wide$pair_matches_expectation
)

wide$pos <- as.numeric(
    wide$pos
)

wide$position_Mb <- (
    wide$pos /
    1e6
)


# ============================================================
# USE THE BIOLOGICAL GENOTYPE CALLS DIRECTLY
#
# No reciprocal-mixture reclassification.
# No fallback to the 115b/142 focal calls.
# Failed DP2/GQ20 calls remain explicit "No call".
# ============================================================

clean_biological_gt <- function(
    gt,
    qc
) {

    gt <- toupper(
        trimws(
            as.character(gt)
        )
    )

    gt[
        is.na(qc) |
        !qc |
        !(gt %in% c(
            "AA",
            "PA",
            "PP"
        ))
    ] <- "NO CALL"

    out <- rep(
        "No call",
        length(gt)
    )

    out[
        gt == "AA"
    ] <- "A/A"

    out[
        gt == "PA"
    ] <- "P/A"

    out[
        gt == "PP"
    ] <- "P/P"

    out
}


wide$left_plot_GT <- clean_biological_gt(
    wide$left_biological_genotype,
    wide$left_QC
)

wide$right_plot_GT <- clean_biological_gt(
    wide$right_biological_genotype,
    wide$right_QC
)

wide$both_plot_callable <- (
    wide$left_plot_GT != "No call" &
    wide$right_plot_GT != "No call"
)


# ============================================================
# EXACT ANALYSED MODULE INTERVALS
# ============================================================

module_intervals <- data.frame(
    module = c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
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

x_limits <- c(
    1.060,
    1.240
)


# ============================================================
# PLOT DATA
# ============================================================

make_focal_long <- function(
    modules
) {

    d <- wide[
        wide$module %in% modules,
        ,
        drop = FALSE
    ]

    d <- d[
        order(
            d$pos,
            d$module
        ),
        ,
        drop = FALSE
    ]

    long <- rbind(
        data.frame(
            site_id = d$site_id,
            module = d$module,
            position_Mb = d$position_Mb,
            side = "TN_Left",
            y = 1,
            genotype = d$left_plot_GT,
            stringsAsFactors = FALSE
        ),
        data.frame(
            site_id = d$site_id,
            module = d$module,
            position_Mb = d$position_Mb,
            side = "TN_Right",
            y = 0,
            genotype = d$right_plot_GT,
            stringsAsFactors = FALSE
        )
    )

    list(
        wide = d,
        long = long
    )
}


dr_focal_data <- make_focal_long(
    c(
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    )
)

band_focal_data <- make_focal_long(
    "Band"
)


# ============================================================
# GENOTYPE COLOURS
#
# Same colour-blind-safe genotype family used in the earlier
# genotype plots.
# ============================================================

genotype_colours <- c(
    "A/A" = "#E69F00",
    "P/A" = "#CC79A7",
    "P/P" = "#0072B2",
    "No call" = "#BDBDBD"
)


# ============================================================
# BUILD FROZEN-SITE FOCAL STRIP
# ============================================================

make_focal_plot <- function(
    focal_data,
    modules,
    title,
    subtitle
) {

    d <- focal_data$wide
    long <- focal_data$long

    intervals <- module_intervals[
        module_intervals$module
            %in% modules,
        ,
        drop = FALSE
    ]

    expected_segments <- d[
        d$both_plot_callable &
        !is.na(d$pair_expected) &
        d$pair_expected,
        ,
        drop = FALSE
    ]

    other_segments <- d[
        d$both_plot_callable &
        (
            is.na(d$pair_expected) |
            !d$pair_expected
        ),
        ,
        drop = FALSE
    ]

    ggplot() +

        geom_rect(
            data = intervals,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = -Inf,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = "#8C8C8C",
            alpha = 0.10
        ) +

        geom_segment(
            data = other_segments,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = "#9E9E9E",
            linewidth = 0.30,
            alpha = 0.55
        ) +

        geom_segment(
            data = expected_segments,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = "#009E73",
            linewidth = 0.42,
            alpha = 0.75
        ) +

        geom_point(
            data = long,
            aes(
                x = position_Mb,
                y = y,
                colour = genotype
            ),
            shape = 124,
            size = 2.5,
            alpha = 0.92
        ) +

        scale_colour_manual(
            values = genotype_colours,
            breaks = c(
                "A/A",
                "P/A",
                "P/P",
                "No call"
            ),
            drop = FALSE,
            name = "Biological genotype"
        ) +

        scale_x_continuous(
            name = "Position across the broad Optix interval (Mb)",
            limits = x_limits,
            breaks = c(
                1.10,
                1.15,
                1.20
            ),
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
            title = title,
            subtitle = subtitle,
            caption = paste(
                "Green connectors mark callable left-right pairs matching",
                "the pre-specified biological expectation;",
                "gray connectors mark other callable pairs.",
                "Failed or unavailable DP2/GQ20 focal calls remain No call."
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


dr_focal_plot <- make_focal_plot(
    dr_focal_data,
    c(
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    ),
    title = paste(
        "TN_Left and TN_Right genotypes at frozen",
        "Ray and Dennis comparative sites"
    ),
    subtitle = paste(
        "All 104 frozen Hmr18 sites are retained:",
        "Ray 1 = 15, Ray 2 = 65, Dennis 1 = 3, Dennis 2 = 21.",
        "Calls are the direct DP2/GQ20 biological genotypes from the frozen-site checkpoint;",
        "no reciprocal-mixture recoding."
    )
)

band_focal_plot <- make_focal_plot(
    band_focal_data,
    "Band",
    title = paste(
        "TN_Left and TN_Right genotypes at frozen",
        "Band comparative sites"
    ),
    subtitle = paste(
        "All 8 frozen Band sites are retained.",
        "Expected state is TN_Left=A/A and TN_Right=P/A.",
        "Calls are the direct DP2/GQ20 biological genotypes from the frozen-site checkpoint;",
        "no reciprocal-mixture recoding."
    )
)


# ============================================================
# PRESERVE THE COMPARATIVE BA PANELS
#
# We alter only wording that referred specifically to the old
# adjusted focal strip.  The comparative data/layers themselves
# are untouched.
# ============================================================

dr_ba_plot <- dr_env$ba_plot +

    labs(
        title = paste(
            "Dennis and Ray comparative association",
            "across the broad Optix interval"
        ),
        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Triangles retain the established Wallbank Dennis/Ray",
            "comparative highlights.",
            "TN_Left and TN_Right were excluded from comparative SNP selection."
        )
    )

band_ba_plot <- band_env$ba_plot +

    labs(
        subtitle = paste(
            "Hmr18:1.060-1.240 Mb.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right;",
            "canonical frozen-site focal genotypes are displayed below."
        )
    )


# ============================================================
# CONVERT TO GROBS AND ALIGN EACH BA / FOCAL PAIR
# ============================================================

dr_ba_grob <- ggplotGrob(
    dr_ba_plot
)

dr_focal_grob <- ggplotGrob(
    dr_focal_plot
)

dr_widths <- grid::unit.pmax(
    dr_ba_grob$widths,
    dr_focal_grob$widths
)

dr_ba_grob$widths <- dr_widths
dr_focal_grob$widths <- dr_widths


band_ba_grob <- ggplotGrob(
    band_ba_plot
)

band_focal_grob <- ggplotGrob(
    band_focal_plot
)

band_widths <- grid::unit.pmax(
    band_ba_grob$widths,
    band_focal_grob$widths
)

band_ba_grob$widths <- band_widths
band_focal_grob$widths <- band_widths


# ============================================================
# DRAWING HELPERS
# ============================================================

draw_pair <- function(
    ba_grob,
    focal_grob,
    panel_label,
    ba_fraction = 0.66
) {

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 2,
                ncol = 1,
                heights = unit(
                    c(
                        ba_fraction,
                        1 - ba_fraction
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
        focal_grob
    )

    popViewport()

    popViewport()
}


# ============================================================
# WRITE VECTOR PDF
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
    dr_ba_grob,
    dr_focal_grob,
    "A",
    ba_fraction = 0.67
)

popViewport()


pushViewport(
    viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
    )
)

draw_pair(
    band_ba_grob,
    band_focal_grob,
    "B",
    ba_fraction = 0.64
)

popViewport()

popViewport()

dev.off()


# ============================================================
# AUDIT EXACTLY WHAT WAS PLOTTED
# ============================================================

audit <- wide[
    wide$module %in% c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    ),
    c(
        "site_id",
        "module",
        "chrom",
        "pos",
        "left_biological_genotype",
        "left_QC_pass_DP2_GQ20",
        "left_expected_genotype",
        "right_biological_genotype",
        "right_QC_pass_DP2_GQ20",
        "right_expected_genotype",
        "both_callable",
        "pair_matches_expectation"
    ),
    drop = FALSE
]

audit$left_plotted_genotype <- wide$left_plot_GT[
    wide$module %in% c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    )
]

audit$right_plotted_genotype <- wide$right_plot_GT[
    wide$module %in% c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    )
]

write.table(
    audit,
    file = audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# TERMINAL SUMMARY
# ============================================================

cat(
    "\n========================================\n",
    "FROZEN-SITE B / D / R FIGURE REMAKE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "Frozen site counts used:\n"
)

print(
    table(
        audit$module
    )
)

cat(
    "\nFocal pair expectation summary:\n"
)

print(
    with(
        audit,
        table(
            module,
            pair_matches_expectation,
            useNA = "ifany"
        )
    )
)

cat(
    "\nPDF:\n",
    output_pdf,
    "\n\nAudit table:\n",
    audit_file,
    "\n",
    sep = ""
)

if (
    file.exists(output_pdf) &&
    file.info(output_pdf)$size > 0
) {
    cat(
        "\nPDF bytes: ",
        file.info(output_pdf)$size,
        "\n",
        sep = ""
    )
} else {
    stop(
        "Output PDF was not created."
    )
}
