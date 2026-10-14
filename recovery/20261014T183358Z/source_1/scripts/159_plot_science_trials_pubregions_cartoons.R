#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_file <- file.path(
    project,
    "scripts",
    "158_plot_science_condensation_trials.R"
)

if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing working Analysis-158 script:",
            source_file
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

# ============================================================
# SOURCE THE WORKING 158 ANALYSIS
#
# This recreates trial1 / trial2 and all the underlying shared
# objects, but leaves the original 158 script unchanged.
# ============================================================

source(
    source_file,
    local = TRUE
)

required_objects <- c(
    "trial1",
    "trial2",
    "science_x_limits",
    "science_x_breaks",
    "trait_colours",
    "text_colour",
    "optix_colour"
)

missing_objects <- setdiff(
    required_objects,
    ls()
)

if (length(missing_objects) > 0) {
    stop(
        paste(
            "After sourcing 158, these objects were missing:",
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
    "science_condensation_trials_pubregions_cartoons"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PREVIOUSLY PUBLISHED REGIONS
#
# IMPORTANT:
# Replace xmin/xmax below with the exact Mb boundaries you want.
#
# These should be in Mb on Hmr18, matching the 158 x-axis.
#
# The layout is deliberately thin and schematic.
# ============================================================

published_regions <- data.frame(
    region = c(
        "B",
        "R1",
        "R2",
        "D1",
        "D2"
    ),

    module = c(
        "B",
        "R",
        "R",
        "D",
        "D"
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

    colour = c(
        trait_colours[["B"]],
        trait_colours[["R"]],
        trait_colours[["R"]],
        trait_colours[["D"]],
        trait_colours[["D"]]
    ),

    stringsAsFactors = FALSE
)

if (any(!is.finite(published_regions$xmin)) ||
    any(!is.finite(published_regions$xmax))) {

    cat(
        "\nFill in the published_regions xmin/xmax values in:\n",
        "  scripts/159_plot_science_trials_pubregions_cartoons.R\n\n",
        "The rest of the script is ready.\n",
        sep = ""
    )

    print(
        published_regions,
        row.names = FALSE
    )

    stop(
        "published_regions table still contains NA bounds."
    )
}


# ============================================================
# TOP-OF-PANEL REGION BAR POSITIONS
# ============================================================

pub_trial1 <- published_regions

pub_trial1$y <- c(
    4.18,   # B
    4.18,   # R1
    4.18,   # R2
    4.18,   # D1
    4.06    # D2; offset because D1/D2 are adjacent
)

pub_trial1$label_y <- pub_trial1$y + 0.045


pub_trial2 <- published_regions

pub_trial2$y <- c(
    3.27,   # B
    3.27,   # R1
    3.27,   # R2
    3.27,   # D1
    3.16    # D2; offset because D1/D2 are adjacent
)

pub_trial2$label_y <- pub_trial2$y + 0.040

# ============================================================
# BUTTERFLY PATTERN CARTOONS
#
# These are intentionally schematic rather than detailed.
# The goal is recognition of the pattern logic, not realism.
# ============================================================

make_butterfly_base <- function() {

    left_fore <- polygonGrob(
        x = unit(
            c(
                0.47, 0.32, 0.16, 0.10, 0.17, 0.29, 0.43
            ),
            "npc"
        ),
        y = unit(
            c(
                0.60, 0.84, 0.88, 0.70, 0.52, 0.42, 0.50
            ),
            "npc"
        ),
        gp = gpar(
            fill = "#F7F7F7",
            col = "#4A4A4A",
            lwd = 1
        )
    )

    right_fore <- polygonGrob(
        x = unit(
            1 - c(
                0.47, 0.32, 0.16, 0.10, 0.17, 0.29, 0.43
            ),
            "npc"
        ),
        y = unit(
            c(
                0.60, 0.84, 0.88, 0.70, 0.52, 0.42, 0.50
            ),
            "npc"
        ),
        gp = gpar(
            fill = "#F7F7F7",
            col = "#4A4A4A",
            lwd = 1
        )
    )

    left_hind <- polygonGrob(
        x = unit(
            c(
                0.45, 0.29, 0.16, 0.18, 0.30, 0.42
            ),
            "npc"
        ),
        y = unit(
            c(
                0.42, 0.34, 0.20, 0.08, 0.12, 0.28
            ),
            "npc"
        ),
        gp = gpar(
            fill = "#FAFAFA",
            col = "#4A4A4A",
            lwd = 1
        )
    )

    right_hind <- polygonGrob(
        x = unit(
            1 - c(
                0.45, 0.29, 0.16, 0.18, 0.30, 0.42
            ),
            "npc"
        ),
        y = unit(
            c(
                0.42, 0.34, 0.20, 0.08, 0.12, 0.28
            ),
            "npc"
        ),
        gp = gpar(
            fill = "#FAFAFA",
            col = "#4A4A4A",
            lwd = 1
        )
    )

    body <- segmentsGrob(
        x0 = unit(0.50, "npc"),
        x1 = unit(0.50, "npc"),
        y0 = unit(0.12, "npc"),
        y1 = unit(0.72, "npc"),
        gp = gpar(
            col = "#333333",
            lwd = 3
        )
    )

    grobTree(
        left_fore,
        right_fore,
        left_hind,
        right_hind,
        body
    )
}


make_band_cartoon <- function() {

    grobTree(
        make_butterfly_base(),

        # Schematic "band" on both forewings
        segmentsGrob(
            x0 = unit(c(0.18, 0.82), "npc"),
            x1 = unit(c(0.43, 0.57), "npc"),
            y0 = unit(c(0.63, 0.63), "npc"),
            y1 = unit(c(0.52, 0.52), "npc"),
            gp = gpar(
                col = trait_colours[["B"]],
                lwd = 5
            )
        ),

        textGrob(
            label = "Band expectation",
            x = unit(0.5, "npc"),
            y = unit(-0.06, "npc"),
            gp = gpar(
                cex = 0.7,
                col = "#333333"
            )
        )
    )
}


make_DR_cartoon <- function() {

    grobTree(
        make_butterfly_base(),

        # Schematic Dennis patch (orange)
        segmentsGrob(
            x0 = unit(c(0.33, 0.67), "npc"),
            x1 = unit(c(0.45, 0.55), "npc"),
            y0 = unit(c(0.77, 0.77), "npc"),
            y1 = unit(c(0.66, 0.66), "npc"),
            gp = gpar(
                col = trait_colours[["D"]],
                lwd = 5
            )
        ),

        # Schematic Ray streaks (magenta)
        segmentsGrob(
            x0 = unit(c(0.22, 0.78, 0.25, 0.75), "npc"),
            x1 = unit(c(0.40, 0.60, 0.35, 0.65), "npc"),
            y0 = unit(c(0.33, 0.33, 0.22, 0.22), "npc"),
            y1 = unit(c(0.46, 0.46, 0.38, 0.38), "npc"),
            gp = gpar(
                col = trait_colours[["R"]],
                lwd = 3.5
            )
        ),

        textGrob(
            label = "Dennis / Ray expectation",
            x = unit(0.5, "npc"),
            y = unit(-0.06, "npc"),
            gp = gpar(
                cex = 0.7,
                col = "#333333"
            )
        )
    )
}

band_cartoon <- make_band_cartoon()
DR_cartoon <- make_DR_cartoon()


# ============================================================
# COMMON REGION-ANNOTATION LAYERS
# ============================================================

add_published_region_layers <- function(p, region_df) {

    p +
        geom_segment(
            data = region_df,
            aes(
                x = xmin,
                xend = xmax,
                y = y,
                yend = y
            ),
            inherit.aes = FALSE,
            linewidth = 2.2,
            lineend = "butt",
            colour = region_df$colour
        ) +
        geom_text(
            data = region_df,
            aes(
                x = (xmin + xmax) / 2,
                y = label_y,
                label = region
            ),
            inherit.aes = FALSE,
            size = 2.25,
            vjust = 0,
            colour = region_df$colour
        )
}


# ============================================================
# TRIAL 1 WITH REGIONS + CARTOONS
# ============================================================

trial1_pub <- trial1

trial1_pub <- add_published_region_layers(
    trial1_pub,
    pub_trial1
)

trial1_pub <- trial1_pub +
    annotation_custom(
        grob = band_cartoon,
        xmin = science_x_limits[1] - 0.026,
        xmax = science_x_limits[1] - 0.006,
        ymin = -0.19,
        ymax = 0.36
    ) +
    annotation_custom(
        grob = DR_cartoon,
        xmin = science_x_limits[1] - 0.026,
        xmax = science_x_limits[1] - 0.006,
        ymin = -0.56,
        ymax = -0.02
    ) +
    labs(
        subtitle = paste(
            "Published intervals are shown as thin bars;",
            "lower cartoons indicate the expected Band or Dennis/Ray pattern class."
        )
    )


# ============================================================
# TRIAL 2 WITH REGIONS + CARTOONS
# ============================================================

trial2_pub <- trial2

trial2_pub <- add_published_region_layers(
    trial2_pub,
    pub_trial2
)

trial2_pub <- trial2_pub +
    annotation_custom(
        grob = band_cartoon,
        xmin = science_x_limits[1] - 0.026,
        xmax = science_x_limits[1] - 0.006,
        ymin = -0.27,
        ymax = 0.28
    ) +
    annotation_custom(
        grob = DR_cartoon,
        xmin = science_x_limits[1] - 0.026,
        xmax = science_x_limits[1] - 0.006,
        ymin = -0.61,
        ymax = -0.08
    ) +
    labs(
        subtitle = paste(
            "Published intervals are shown as thin bars;",
            "lower cartoons indicate the expected Band or Dennis/Ray pattern class."
        )
    )


# ============================================================
# SAVE
# ============================================================

trial1_pdf <- file.path(
    output_dir,
    "Science_trial_1_three_compact_BDR_BA_tracks.pubregions_cartoons.pdf"
)

trial1_png <- file.path(
    output_dir,
    "Science_trial_1_three_compact_BDR_BA_tracks.pubregions_cartoons.png"
)

trial2_pdf <- file.path(
    output_dir,
    "Science_trial_2_BDR_BA_rug_tracks.pubregions_cartoons.pdf"
)

trial2_png <- file.path(
    output_dir,
    "Science_trial_2_BDR_BA_rug_tracks.pubregions_cartoons.png"
)

comparison_pdf <- file.path(
    output_dir,
    "Science_trials_1_2_pubregions_cartoons.comparison.pdf"
)

ggsave(
    filename = trial1_pdf,
    plot = trial1_pub,
    width = 7.4,
    height = 4.6,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = trial1_png,
    plot = trial1_pub,
    width = 7.4,
    height = 4.6,
    units = "in",
    dpi = 450,
    limitsize = FALSE
)

ggsave(
    filename = trial2_pdf,
    plot = trial2_pub,
    width = 7.4,
    height = 3.8,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
)

ggsave(
    filename = trial2_png,
    plot = trial2_pub,
    width = 7.4,
    height = 3.8,
    units = "in",
    dpi = 450,
    limitsize = FALSE
)

grDevices::cairo_pdf(
    filename = comparison_pdf,
    width = 7.4,
    height = 4.8,
    onefile = TRUE
)

print(trial1_pub)
print(trial2_pub)

grDevices::dev.off()


# ============================================================
# SUMMARY
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 159 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "Trial 1:\n",
    trial1_pdf,
    "\n\n",
    "Trial 2:\n",
    trial2_pdf,
    "\n\n",
    "Comparison:\n",
    comparison_pdf,
    "\n",
    sep = ""
)
