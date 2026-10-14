#!/usr/bin/env Rscript

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_file <- file.path(
    project,
    "scripts",
    "159_plot_science_trials_pubregions_cartoons.R"
)

if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing working Analysis-159 script:",
            source_file
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
})


# ============================================================
# LOAD THE WORKING ANALYSIS-159 OBJECTS
#
# Suppress the inherited 149/158 missing-point warnings here;
# they were already audited and are unrelated to this layout
# derivative.
# ============================================================

suppressWarnings(
    source(
        source_file,
        local = TRUE
    )
)

required_objects <- c(
    "trial1_pub",
    "trial2_pub",
    "trait_colours",
    "make_butterfly_base"
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
    "science_condensation_trials_cartoon_sidecar"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# REMOVE THE OLD OUT-OF-BOUNDS CARTOONS
#
# Also remove the two focal-pattern text annotations that 158
# deliberately placed outside x=1.06. Those labels will move
# into the new fixed-width sidecar.
# ============================================================

get_layer_labels <- function(layer) {

    labels <- character(0)

    if (
        is.data.frame(layer$data) &&
        "label" %in% names(layer$data)
    ) {
        labels <- c(
            labels,
            as.character(
                layer$data$label
            )
        )
    }

    if (
        !is.null(
            layer$aes_params$label
        )
    ) {
        labels <- c(
            labels,
            as.character(
                layer$aes_params$label
            )
        )
    }

    if (
        !is.null(
            layer$geom_params$label
        )
    ) {
        labels <- c(
            labels,
            as.character(
                layer$geom_params$label
            )
        )
    }

    labels
}


is_old_cartoon_or_focal_label <- function(layer) {

    geom_class <- class(
        layer$geom
    )[1]

    # Both old butterfly cartoons came from annotation_custom().
    if (
        identical(
            geom_class,
            "GeomCustomAnn"
        )
    ) {
        return(TRUE)
    }

    labels <- get_layer_labels(
        layer
    )

    if (length(labels) == 0L) {
        return(FALSE)
    }

    any(
        grepl(
            "^B\\s+A/A\\s*->\\s*P/A",
            labels
        )
    ) ||
    any(
        grepl(
            "^D/R\\s+P/P\\s*->\\s*P/A",
            labels
        )
    )
}


strip_old_side_annotations <- function(p) {

    remove <- vapply(
        p$layers,
        is_old_cartoon_or_focal_label,
        logical(1)
    )

    cat(
        "Removing ",
        sum(remove),
        " old out-of-bounds layer(s).\n",
        sep = ""
    )

    p$layers <- p$layers[
        !remove
    ]

    p
}


trial1_clean <- strip_old_side_annotations(
    trial1_pub
)

trial2_clean <- strip_old_side_annotations(
    trial2_pub
)


# ============================================================
# MAKE SMALL CARTOONS FOR THE SIDECAR
#
# Everything is drawn INSIDE the sidecar's own viewport.
# Nothing uses genomic x coordinates.
# ============================================================

make_band_icon <- function() {

    grobTree(

        make_butterfly_base(),

        segmentsGrob(
            x0 = unit(
                c(
                    0.18,
                    0.82
                ),
                "npc"
            ),

            x1 = unit(
                c(
                    0.43,
                    0.57
                ),
                "npc"
            ),

            y0 = unit(
                c(
                    0.63,
                    0.63
                ),
                "npc"
            ),

            y1 = unit(
                c(
                    0.52,
                    0.52
                ),
                "npc"
            ),

            gp = gpar(
                col = trait_colours[["B"]],
                lwd = 5
            )
        )
    )
}


make_DR_icon <- function() {

    grobTree(

        make_butterfly_base(),

        # Dennis forewing patch
        segmentsGrob(
            x0 = unit(
                c(
                    0.33,
                    0.67
                ),
                "npc"
            ),

            x1 = unit(
                c(
                    0.45,
                    0.55
                ),
                "npc"
            ),

            y0 = unit(
                c(
                    0.77,
                    0.77
                ),
                "npc"
            ),

            y1 = unit(
                c(
                    0.66,
                    0.66
                ),
                "npc"
            ),

            gp = gpar(
                col = trait_colours[["D"]],
                lwd = 5
            )
        ),

        # Ray streaks
        segmentsGrob(
            x0 = unit(
                c(
                    0.22,
                    0.78,
                    0.25,
                    0.75
                ),
                "npc"
            ),

            x1 = unit(
                c(
                    0.40,
                    0.60,
                    0.35,
                    0.65
                ),
                "npc"
            ),

            y0 = unit(
                c(
                    0.33,
                    0.33,
                    0.22,
                    0.22
                ),
                "npc"
            ),

            y1 = unit(
                c(
                    0.46,
                    0.46,
                    0.38,
                    0.38
                ),
                "npc"
            ),

            gp = gpar(
                col = trait_colours[["R"]],
                lwd = 3.5
            )
        )
    )
}


band_icon <- make_band_icon()
DR_icon <- make_DR_icon()


# ============================================================
# SIDECAR
#
# The y positions are chosen to line up with the bottom B and
# D/R focal-genotype rugs in BOTH trial 1 and trial 2.
#
# These are panel-relative coordinates, not genomic coordinates.
# ============================================================

make_focal_sidecar <- function() {

    grobTree(

        # ----------------------------------------------------
        # BAND
        # ----------------------------------------------------

        grobTree(
            band_icon,

            vp = viewport(
                x = 0.18,
                y = 0.119,
                width = 0.27,
                height = 0.090
            )
        ),

        textGrob(
            "Band",
            x = unit(
                0.36,
                "npc"
            ),
            y = unit(
                0.132,
                "npc"
            ),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "A/A -> P/A",
            x = unit(
                0.36,
                "npc"
            ),
            y = unit(
                0.102,
                "npc"
            ),
            hjust = 0,
            gp = gpar(
                fontsize = 6.3,
                col = "#4A4A4A"
            )
        ),


        # ----------------------------------------------------
        # DENNIS / RAY
        # ----------------------------------------------------

        grobTree(
            DR_icon,

            vp = viewport(
                x = 0.18,
                y = 0.050,
                width = 0.27,
                height = 0.090
            )
        ),

        textGrob(
            "Dennis / Ray",
            x = unit(
                0.36,
                "npc"
            ),
            y = unit(
                0.063,
                "npc"
            ),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "P/P -> P/A",
            x = unit(
                0.36,
                "npc"
            ),
            y = unit(
                0.033,
                "npc"
            ),
            hjust = 0,
            gp = gpar(
                fontsize = 6.3,
                col = "#4A4A4A"
            )
        )
    )
}


# ============================================================
# ATTACH THE SIDECAR TO THE GGPLOT PANEL
#
# Critically:
#
#   * the genomic panel itself is not widened or rescaled;
#   * the sidecar receives its own fixed physical width;
#   * it spans exactly the same gtable rows as the plotting
#     panel, so the cartoons align with the lower rug tracks;
#   * no content depends on clip="off".
# ============================================================

add_panel_sidecar <- function(
    p,
    sidecar_width = unit(
        1.38,
        "in"
    )
) {

    g <- ggplotGrob(
        p
    )

    panel_rows <- which(
        grepl(
            "^panel",
            g$layout$name
        )
    )

    if (length(panel_rows) == 0L) {
        stop(
            "Could not find panel in ggplot gtable."
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


    # Insert a new column immediately before the panel.
    g <- gtable_add_cols(
        g,
        widths = sidecar_width,
        pos = panel_l - 1L
    )


    # Newly inserted column occupies the old panel_l position.
    sidecar_column <- panel_l

    g <- gtable_add_grob(
        g,
        grobs = make_focal_sidecar(),
        t = panel_t,
        b = panel_b,
        l = sidecar_column,
        r = sidecar_column,
        clip = "off",
        name = "focal-pattern-sidecar"
    )

    g
}


# ============================================================
# CLEAN SUBTITLES
# ============================================================

trial1_clean <- trial1_clean +
    labs(
        subtitle = paste(
            "Published B, R1/R2 and D1/D2 intervals are shown above;",
            "focal pattern expectations are shown at left."
        )
    )

trial2_clean <- trial2_clean +
    labs(
        subtitle = paste(
            "Published B, R1/R2 and D1/D2 intervals are shown above;",
            "focal pattern expectations are shown at left."
        )
    )


trial1_grob <- add_panel_sidecar(
    trial1_clean
)

trial2_grob <- add_panel_sidecar(
    trial2_clean
)


# ============================================================
# SAVE HELPER
# ============================================================

draw_grob_pdf <- function(
    grob,
    filename,
    width,
    height
) {

    grDevices::cairo_pdf(
        filename = filename,
        width = width,
        height = height,
        onefile = FALSE
    )

    grid.newpage()

    grid.draw(
        grob
    )

    grDevices::dev.off()
}


draw_grob_png <- function(
    grob,
    filename,
    width,
    height,
    dpi = 450
) {

    grDevices::png(
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

    grDevices::dev.off()
}


# ============================================================
# SAVE INDIVIDUAL FIGURES
#
# Width is slightly larger than 159 because the cartoon column
# is now REAL page space rather than content hanging outside it.
#
# The genomic plotting width remains essentially the same.
# ============================================================

trial1_pdf <- file.path(
    output_dir,
    "Science_trial_1_BDR_BA_tracks.pubregions_cartoon_sidecar.pdf"
)

trial1_png <- file.path(
    output_dir,
    "Science_trial_1_BDR_BA_tracks.pubregions_cartoon_sidecar.png"
)

trial2_pdf <- file.path(
    output_dir,
    "Science_trial_2_BDR_BA_rugs.pubregions_cartoon_sidecar.pdf"
)

trial2_png <- file.path(
    output_dir,
    "Science_trial_2_BDR_BA_rugs.pubregions_cartoon_sidecar.png"
)


draw_grob_pdf(
    trial1_grob,
    trial1_pdf,
    width = 8.55,
    height = 4.60
)

draw_grob_png(
    trial1_grob,
    trial1_png,
    width = 8.55,
    height = 4.60
)

draw_grob_pdf(
    trial2_grob,
    trial2_pdf,
    width = 8.55,
    height = 3.80
)

draw_grob_png(
    trial2_grob,
    trial2_png,
    width = 8.55,
    height = 3.80
)


# ============================================================
# MULTIPAGE COMPARISON
# ============================================================

comparison_pdf <- file.path(
    output_dir,
    "Science_trials_1_2_cartoon_sidecar.comparison.pdf"
)

grDevices::cairo_pdf(
    filename = comparison_pdf,
    width = 8.55,
    height = 4.75,
    onefile = TRUE
)

grid.newpage()
grid.draw(
    trial1_grob
)

grid.newpage()
grid.draw(
    trial2_grob
)

grDevices::dev.off()


# ============================================================
# SUMMARY
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 160 COMPLETE\n",
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
