#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
})


# ============================================================
# ANALYSIS 170
#
# Correct focal biological model:
#
#       LEFT  = DR/DR
#       RIGHT = B/DR
#
# This is the ONLY focal genotype pattern treated as supportive.
#
# No B/B focal state.
# No reverse orientation.
# No alternative "Band-compatible" focal genotype definition.
#
# Upper panel:
#   established separate B / D / R comparative BA rugs.
#
# Lower panel:
#   one L row
#   one R row
#
# 170A:
#   only focal-pattern sites that also have BA > 0.99
#   for at least one of B, D, or R.
#
# 170B:
#   all exact focal-pattern sites faintly;
#   BA > 0.99 subset emphasized.
#
# Lower-panel genotype grammar:
#
#   L DR/DR = one DR-coloured tick
#   R B/DR  = paired B + DR ticks
#
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


# ============================================================
# LOAD ESTABLISHED DATA OBJECTS
# ============================================================

suppressWarnings(
    source(
        source_file,
        local = TRUE
    )
)


required_objects <- c(
    "all_ba_science",
    "dr_correct",
    "trait_colours",
    "science_x_limits",
    "science_x_breaks",
    "optix_start_science",
    "optix_end_science",
    "band_cartoon",
    "DR_cartoon"
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
    "science_exact_L_DRDR_R_BDR_homhetcolours"
)


dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


pdf_A <- file.path(
    output_dir,
    "Science_170A_exact_L_DRDR_R_BDR_highBA_only.pdf"
)


png_A <- file.path(
    output_dir,
    "Science_170A_exact_L_DRDR_R_BDR_highBA_only.png"
)


pdf_B <- file.path(
    output_dir,
    "Science_170B_exact_L_DRDR_R_BDR_all_plus_highBA.pdf"
)


png_B <- file.path(
    output_dir,
    "Science_170B_exact_L_DRDR_R_BDR_all_plus_highBA.png"
)


comparison_pdf <- file.path(
    output_dir,
    "Science_170A_170B_comparison.pdf"
)


audit_file <- file.path(
    output_dir,
    "Science_170_exact_L_DRDR_R_BDR.audit.tsv"
)


summary_file <- file.path(
    output_dir,
    "Science_170_exact_L_DRDR_R_BDR.summary.tsv"
)


# ============================================================
# HELPERS
# ============================================================

normalise_name <- function(x) {

    gsub(
        "[^a-z0-9]+",
        "",
        tolower(
            as.character(x)
        ),
        perl = TRUE
    )
}


resolve_column <- function(
    x,
    candidates
) {

    existing <- normalise_name(
        names(x)
    )


    for (candidate in candidates) {

        target <- normalise_name(
            candidate
        )

        hit <- which(
            existing == target
        )


        if (length(hit) == 1L) {

            return(
                names(x)[hit]
            )
        }
    }


    stop(
        paste0(
            "Could not resolve any of: ",
            paste(
                candidates,
                collapse = ", "
            ),
            "\nAvailable columns:\n",
            paste(
                names(x),
                collapse = "\n"
            )
        )
    )
}


extract_position_Mb <- function(x) {

    position_column <- resolve_column(
        x,
        c(
            "position_Mb",
            "position_mb",
            "pos_Mb"
        )
    )


    suppressWarnings(
        as.numeric(
            x[[position_column]]
        )
    )
}


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


# ============================================================
# COLOURS
# ============================================================

# Supplement genoplot genotype colours
HOM_AA_colour <- "#2C9BCB"
HET_colour    <- "#F28E2B"
HOM_PP_colour <- "#B61F72"

# Focal biological states used here
DR_colour <- HOM_PP_colour

neutral_dark <- "#66707A"

neutral_mid <- "#AEB6BE"

neutral_light <- "#D9DEE3"

neutral_faint <- "#C9CED4"

text_colour <- "#242424"

optix_colour <- "#3A7CA5"


# ============================================================
# UPPER B / D / R BA TRACKS
# ============================================================

ba_trait_column <- resolve_column(
    all_ba_science,
    c(
        "trait"
    )
)


ba_position_column <- resolve_column(
    all_ba_science,
    c(
        "position_Mb",
        "position_mb"
    )
)


ba_value_column <- resolve_column(
    all_ba_science,
    c(
        "BA",
        "balanced_accuracy"
    )
)


ba <- data.frame(

    trait = as.character(
        all_ba_science[[ba_trait_column]]
    ),

    position_Mb = suppressWarnings(
        as.numeric(
            all_ba_science[[ba_position_column]]
        )
    ),

    BA = suppressWarnings(
        as.numeric(
            all_ba_science[[ba_value_column]]
        )
    ),

    stringsAsFactors = FALSE
)


ba <- ba[
    is.finite(
        ba$position_Mb
    ) &
    is.finite(
        ba$BA
    ) &
    ba$position_Mb >=
        science_x_limits[1] &
    ba$position_Mb <=
        science_x_limits[2],
    ,
    drop = FALSE
]


ba$pos <- as.integer(
    round(
        ba$position_Mb *
            1e6
    )
)


ba$trait <- factor(
    ba$trait,
    levels = c(
        "B",
        "D",
        "R"
    ),
    ordered = TRUE
)


rug_bases <- c(
    B = 5.04,
    D = 4.27,
    R = 3.50
)


ba$rug_base <- unname(
    rug_bases[
        as.character(
            ba$trait
        )
    ]
)


ba$tick_height <- ifelse(

    ba$BA > 0.99,
    0.53,

    ifelse(

        ba$BA >= 0.90,
        0.37,

        ifelse(

            ba$BA >= 0.80,
            0.23,

            0.085
        )
    )
)


ba$rug_top <- (
    ba$rug_base +
    ba$tick_height
)


high_ba <- ba[
    ba$BA > 0.99,
    ,
    drop = FALSE
]


# ============================================================
# EXACT FOCAL SUPPORTIVE POSITIONS
#
# IMPORTANT:
#
# dr_correct is the established exact P/P -> P/A focal set.
#
# Under the corrected biological naming requested here:
#
#       P/P -> P/A
#
# is represented in this focal model as:
#
#       LEFT DR/DR -> RIGHT B/DR
#
# This is the ONLY focal-support definition used below.
#
# We deliberately do NOT use band_correct.
# ============================================================

exact_position_Mb <- extract_position_Mb(
    dr_correct
)


exact_position_Mb <- exact_position_Mb[
    is.finite(
        exact_position_Mb
    ) &
    exact_position_Mb >=
        science_x_limits[1] &
    exact_position_Mb <=
        science_x_limits[2]
]


exact_focal <- data.frame(

    position_Mb = exact_position_Mb,

    stringsAsFactors = FALSE
)


exact_focal$pos <- as.integer(
    round(
        exact_focal$position_Mb *
            1e6
    )
)


exact_focal <- exact_focal[
    order(
        exact_focal$position_Mb
    ),
    ,
    drop = FALSE
]


exact_focal <- exact_focal[
    !duplicated(
        exact_focal$pos
    ),
    ,
    drop = FALSE
]


if (nrow(exact_focal) == 0L) {

    stop(
        paste(
            "No exact Left DR/DR -> Right B/DR",
            "focal positions were recovered."
        )
    )
}


# ============================================================
# HIGH-BA SUPPORT AT EXACT FOCAL POSITIONS
#
# Do NOT inherit high_B/high_D/high_R focal definitions here.
#
# Instead:
#
#   1. take exact focal pattern positions from dr_correct;
#   2. take BA > 0.99 independently from each comparative track;
#   3. intersect by genomic position.
#
# Thus the focal genotype rule is identical for B, D, and R.
# ============================================================

exact_key <- exact_focal$pos


high_exact <- high_ba[
    high_ba$pos %in%
        exact_key,
    ,
    drop = FALSE
]


high_exact <- high_exact[
    order(
        high_exact$position_Mb,
        high_exact$trait
    ),
    ,
    drop = FALSE
]


B_high <- high_exact[
    as.character(
        high_exact$trait
    ) == "B",
    ,
    drop = FALSE
]


D_high <- high_exact[
    as.character(
        high_exact$trait
    ) == "D",
    ,
    drop = FALSE
]


R_high <- high_exact[
    as.character(
        high_exact$trait
    ) == "R",
    ,
    drop = FALSE
]


high_union_positions <- sort(
    unique(
        high_exact$pos
    )
)


high_union <- exact_focal[
    exact_focal$pos %in%
        high_union_positions,
    ,
    drop = FALSE
]


# ============================================================
# SUPPORT CLASS PER EXACT POSITION
#
# This is written to the audit table.
#
# It is NOT needed as a separate lower-panel colour encoding;
# the B/D/R identity can already be read directly from the
# aligned upper BA tracks.
# ============================================================

support_for_position <- function(position) {

    traits <- as.character(
        high_exact$trait[
            high_exact$pos ==
                position
        ]
    )


    traits <- unique(
        traits
    )


    if (length(traits) == 0L) {

        return(
            "none"
        )
    }


    paste(
        traits,
        collapse = "+"
    )
}


exact_focal$high_BA_support <- vapply(
    exact_focal$pos,
    support_for_position,
    character(1)
)


exact_focal$high_BA_any <- (
    exact_focal$high_BA_support !=
        "none"
)


# ============================================================
# LOWER L / R GEOMETRY
#
# One physical L row.
# One physical R row.
#
# LEFT DR/DR:
#     single magenta tick
#
# RIGHT B/DR:
#     single heterozygote-coloured tick
#
# No connector forest.
# No midpoint symbols.
# No B/B state.
# ============================================================

left_y <- 1.78

right_y <- 1.24


background_half_height <- 0.060

strong_half_height <- 0.125

pair_offset <- 0


# ------------------------------------------------------------
# ALL exact-pattern positions
# ------------------------------------------------------------

all_left <- data.frame(

    position_Mb = exact_focal$position_Mb,

    y0 = left_y -
        background_half_height,

    y1 = left_y +
        background_half_height,

    stringsAsFactors = FALSE
)


all_right_B <- data.frame(

    position_Mb = exact_focal$position_Mb -
        pair_offset,

    y0 = right_y -
        background_half_height,

    y1 = right_y +
        background_half_height,

    stringsAsFactors = FALSE
)


all_right_DR <- data.frame(

    position_Mb = exact_focal$position_Mb +
        pair_offset,

    y0 = right_y -
        background_half_height,

    y1 = right_y +
        background_half_height,

    stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# HIGH-BA subset
# ------------------------------------------------------------

high_left <- data.frame(

    position_Mb = high_union$position_Mb,

    y0 = left_y -
        strong_half_height,

    y1 = left_y +
        strong_half_height,

    stringsAsFactors = FALSE
)


high_right_B <- data.frame(

    position_Mb = high_union$position_Mb -
        pair_offset,

    y0 = right_y -
        strong_half_height,

    y1 = right_y +
        strong_half_height,

    stringsAsFactors = FALSE
)


high_right_DR <- data.frame(

    position_Mb = high_union$position_Mb +
        pair_offset,

    y0 = right_y -
        strong_half_height,

    y1 = right_y +
        strong_half_height,

    stringsAsFactors = FALSE
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


published_regions$y <- 5.86

published_regions$label_y <- 5.97


# ============================================================
# AUDIT
# ============================================================

audit <- exact_focal[
    ,
    c(
        "pos",
        "position_Mb",
        "high_BA_support",
        "high_BA_any"
    ),
    drop = FALSE
]


audit$left_state <- "DR/DR"

audit$right_state <- "B/DR"


audit <- audit[
    ,
    c(
        "pos",
        "position_Mb",
        "left_state",
        "right_state",
        "high_BA_support",
        "high_BA_any"
    ),
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

    statistic = c(
        "exact_L_DRDR_R_BDR_sites",
        "exact_sites_with_any_BA_gt_0.99",
        "B_BA_gt_0.99_at_exact_pattern",
        "D_BA_gt_0.99_at_exact_pattern",
        "R_BA_gt_0.99_at_exact_pattern",
        "unique_high_BA_positions"
    ),

    value = c(
        nrow(
            exact_focal
        ),
        sum(
            exact_focal$high_BA_any
        ),
        nrow(
            B_high
        ),
        nrow(
            D_high
        ),
        nrow(
            R_high
        ),
        length(
            high_union_positions
        )
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
# BASE UPPER PANEL + LOWER FRAME
# ============================================================

build_base_plot <- function(
    title_text,
    subtitle_text,
    caption_text
) {

    ggplot() +


        # ----------------------------------------------------
        # Published intervals
        # ----------------------------------------------------

        geom_segment(
            data = published_regions,

            aes(
                x = xmin,
                xend = xmax,
                y = y,
                yend = y,
                colour = trait
            ),

            linewidth = 2.30,
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

            size = 2.28,
            fontface = "bold",
            vjust = 0,
            show.legend = FALSE
        ) +


        # ----------------------------------------------------
        # optix
        # ----------------------------------------------------

        geom_segment(

            aes(
                x = optix_start_science,
                xend = optix_end_science,
                y = 5.53,
                yend = 5.53
            ),

            inherit.aes = FALSE,

            colour = optix_colour,

            linewidth = 2.90,

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

            y = 5.62,

            label = "optix",

            colour = optix_colour,

            size = 2.42,

            fontface = "italic",

            vjust = 0
        ) +


        # ----------------------------------------------------
        # BA baselines
        # ----------------------------------------------------

        geom_segment(
            data = data.frame(
                trait = c(
                    "B",
                    "D",
                    "R"
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

            colour = neutral_mid,

            linewidth = 0.28
        ) +


        # ----------------------------------------------------
        # All comparative BA
        # ----------------------------------------------------

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

            alpha = 0.32,

            linewidth = 0.23
        ) +


        # ----------------------------------------------------
        # BA > 0.99
        # ----------------------------------------------------

        geom_segment(
            data = high_ba,

            aes(
                x = position_Mb,
                xend = position_Mb,
                y = rug_base,
                yend = rug_top,
                colour = trait
            ),

            inherit.aes = FALSE,

            linewidth = 0.75,

            alpha = 1,

            lineend = "butt",

            show.legend = FALSE
        ) +


        # ----------------------------------------------------
        # B / D / R labels
        # ----------------------------------------------------

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
                    0.25,

                stringsAsFactors = FALSE
            ),

            aes(
                x = x,
                y = y,
                label = trait,
                colour = trait
            ),

            hjust = 0,

            fontface = "bold",

            size = 3.0,

            show.legend = FALSE
        ) +


        # ----------------------------------------------------
        # LOWER-PANEL SEPARATOR
        # ----------------------------------------------------

        geom_segment(

            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = 2.53,
                yend = 2.53
            ),

            inherit.aes = FALSE,

            colour = "#D4D7DA",

            linewidth = 0.35
        ) +


        annotate(
            "text",

            x = science_x_limits[1] +
                0.0015,

            y = 2.64,

            label = "Focal L/R genotype",

            hjust = 0,

            vjust = 0,

            size = 2.40,

            fontface = "bold",

            colour = "#505050"
        ) +


        annotate(
            "text",

            x = science_x_limits[2] -
                0.0010,

            y = 2.64,

            label = "supportive pattern: L = DR/DR, R = B/DR",

            hjust = 1,

            vjust = 0,

            size = 2.05,

            colour = "#666666"
        ) +


        # ----------------------------------------------------
        # L / R ROWS
        # ----------------------------------------------------

        geom_segment(

            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = left_y,
                yend = left_y
            ),

            inherit.aes = FALSE,

            colour = neutral_light,

            linewidth = 0.30
        ) +


        geom_segment(

            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = right_y,
                yend = right_y
            ),

            inherit.aes = FALSE,

            colour = neutral_light,

            linewidth = 0.30
        ) +


        annotate(
            "text",

            x = science_x_limits[1] +
                0.0013,

            y = left_y +
                0.17,

            label = "L",

            hjust = 0,

            fontface = "bold",

            size = 3.45,

            colour = "#333333"
        ) +


        annotate(
            "text",

            x = science_x_limits[1] +
                0.0013,

            y = right_y +
                0.17,

            label = "R",

            hjust = 0,

            fontface = "bold",

            size = 3.45,

            colour = "#333333"
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
                0.95,
                6.15
            ),

            clip = "off"
        ) +


        labs(
            title = title_text,
            subtitle = subtitle_text,
            caption = caption_text
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
                size = 7.0,
                colour = "#505050",
                margin = margin(
                    b = 4
                )
            ),

            plot.caption = element_text(
                size = 6.3,
                colour = "#555555",
                hjust = 0,
                margin = margin(
                    t = 4
                )
            ),

            axis.title.x = element_text(
                size = 7.4,
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
}


# ============================================================
# CARTOON / GLYPH SIDECAR
#
# Only ONE supportive focal pattern is shown:
#
#       L  DR/DR
#       R  B/DR
#
# For B/DR we show the B cartoon + DR cartoon side by side,
# rather than pretending there is a separate B/B focal class.
# ============================================================

make_supportive_key <- function() {

    grobTree(

        textGrob(
            "Supportive focal SNP",
            x = unit(0.06, "npc"),
            y = unit(0.94, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8.2,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "Only this ordered L/R genotype pattern is highlighted",
            x = unit(0.06, "npc"),
            y = unit(0.875, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.7,
                col = "#666666"
            )
        ),

        # ====================================================
        # LEFT
        # ====================================================

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.69, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 8,
                fontface = "bold",
                col = "#333333"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.64, "npc"),
            y1 = unit(0.74, "npc"),
            gp = gpar(
                col = DR_colour,
                lwd = 3.5
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.26,
                y = 0.69,
                width = 0.20,
                height = 0.16
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.39, "npc"),
            y = unit(0.71, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = DR_colour
            )
        ),

        textGrob(
            "DR only",
            x = unit(0.39, "npc"),
            y = unit(0.65, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        # ====================================================
        # Same SNP / transition
        # ====================================================

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.58, "npc"),
            y1 = unit(0.49, "npc"),
            arrow = arrow(
                length = unit(0.055, "in"),
                type = "closed"
            ),
            gp = gpar(
                col = "#808080",
                lwd = 0.9
            )
        ),

        textGrob(
            "same SNP",
            x = unit(0.18, "npc"),
            y = unit(0.535, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.7,
                col = "#777777"
            )
        ),

        # ====================================================
        # RIGHT
        # ====================================================

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.39, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 8,
                fontface = "bold",
                col = "#333333"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.34, "npc"),
            y1 = unit(0.44, "npc"),
            gp = gpar(
                col = HET_colour,
                lwd = 3.5
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.22,
                y = 0.39,
                width = 0.14,
                height = 0.11
            )
        ),

        textGrob(
            "+",
            x = unit(0.285, "npc"),
            y = unit(0.39, "npc"),
            gp = gpar(
                fontsize = 7.5,
                fontface = "bold",
                col = "#555555"
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.35,
                y = 0.39,
                width = 0.14,
                height = 0.11
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.45, "npc"),
            y = unit(0.41, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = HET_colour
            )
        ),

        textGrob(
            "Band + DR",
            x = unit(0.45, "npc"),
            y = unit(0.35, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        # ====================================================
        # Colour grammar
        # ====================================================

        textGrob(
            "Supplement genoplot colours",
            x = unit(0.06, "npc"),
            y = unit(0.20, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.7,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(0.10, "npc"),
            x1 = unit(0.10, "npc"),
            y0 = unit(0.11, "npc"),
            y1 = unit(0.16, "npc"),
            gp = gpar(
                col = DR_colour,
                lwd = 3.0
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.16, "npc"),
            y = unit(0.135, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#555555"
            )
        ),

        segmentsGrob(
            x0 = unit(0.37, "npc"),
            x1 = unit(0.37, "npc"),
            y0 = unit(0.11, "npc"),
            y1 = unit(0.16, "npc"),
            gp = gpar(
                col = HET_colour,
                lwd = 3.0
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.43, "npc"),
            y = unit(0.135, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#555555"
            )
        )
    )
}


add_sidecar <- function(
    plot,
    sidecar_width = unit(
        1.85,
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
        grobs = make_supportive_key(),
        t = panel_t,
        b = panel_b,
        l = panel_l,
        r = panel_l,
        clip = "off",
        name = "supportive-focal-SNP-key"
    )


    g
}


# ============================================================
# 170A
#
# ONLY exact-pattern sites with BA > 0.99 somewhere.
#
# Nothing else in lower panel.
# ============================================================

plot_A <- build_base_plot(

    title_text = paste(
        "B, Dennis and Ray association with",
        "focal Left–Right genotype support"
    ),

    subtitle_text = paste(
        "170A: lower panel shows only SNPs matching",
        "L = DR/DR and R = B/DR that also have BA > 0.99",
        "for B, Dennis, or Ray."
    ),

    caption_text = paste(
        "Lower panel: single magenta L tick = DR/DR;",
        "single orange R tick = B/DR.",
        "B/D/R comparative support is shown in the aligned BA tracks above."
    )
) +


    geom_segment(
        data = high_left,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.95,

        alpha = 1
    ) +


    geom_segment(
        data = high_right_B,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = HET_colour,

        linewidth = 0.92,

        alpha = 1
    ) +


    geom_segment(
        data = high_right_DR,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.92,

        alpha = 0
    )


# ============================================================
# 170B
#
# ALL exact L DR/DR -> R B/DR sites are present faintly.
#
# High-BA subset is stronger / taller.
# ============================================================

plot_B <- build_base_plot(

    title_text = paste(
        "B, Dennis and Ray association with",
        "focal Left–Right genotype support"
    ),

    subtitle_text = paste(
        "170B: all SNPs matching L = DR/DR and R = B/DR",
        "are shown faintly; the subset with BA > 0.99",
        "for B, Dennis, or Ray is emphasized."
    ),

    caption_text = paste(
        "Lower panel: L = DR/DR and R = B/DR is the only",
        "focal genotype pattern treated as supportive.",
        "High-BA sites are taller and fully saturated."
    )
) +


    # --------------------------------------------------------
    # Faint exact-pattern context
    # --------------------------------------------------------

    geom_segment(
        data = all_left,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.28,

        alpha = 0.23
    ) +


    geom_segment(
        data = all_right_B,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = HET_colour,

        linewidth = 0.28,

        alpha = 0.23
    ) +


    geom_segment(
        data = all_right_DR,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.28,

        alpha = 0
    ) +


    # --------------------------------------------------------
    # Strong high-BA subset
    # --------------------------------------------------------

    geom_segment(
        data = high_left,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.95,

        alpha = 1
    ) +


    geom_segment(
        data = high_right_B,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = HET_colour,

        linewidth = 0.92,

        alpha = 1
    ) +


    geom_segment(
        data = high_right_DR,

        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),

        inherit.aes = FALSE,

        colour = DR_colour,

        linewidth = 0.92,

        alpha = 0
    )


# ============================================================
# SIDECARS
# ============================================================

grob_A <- add_sidecar(
    plot_A
)


grob_B <- add_sidecar(
    plot_B
)


# ============================================================
# SAVE INDIVIDUAL VERSIONS
# ============================================================

draw_pdf(
    grob_A,
    pdf_A,
    width = 8.75,
    height = 5.10
)


draw_png(
    grob_A,
    png_A,
    width = 8.75,
    height = 5.10
)


draw_pdf(
    grob_B,
    pdf_B,
    width = 8.75,
    height = 5.10
)


draw_png(
    grob_B,
    png_B,
    width = 8.75,
    height = 5.10
)


# ============================================================
# TWO-PAGE COMPARISON
# ============================================================

cairo_pdf(
    filename = comparison_pdf,
    width = 8.75,
    height = 5.10,
    onefile = TRUE
)


grid.newpage()

grid.draw(
    grob_A
)


grid.newpage()

grid.draw(
    grob_B
)


dev.off()


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n"
)

cat(
    "ANALYSIS 170 COMPLETE\n"
)

cat(
    "========================================\n\n"
)


cat(
    "FOCAL RULE USED:\n"
)

cat(
    "  Left  = DR/DR\n"
)

cat(
    "  Right = B/DR\n\n"
)

cat(
    "No B/B focal state is used.\n"
)

cat(
    "No reverse orientation is accepted.\n\n"
)


print(
    summary_table,
    row.names = FALSE
)


cat(
    "\n170A:\n",
    pdf_A,
    "\n",
    sep = ""
)


cat(
    "\n170B:\n",
    pdf_B,
    "\n",
    sep = ""
)


cat(
    "\nComparison:\n",
    comparison_pdf,
    "\n",
    sep = ""
)


cat(
    "\nAudit:\n",
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
