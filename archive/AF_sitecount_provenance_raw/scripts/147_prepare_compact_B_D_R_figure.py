#!/usr/bin/env python3

from pathlib import Path

project = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

source = (
    project
    / "scripts"
    / "146_plot_B_D_R_correct_pattern_BA099_ribbons.R"
)

target = (
    project
    / "scripts"
    / "147_plot_B_D_R_compact_BA099_ribbons.R"
)

if not source.exists():
    raise FileNotFoundError(
        f"Missing source script: {source}"
    )

text = source.read_text()


def replace_section(
    text,
    start_heading,
    end_heading,
    replacement,
):
    start = text.index(start_heading)
    end = text.index(
        end_heading,
        start,
    )

    return (
        text[:start]
        + replacement.rstrip()
        + "\n\n"
        + text[end:]
    )


# ============================================================
# OUTPUT NAMES
# ============================================================

text = text.replace(
    "Wallbank_Dennis_Ray_plus_restricted_Band_correct_pattern_BAgt099_ribbons.pdf",
    "Wallbank_Dennis_Ray_plus_restricted_Band_compact_BAgt099_ribbons.pdf",
)

text = text.replace(
    "ANALYSIS 146 FIGURE COMPLETE",
    "ANALYSIS 147 FIGURE COMPLETE",
)


# ============================================================
# INSERT CLEANER UPPER-PANEL STYLING
#
# Comparative BA becomes deliberately neutral.
# The focal genotype ribbons carry the biological colours.
# ============================================================

marker = "# ============================================================\n# BROAD REGION"

style_block = r'''
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


'''

if marker not in text:
    raise RuntimeError(
        "Could not find BROAD REGION marker."
    )

text = text.replace(
    marker,
    style_block + marker,
    1,
)


# ============================================================
# REPLACE THE LOWER TRACKS WITH A COMPACT VERSION
# ============================================================

compact_ribbon_block = r'''
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
'''

text = replace_section(
    text,
    "# SIMPLE CORRECT-PATTERN RIBBON",
    "# ALIGN GROBS",
    compact_ribbon_block,
)


# ============================================================
# COMPACT PAGE SIZE
# ============================================================

old_device = '''grDevices::cairo_pdf(
    filename = output_pdf,
    width = 11.7,
    height = 18.6,
    onefile = TRUE
)'''

new_device = '''grDevices::cairo_pdf(
    filename = output_pdf,
    width = 10.8,
    height = 12.8,
    onefile = TRUE
)'''

if old_device not in text:
    raise RuntimeError(
        "Could not find original PDF device block."
    )

text = text.replace(
    old_device,
    new_device,
)


# Overall A/B balance.
text = text.replace(
    '''c(
                    10.2,
                    8.4
                )''',
    '''c(
                    6.6,
                    6.2
                )''',
)


# Make the focal strips much thinner relative to BA panels.
text = text.replace(
    "upper_fraction = 0.67",
    "upper_fraction = 0.79",
)

text = text.replace(
    "upper_fraction = 0.64",
    "upper_fraction = 0.77",
)


target.write_text(
    text
)

print(
    f"Created: {target}"
)
