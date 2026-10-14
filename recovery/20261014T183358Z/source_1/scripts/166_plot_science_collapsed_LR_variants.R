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
# PATHS
# ============================================================

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

source_file <- file.path(
    project,
    "scripts",
    "159_plot_science_trials_pubregions_cartoons.R"
)

if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing required source:",
            source_file
        )
    )
}


# ============================================================
# LOAD ESTABLISHED OBJECTS
# ============================================================

suppressWarnings(
    source(
        source_file,
        local = TRUE
    )
)

required_objects <- c(
    "all_ba_science",
    "band_correct",
    "dr_correct",
    "high_B",
    "high_D",
    "high_R",
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
    "science_collapsed_LR_variants"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

audit_file <- file.path(
    output_dir,
    "Science_collapsed_LR_variants.audit.tsv"
)

summary_file <- file.path(
    output_dir,
    "Science_collapsed_LR_variants.summary.tsv"
)

comparison_pdf <- file.path(
    output_dir,
    "Science_collapsed_LR_variants.comparison.pdf"
)


# ============================================================
# HELPERS
# ============================================================

normalize_name <- function(x) {
    x <- tolower(
        as.character(x)
    )

    gsub(
        "[^a-z0-9]+",
        "",
        x,
        perl = TRUE
    )
}

resolve_column <- function(x, candidates) {
    normalized_columns <- normalize_name(
        names(x)
    )

    for (candidate in candidates) {
        target <- normalize_name(
            candidate
        )

        hit <- which(
            normalized_columns == target
        )

        if (length(hit) == 1L) {
            return(
                names(x)[hit]
            )
        }
    }

    stop(
        paste0(
            "Could not resolve column from: ",
            paste(candidates, collapse = ", "),
            "\nAvailable columns:\n",
            paste(names(x), collapse = "\n")
        )
    )
}

optional_column <- function(x, candidates) {
    normalized_columns <- normalize_name(
        names(x)
    )

    for (candidate in candidates) {
        target <- normalize_name(
            candidate
        )

        hit <- which(
            normalized_columns == target
        )

        if (length(hit) == 1L) {
            return(
                names(x)[hit]
            )
        }
    }

    NA_character_
}

deduplicate_positions <- function(d, key_cols) {
    if (nrow(d) == 0L) {
        return(d)
    }

    keys <- do.call(
        paste,
        c(
            d[
                ,
                key_cols,
                drop = FALSE
            ],
            sep = ":"
        )
    )

    d[
        !duplicated(keys),
        ,
        drop = FALSE
    ]
}

make_exact_frame <- function(
    x,
    compatibility_class,
    left_state,
    right_state
) {
    if (nrow(x) == 0L) {
        return(
            data.frame(
                compatibility_class = character(),
                position_Mb = numeric(),
                left_state = character(),
                right_state = character(),
                stringsAsFactors = FALSE
            )
        )
    }

    pos_col <- resolve_column(
        x,
        c(
            "position_Mb",
            "position_mb",
            "pos_Mb"
        )
    )

    d <- data.frame(
        compatibility_class = compatibility_class,
        position_Mb = suppressWarnings(
            as.numeric(
                x[[pos_col]]
            )
        ),
        left_state = left_state,
        right_state = right_state,
        stringsAsFactors = FALSE
    )

    d <- d[
        is.finite(d$position_Mb) &
        d$position_Mb >= science_x_limits[1] &
        d$position_Mb <= science_x_limits[2],
        ,
        drop = FALSE
    ]

    d$pos <- as.integer(
        round(d$position_Mb * 1e6)
    )

    d <- d[
        order(d$position_Mb),
        ,
        drop = FALSE
    ]

    deduplicate_positions(
        d,
        key_cols = c(
            "compatibility_class",
            "pos"
        )
    )
}

make_highlight_keys <- function(
    x,
    compatibility_class
) {
    if (nrow(x) == 0L) {
        return(
            integer(0)
        )
    }

    pos_col <- resolve_column(
        x,
        c(
            "position_Mb",
            "position_mb",
            "pos_Mb"
        )
    )

    d <- data.frame(
        compatibility_class = compatibility_class,
        position_Mb = suppressWarnings(
            as.numeric(
                x[[pos_col]]
            )
        ),
        stringsAsFactors = FALSE
    )

    d <- d[
        is.finite(d$position_Mb) &
        d$position_Mb >= science_x_limits[1] &
        d$position_Mb <= science_x_limits[2],
        ,
        drop = FALSE
    ]

    as.integer(
        unique(
            round(d$position_Mb * 1e6)
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
    grid.draw(grob)
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
    grid.draw(grob)
    dev.off()
}


# ============================================================
# COLOURS
# ============================================================

state_colours <- c(
    "B/B" = "#3A7CA5",
    "B/DR" = "#8B8B8B",
    "DR/DR" = "#C46A97"
)

support_colours <- c(
    "background" = "#C4C9CE",
    "B" = unname(trait_colours["B"]),
    "D" = unname(trait_colours["D"]),
    "R" = unname(trait_colours["R"]),
    "D+R" = "#8F63B8"
)

combined_colour_scale <- c(
    trait_colours,
    "background" = support_colours["background"],
    "D+R" = support_colours["D+R"]
)

neutral_light <- "#E0E4E8"
neutral_mid <- "#AEB6BE"
neutral_dark <- "#66707A"
text_dark <- "#242424"
optix_colour <- "#3A7CA5"


# ============================================================
# UPPER PANEL DATA
# ============================================================

ba_trait_col <- resolve_column(
    all_ba_science,
    c("trait")
)

ba_pos_col <- resolve_column(
    all_ba_science,
    c("position_Mb", "position_mb")
)

ba_value_col <- resolve_column(
    all_ba_science,
    c("BA", "balanced_accuracy")
)

ba <- data.frame(
    trait = as.character(
        all_ba_science[[ba_trait_col]]
    ),
    position_Mb = suppressWarnings(
        as.numeric(
            all_ba_science[[ba_pos_col]]
        )
    ),
    BA = suppressWarnings(
        as.numeric(
            all_ba_science[[ba_value_col]]
        )
    ),
    stringsAsFactors = FALSE
)

ba <- ba[
    is.finite(ba$position_Mb) &
    is.finite(ba$BA) &
    ba$position_Mb >= science_x_limits[1] &
    ba$position_Mb <= science_x_limits[2],
    ,
    drop = FALSE
]

ba$trait <- factor(
    ba$trait,
    levels = c("B", "D", "R"),
    ordered = TRUE
)

rug_bases <- c(
    B = 5.10,
    D = 4.32,
    R = 3.54
)

ba$rug_base <- unname(
    rug_bases[
        as.character(ba$trait)
    ]
)

ba$tick_height <- ifelse(
    ba$BA > 0.99,
    0.54,
    ifelse(
        ba$BA >= 0.90,
        0.38,
        ifelse(
            ba$BA >= 0.80,
            0.24,
            0.09
        )
    )
)

ba$rug_top <- ba$rug_base + ba$tick_height

high_all <- ba[
    ba$BA > 0.99,
    ,
    drop = FALSE
]


# ============================================================
# LOWER PANEL CORE DATA
#
# Collapsed across compatibility classes into one Left row
# and one Right row.
#
# Band-compatible:
#   Left  = B/B
#   Right = B/DR
#
# Dennis/Ray-compatible:
#   Left  = DR/DR
#   Right = B/DR
# ============================================================

band_bg <- make_exact_frame(
    band_correct,
    compatibility_class = "Band-compatible",
    left_state = "B/B",
    right_state = "B/DR"
)

dr_bg <- make_exact_frame(
    dr_correct,
    compatibility_class = "Dennis-Ray-compatible",
    left_state = "DR/DR",
    right_state = "B/DR"
)

highlight_B_keys <- make_highlight_keys(
    high_B,
    "Band-compatible"
)

highlight_D_keys <- make_highlight_keys(
    high_D,
    "Dennis-Ray-compatible"
)

highlight_R_keys <- make_highlight_keys(
    high_R,
    "Dennis-Ray-compatible"
)

band_bg$support_class <- ifelse(
    band_bg$pos %in% highlight_B_keys,
    "B",
    "background"
)

dr_bg$support_class <- ifelse(
    dr_bg$pos %in% highlight_D_keys &
    dr_bg$pos %in% highlight_R_keys,
    "D+R",
    ifelse(
        dr_bg$pos %in% highlight_D_keys,
        "D",
        ifelse(
            dr_bg$pos %in% highlight_R_keys,
            "R",
            "background"
        )
    )
)

sites <- rbind(
    band_bg,
    dr_bg
)

sites$compatibility_class <- factor(
    sites$compatibility_class,
    levels = c(
        "Band-compatible",
        "Dennis-Ray-compatible"
    )
)

sites$support_class <- factor(
    sites$support_class,
    levels = c(
        "background",
        "B",
        "D",
        "R",
        "D+R"
    )
)

sites <- sites[
    order(
        sites$position_Mb,
        sites$compatibility_class
    ),
    ,
    drop = FALSE
]

if (nrow(sites) == 0L) {
    stop(
        "No collapsed Left/Right sites were recovered."
    )
}


# ============================================================
# POINT / CONNECTOR DATA FOR COMMON L/R ROWS
# ============================================================

left_y <- 1.78
right_y <- 1.28

points_bg <- rbind(
    data.frame(
        side = "Left",
        y = left_y,
        position_Mb = sites$position_Mb,
        pos = sites$pos,
        compatibility_class = sites$compatibility_class,
        state = sites$left_state,
        support_class = sites$support_class,
        stringsAsFactors = FALSE
    ),
    data.frame(
        side = "Right",
        y = right_y,
        position_Mb = sites$position_Mb,
        pos = sites$pos,
        compatibility_class = sites$compatibility_class,
        state = sites$right_state,
        support_class = sites$support_class,
        stringsAsFactors = FALSE
    )
)

points_bg$state <- factor(
    points_bg$state,
    levels = c(
        "B/B",
        "B/DR",
        "DR/DR"
    )
)

points_bg$side <- factor(
    points_bg$side,
    levels = c(
        "Left",
        "Right"
    )
)

points_hi <- points_bg[
    points_bg$support_class != "background",
    ,
    drop = FALSE
]

connectors_bg <- data.frame(
    position_Mb = sites$position_Mb,
    pos = sites$pos,
    compatibility_class = sites$compatibility_class,
    support_class = sites$support_class,
    y_left = left_y,
    y_right = right_y,
    stringsAsFactors = FALSE
)

connectors_hi <- connectors_bg[
    connectors_bg$support_class != "background",
    ,
    drop = FALSE
]


# ============================================================
# PUBLISHED REGIONS
# ============================================================

published_regions <- data.frame(
    region = c("B", "R1", "R2", "D1", "D2"),
    xmin = c(1112651, 1156746, 1174547, 1215210, 1218557) / 1e6,
    xmax = c(1114314, 1160220, 1187833, 1216124, 1221494) / 1e6,
    trait = c("B", "R", "R", "D", "D"),
    stringsAsFactors = FALSE
)

published_regions$y <- 5.92
published_regions$label_y <- 6.03


# ============================================================
# SUMMARY / AUDIT
# ============================================================

DR_overlap_n <- sum(
    highlight_D_keys %in% highlight_R_keys
)

summary_table <- data.frame(
    quantity = c(
        "Band-compatible background sites",
        "Dennis-Ray-compatible background sites",
        "B_high_BA_sites",
        "D_high_BA_sites",
        "R_high_BA_sites",
        "D_R_overlap_sites"
    ),
    value = c(
        nrow(band_bg),
        nrow(dr_bg),
        length(highlight_B_keys),
        length(highlight_D_keys),
        length(highlight_R_keys),
        DR_overlap_n
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

audit <- sites[
    ,
    c(
        "compatibility_class",
        "pos",
        "position_Mb",
        "left_state",
        "right_state",
        "support_class"
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


# ============================================================
# LOWER-PANEL GUIDE SIDECAR
# ============================================================

make_sidecar <- function() {
    grobTree(

        rectGrob(
            x = 0.5,
            y = 0.5,
            width = 0.98,
            height = 0.98,
            gp = gpar(
                fill = "white",
                col = NA
            )
        ),

        textGrob(
            "Lower panel",
            x = unit(0.06, "npc"),
            y = unit(0.965, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8.0,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "What is shown",
            x = unit(0.06, "npc"),
            y = unit(0.90, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "Each vertical pair is one SNP",
            x = unit(0.08, "npc"),
            y = unit(0.85, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "comparing Left vs Right.",
            x = unit(0.08, "npc"),
            y = unit(0.81, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Compatibility classes",
            x = unit(0.06, "npc"),
            y = unit(0.73, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "Band-compatible:",
            x = unit(0.08, "npc"),
            y = unit(0.68, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "Left = B/B, Right = B/DR",
            x = unit(0.14, "npc"),
            y = unit(0.64, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Dennis/Ray-compatible:",
            x = unit(0.08, "npc"),
            y = unit(0.58, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "Left = DR/DR, Right = B/DR",
            x = unit(0.14, "npc"),
            y = unit(0.54, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Support colours",
            x = unit(0.06, "npc"),
            y = unit(0.46, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.41, "npc"),
            pch = 16,
            size = unit(1.8, "mm"),
            gp = gpar(
                col = support_colours["background"]
            )
        ),

        textGrob(
            "grey = all exact L/R-compatible sites",
            x = unit(0.17, "npc"),
            y = unit(0.41, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#4F4F4F"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.36, "npc"),
            pch = 16,
            size = unit(1.8, "mm"),
            gp = gpar(
                col = support_colours["B"]
            )
        ),

        textGrob(
            "blue = B high-BA support",
            x = unit(0.17, "npc"),
            y = unit(0.36, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#4F4F4F"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.31, "npc"),
            pch = 16,
            size = unit(1.8, "mm"),
            gp = gpar(
                col = support_colours["D"]
            )
        ),

        textGrob(
            "orange = Dennis high-BA support",
            x = unit(0.17, "npc"),
            y = unit(0.31, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#4F4F4F"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.26, "npc"),
            pch = 16,
            size = unit(1.8, "mm"),
            gp = gpar(
                col = support_colours["R"]
            )
        ),

        textGrob(
            "magenta = Ray high-BA support",
            x = unit(0.17, "npc"),
            y = unit(0.26, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Genotype state",
            x = unit(0.06, "npc"),
            y = unit(0.18, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.12,
                y = 0.12,
                width = 0.13,
                height = 0.09
            )
        ),

        textGrob(
            "B/B = B only",
            x = unit(0.21, "npc"),
            y = unit(0.12, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = state_colours["B/B"]
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.10,
                y = 0.07,
                width = 0.08,
                height = 0.06
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.15,
                y = 0.07,
                width = 0.08,
                height = 0.06
            )
        ),

        textGrob(
            "B/DR = both",
            x = unit(0.21, "npc"),
            y = unit(0.07, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = state_colours["B/DR"]
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.12,
                y = 0.02,
                width = 0.13,
                height = 0.09
            )
        ),

        textGrob(
            "DR/DR = DR only",
            x = unit(0.21, "npc"),
            y = unit(0.02, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = state_colours["DR/DR"]
            )
        )
    )
}

add_sidecar <- function(
    plot,
    sidecar_width = unit(2.20, "in")
) {
    g <- ggplotGrob(plot)

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
        grobs = make_sidecar(),
        t = panel_t,
        b = panel_b,
        l = panel_l,
        r = panel_l,
        clip = "off",
        name = "collapsed-LR-sidecar"
    )

    g
}


# ============================================================
# BASE PLOT FUNCTION
# ============================================================

build_base_plot <- function(
    variant_title,
    variant_subtitle,
    inline_note = NULL
) {
    p <- ggplot() +

        # published intervals
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
                x = (xmin + xmax) / 2,
                y = label_y,
                label = region,
                colour = trait
            ),
            size = 2.30,
            fontface = "bold",
            vjust = 0,
            show.legend = FALSE
        ) +

        # optix
        geom_segment(
            aes(
                x = optix_start_science,
                xend = optix_end_science,
                y = 5.58,
                yend = 5.58
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
            y = 5.67,
            label = "optix",
            colour = optix_colour,
            size = 2.45,
            fontface = "italic",
            vjust = 0
        ) +

        # upper baselines
        geom_segment(
            data = data.frame(
                trait = c("B", "D", "R"),
                y = c(5.10, 4.32, 3.54),
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

        # all BA sites
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

        # BA > 0.99 sites
        geom_segment(
            data = high_all,
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

        # trait labels
        geom_text(
            data = data.frame(
                trait = c("B", "D", "R"),
                x = science_x_limits[1] + 0.0015,
                y = c(5.34, 4.56, 3.78),
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

        # separator
        geom_segment(
            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = 2.62,
                yend = 2.62
            ),
            inherit.aes = FALSE,
            colour = "#D4D7DA",
            linewidth = 0.35
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = 2.72,
            label = "Collapsed focal Left vs Right panel",
            hjust = 0,
            vjust = 0,
            size = 2.35,
            fontface = "bold",
            colour = "#505050"
        ) +

        # lower guide rows
        geom_segment(
            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = left_y,
                yend = left_y
            ),
            inherit.aes = FALSE,
            colour = neutral_light,
            linewidth = 0.26
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
            linewidth = 0.26
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = left_y + 0.12,
            label = "Left",
            hjust = 0,
            fontface = "bold",
            size = 3.0,
            colour = "#444444"
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = right_y + 0.12,
            label = "Right",
            hjust = 0,
            fontface = "bold",
            size = 3.0,
            colour = "#444444"
        ) +

        scale_fill_manual(
            values = state_colours,
            drop = FALSE
        ) +

        scale_colour_manual(
            values = combined_colour_scale,
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
                mult = c(0, 0)
            )
        ) +

        coord_cartesian(
            xlim = science_x_limits,
            ylim = c(1.05, 6.22),
            clip = "off"
        ) +

        labs(
            title = variant_title,
            subtitle = variant_subtitle,
            caption = paste(
                "Point fill = genotype state.",
                "Grey marks = all exact focal-compatible sites.",
                "Colour marks high-BA support for B, D, or R."
            )
        ) +

        theme_classic(
            base_size = 8
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 9.5,
                colour = text_dark,
                margin = margin(b = 2)
            ),
            plot.subtitle = element_text(
                size = 7.0,
                colour = "#505050",
                margin = margin(b = 4)
            ),
            plot.caption = element_text(
                size = 6.3,
                colour = "#555555",
                hjust = 0,
                margin = margin(t = 4)
            ),
            axis.title.x = element_text(
                size = 7.4,
                margin = margin(t = 4)
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

    if (!is.null(inline_note)) {
        p <- p +
            annotate(
                "text",
                x = science_x_limits[2] - 0.001,
                y = 2.70,
                label = inline_note,
                hjust = 1,
                vjust = 0,
                size = 2.0,
                colour = "#666666"
            )
    }

    p
}


# ============================================================
# VARIANT 1: CONNECTOR EMPHASIS
# ============================================================

plot_variant_1 <- build_base_plot(
    variant_title = "Collapsed Left vs Right focal panel",
    variant_subtitle = paste(
        "Variant 1: high-BA support is emphasized mainly by coloured connectors;"
        ,
        "genotype state is shown by point fill."
    )
) +

    geom_segment(
        data = connectors_bg,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left
        ),
        colour = support_colours["background"],
        linewidth = 0.24,
        alpha = 0.72
    ) +

    geom_point(
        data = points_bg,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.55,
        stroke = 0.20,
        colour = "#666666",
        alpha = 0.93
    ) +

    geom_segment(
        data = connectors_hi,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left,
            colour = support_class
        ),
        linewidth = 0.72,
        alpha = 0.97,
        show.legend = FALSE
    ) +

    geom_point(
        data = points_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 2.15,
        stroke = 0.25,
        colour = "#333333",
        alpha = 1
    )


# ============================================================
# VARIANT 2: OUTLINE EMPHASIS
# ============================================================

plot_variant_2 <- build_base_plot(
    variant_title = "Collapsed Left vs Right focal panel",
    variant_subtitle = paste(
        "Variant 2: high-BA support is emphasized mainly by coloured point outlines;"
        ,
        "connectors remain grey."
    )
) +

    geom_segment(
        data = connectors_bg,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left
        ),
        colour = support_colours["background"],
        linewidth = 0.26,
        alpha = 0.76
    ) +

    geom_point(
        data = points_bg,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.65,
        stroke = 0.22,
        colour = "#777777",
        alpha = 0.95
    ) +

    geom_point(
        data = points_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = state,
            colour = support_class
        ),
        shape = 21,
        size = 2.45,
        stroke = 0.85,
        alpha = 1,
        show.legend = FALSE
    )


# ============================================================
# VARIANT 3: DUAL EMPHASIS + INLINE NOTE
# ============================================================

plot_variant_3 <- build_base_plot(
    variant_title = "Collapsed Left vs Right focal panel",
    variant_subtitle = paste(
        "Variant 3: high-BA support is emphasized by both coloured connectors"
        ,
        "and coloured point outlines."
    ),
    inline_note = "grey = all exact-compatible sites; colour = BA > 0.99"
) +

    geom_segment(
        data = connectors_bg,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left
        ),
        colour = support_colours["background"],
        linewidth = 0.24,
        alpha = 0.72
    ) +

    geom_point(
        data = points_bg,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.55,
        stroke = 0.20,
        colour = "#6C6C6C",
        alpha = 0.93
    ) +

    geom_segment(
        data = connectors_hi,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left,
            colour = support_class
        ),
        linewidth = 0.70,
        alpha = 0.96,
        show.legend = FALSE
    ) +

    geom_point(
        data = points_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = state,
            colour = support_class
        ),
        shape = 21,
        size = 2.40,
        stroke = 0.82,
        alpha = 1,
        show.legend = FALSE
    )


# ============================================================
# ADD SIDECAR AND SAVE
# ============================================================

plots <- list(
    variant_1_connector = plot_variant_1,
    variant_2_outline = plot_variant_2,
    variant_3_dual = plot_variant_3
)

final_grobs <- lapply(
    plots,
    add_sidecar
)

for (nm in names(final_grobs)) {
    pdf_file <- file.path(
        output_dir,
        paste0(
            "Science_",
            nm,
            ".pdf"
        )
    )

    png_file <- file.path(
        output_dir,
        paste0(
            "Science_",
            nm,
            ".png"
        )
    )

    draw_pdf(
        final_grobs[[nm]],
        pdf_file,
        width = 9.20,
        height = 5.30
    )

    draw_png(
        final_grobs[[nm]],
        png_file,
        width = 9.20,
        height = 5.30
    )
}

cairo_pdf(
    filename = comparison_pdf,
    width = 9.20,
    height = 5.30,
    onefile = TRUE
)

for (nm in names(final_grobs)) {
    grid.newpage()
    grid.draw(
        final_grobs[[nm]]
    )
}

dev.off()


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 166 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nOutput directory:\n",
    output_dir,
    "\n",
    sep = ""
)

cat(
    "\nComparison PDF:\n",
    comparison_pdf,
    "\n",
    sep = ""
)

cat(
    "\nAudit table:\n",
    audit_file,
    "\n",
    sep = ""
)

cat(
    "\nSummary table:\n",
    summary_file,
    "\n",
    sep = ""
)

cat(
    "\nIndividual files:\n",
    paste(
        file.path(
            output_dir,
            c(
                "Science_variant_1_connector.pdf",
                "Science_variant_1_connector.png",
                "Science_variant_2_outline.pdf",
                "Science_variant_2_outline.png",
                "Science_variant_3_dual.pdf",
                "Science_variant_3_dual.png"
            )
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)
