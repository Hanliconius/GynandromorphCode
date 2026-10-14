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
    "science_collapsed_LR_difference_emphasis"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

audit_file <- file.path(
    output_dir,
    "Science_collapsed_LR_difference_emphasis.audit.tsv"
)

summary_file <- file.path(
    output_dir,
    "Science_collapsed_LR_difference_emphasis.summary.tsv"
)

comparison_pdf <- file.path(
    output_dir,
    "Science_collapsed_LR_difference_emphasis.comparison.pdf"
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

deduplicate_positions <- function(d, key_cols) {
    if (nrow(d) == 0L) {
        return(d)
    }

    key <- do.call(
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
        !duplicated(key),
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
        round(
            d$position_Mb * 1e6
        )
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

make_highlight_keys <- function(x) {
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

    d <- suppressWarnings(
        as.numeric(
            x[[pos_col]]
        )
    )

    d <- d[
        is.finite(d) &
        d >= science_x_limits[1] &
        d <= science_x_limits[2]
    ]

    as.integer(
        unique(
            round(d * 1e6)
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

state_fill <- c(
    "B/B" = "#3A7CA5",
    "B/DR" = "#B8BDC4",
    "DR/DR" = "#C46A97"
)

state_line <- c(
    "B/B" = "#2B648C",
    "B/DR" = "#7B818A",
    "DR/DR" = "#A45780"
)

support_colours <- c(
    "background" = "#C4C9CE",
    "B" = unname(trait_colours["B"]),
    "D" = unname(trait_colours["D"]),
    "R" = unname(trait_colours["R"]),
    "D+R" = "#8F63B8"
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
    B = 5.08,
    D = 4.31,
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
    high_B
)

highlight_D_keys <- make_highlight_keys(
    high_D
)

highlight_R_keys <- make_highlight_keys(
    high_R
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
        "No collapsed L/R sites were recovered."
    )
}


# ============================================================
# LOWER PANEL GEOMETRY
# ============================================================

left_y <- 1.82
right_y <- 1.26
mid_y <- mean(
    c(
        left_y,
        right_y
    )
)

row_label_x <- science_x_limits[1] + 0.0015
row_note_x <- science_x_limits[1] + 0.0120
mini_offset <- 0.00042


# ============================================================
# POINT / CONNECTOR DATA
# ============================================================

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

midpoints_hi <- data.frame(
    position_Mb = connectors_hi$position_Mb,
    y = mid_y,
    support_class = connectors_hi$support_class,
    stringsAsFactors = FALSE
)

# split connector layers for variant B
split_upper <- data.frame(
    position_Mb = connectors_hi$position_Mb,
    y0 = mid_y,
    y1 = left_y,
    left_state = sites$left_state[
        match(
            connectors_hi$pos,
            sites$pos
        )
    ],
    support_class = connectors_hi$support_class,
    stringsAsFactors = FALSE
)

split_lower <- data.frame(
    position_Mb = connectors_hi$position_Mb,
    y0 = right_y,
    y1 = mid_y,
    right_state = sites$right_state[
        match(
            connectors_hi$pos,
            sites$pos
        )
    ],
    support_class = connectors_hi$support_class,
    stringsAsFactors = FALSE
)

# little two-colour helper marks for highlighted B/DR points
right_hi_bdr <- points_hi[
    points_hi$side == "Right" &
    points_hi$state == "B/DR",
    ,
    drop = FALSE
]

right_hi_bdr_blue <- right_hi_bdr
right_hi_bdr_blue$x2 <- right_hi_bdr_blue$position_Mb - mini_offset

right_hi_bdr_magenta <- right_hi_bdr
right_hi_bdr_magenta$x2 <- right_hi_bdr_magenta$position_Mb + mini_offset


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

published_regions$y <- 5.90
published_regions$label_y <- 6.01


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
# SIDECAR
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
            "Meaning",
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
            "Each vertical pair is one SNP.",
            x = unit(0.08, "npc"),
            y = unit(0.85, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Endpoint colour/fill = focal genotype.",
            x = unit(0.08, "npc"),
            y = unit(0.81, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Grey connector = exact L/R-compatible site.",
            x = unit(0.08, "npc"),
            y = unit(0.77, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Coloured midpoint diamond = BA > 0.99",
            x = unit(0.08, "npc"),
            y = unit(0.73, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "support class (B / D / R).",
            x = unit(0.08, "npc"),
            y = unit(0.69, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.0,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Compatibility classes",
            x = unit(0.06, "npc"),
            y = unit(0.61, "npc"),
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
            y = unit(0.56, "npc"),
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
            y = unit(0.52, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Dennis/Ray-compatible:",
            x = unit(0.08, "npc"),
            y = unit(0.46, "npc"),
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
            y = unit(0.42, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Support colours",
            x = unit(0.06, "npc"),
            y = unit(0.34, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.29, "npc"),
            pch = 23,
            size = unit(2.0, "mm"),
            gp = gpar(
                fill = support_colours["B"],
                col = support_colours["B"]
            )
        ),

        textGrob(
            "blue = B support",
            x = unit(0.17, "npc"),
            y = unit(0.29, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.24, "npc"),
            pch = 23,
            size = unit(2.0, "mm"),
            gp = gpar(
                fill = support_colours["D"],
                col = support_colours["D"]
            )
        ),

        textGrob(
            "orange = Dennis support",
            x = unit(0.17, "npc"),
            y = unit(0.24, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        pointsGrob(
            x = unit(0.10, "npc"),
            y = unit(0.19, "npc"),
            pch = 23,
            size = unit(2.0, "mm"),
            gp = gpar(
                fill = support_colours["R"],
                col = support_colours["R"]
            )
        ),

        textGrob(
            "magenta = Ray support",
            x = unit(0.17, "npc"),
            y = unit(0.19, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Genotype state",
            x = unit(0.06, "npc"),
            y = unit(0.11, "npc"),
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
                y = 0.07,
                width = 0.12,
                height = 0.08
            )
        ),

        textGrob(
            "B/B = B only",
            x = unit(0.21, "npc"),
            y = unit(0.07, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = state_fill["B/B"]
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.10,
                y = 0.035,
                width = 0.08,
                height = 0.05
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.15,
                y = 0.035,
                width = 0.08,
                height = 0.05
            )
        ),

        textGrob(
            "B/DR = both",
            x = unit(0.21, "npc"),
            y = unit(0.035, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = state_line["B/DR"]
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.12,
                y = 0.005,
                width = 0.12,
                height = 0.08
            )
        ),

        textGrob(
            "DR/DR = DR only",
            x = unit(0.21, "npc"),
            y = unit(0.005, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = state_fill["DR/DR"]
            )
        )
    )
}

add_sidecar <- function(
    plot,
    sidecar_width = unit(2.35, "in")
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
        g$layout$l[panel_rows]
    )

    panel_t <- min(
        g$layout$t[panel_rows]
    )

    panel_b <- max(
        g$layout$b[panel_rows]
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
# BASE PLOT
# ============================================================

build_base_plot <- function(
    title_text,
    subtitle_text,
    row_style = c(
        "simple",
        "explicit"
    )
) {
    row_style <- match.arg(
        row_style
    )

    if (row_style == "simple") {
        left_label <- "Left"
        right_label <- "Right"
        left_note <- NULL
        right_note <- NULL
    } else {
        left_label <- "Left"
        right_label <- "Right"
        left_note <- "B/B at Band-compatible sites; DR/DR at Dennis-Ray-compatible sites"
        right_note <- "B/DR at both compatibility classes"
    }

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
            linewidth = 2.28,
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
            size = 2.28,
            fontface = "bold",
            vjust = 0,
            show.legend = FALSE
        ) +

        # optix
        geom_segment(
            aes(
                x = optix_start_science,
                xend = optix_end_science,
                y = 5.57,
                yend = 5.57
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
            y = 5.66,
            label = "optix",
            colour = optix_colour,
            size = 2.42,
            fontface = "italic",
            vjust = 0
        ) +

        # upper baselines
        geom_segment(
            data = data.frame(
                trait = c("B", "D", "R"),
                y = c(5.08, 4.31, 3.54),
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
            linewidth = 0.74,
            alpha = 1,
            lineend = "butt",
            show.legend = FALSE
        ) +

        # trait labels
        geom_text(
            data = data.frame(
                trait = c("B", "D", "R"),
                x = science_x_limits[1] + 0.0015,
                y = c(5.33, 4.56, 3.79),
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
                y = 2.60,
                yend = 2.60
            ),
            inherit.aes = FALSE,
            colour = "#D4D7DA",
            linewidth = 0.35
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = 2.70,
            label = "Collapsed focal Left vs Right panel",
            hjust = 0,
            vjust = 0,
            size = 2.36,
            fontface = "bold",
            colour = "#505050"
        ) +

        annotate(
            "text",
            x = science_x_limits[2] - 0.0012,
            y = 2.70,
            label = "endpoint = genotype; midpoint diamond = comparative support",
            hjust = 1,
            vjust = 0,
            size = 2.0,
            colour = "#666666"
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
            linewidth = 0.27
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
            linewidth = 0.27
        ) +

        annotate(
            "text",
            x = row_label_x,
            y = left_y + 0.12,
            label = left_label,
            hjust = 0,
            fontface = "bold",
            size = 3.0,
            colour = "#444444"
        ) +

        annotate(
            "text",
            x = row_label_x,
            y = right_y + 0.12,
            label = right_label,
            hjust = 0,
            fontface = "bold",
            size = 3.0,
            colour = "#444444"
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
            ylim = c(1.02, 6.20),
            clip = "off"
        ) +

        labs(
            title = title_text,
            subtitle = subtitle_text,
            caption = paste(
                "Grey marks show all exact-compatible focal sites.",
                "Coloured midpoint diamonds show BA > 0.99 support for B, D, or R.",
                "Point colour/fill shows the focal genotype state."
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

    if (!is.null(left_note)) {
        p <- p +
            annotate(
                "text",
                x = row_note_x,
                y = left_y + 0.01,
                label = left_note,
                hjust = 0,
                size = 1.95,
                colour = "#666666"
            )
    }

    if (!is.null(right_note)) {
        p <- p +
            annotate(
                "text",
                x = row_note_x,
                y = right_y + 0.01,
                label = right_note,
                hjust = 0,
                size = 1.95,
                colour = "#666666"
            )
    }

    p
}


# ============================================================
# VARIANT 167A: MIDPOINT-SUPPORT
# ============================================================

plot_167A <- build_base_plot(
    title_text = "Collapsed Left vs Right focal panel",
    subtitle_text = paste(
        "167A: focal genotype difference is shown at the Left and Right endpoints;",
        "comparative B / D / R support is shown separately at the midpoint."
    ),
    row_style = "simple"
) +

    # background connectors
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

    # background endpoints
    geom_point(
        data = points_bg,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.65,
        stroke = 0.20,
        colour = "#6F747C",
        alpha = 0.58
    ) +

    # highlighted endpoints
    geom_point(
        data = points_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 2.55,
        stroke = 0.48,
        colour = "#1F1F1F",
        alpha = 1
    ) +

    # tiny blue+magenta flanking marks for highlighted B/DR right-side sites
    geom_point(
        data = right_hi_bdr_blue,
        aes(
            x = x2,
            y = y
        ),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.80,
        colour = state_fill["B/B"],
        alpha = 0.95
    ) +

    geom_point(
        data = right_hi_bdr_magenta,
        aes(
            x = x2,
            y = y
        ),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.80,
        colour = state_fill["DR/DR"],
        alpha = 0.95
    ) +

    # midpoint support diamonds
    geom_point(
        data = midpoints_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = support_class,
            colour = support_class
        ),
        shape = 23,
        size = 2.15,
        stroke = 0.72,
        alpha = 1,
        show.legend = FALSE
    ) +

    scale_fill_manual(
        values = c(
            state_fill,
            support_colours
        ),
        drop = FALSE
    ) +

    scale_colour_manual(
        values = c(
            state_line,
            support_colours,
            trait_colours
        ),
        drop = FALSE
    )


# ============================================================
# VARIANT 167B: SPLIT-CONNECTOR + MIDPOINT-SUPPORT
# ============================================================

plot_167B <- build_base_plot(
    title_text = "Collapsed Left vs Right focal panel",
    subtitle_text = paste(
        "167B: as in 167A, but highlighted sites also show a split connector",
        "to emphasize the genotype change from Left to Right."
    ),
    row_style = "explicit"
) +

    # background connectors
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
        alpha = 0.70
    ) +

    # background endpoints
    geom_point(
        data = points_bg,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.58,
        stroke = 0.20,
        colour = "#717780",
        alpha = 0.54
    ) +

    # split upper connector = left genotype state
    geom_segment(
        data = split_upper,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1,
            colour = left_state
        ),
        linewidth = 0.72,
        alpha = 0.96,
        show.legend = FALSE
    ) +

    # split lower connector = right genotype state
    geom_segment(
        data = split_lower,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1,
            colour = right_state
        ),
        linewidth = 0.72,
        alpha = 0.96,
        show.legend = FALSE
    ) +

    # highlighted endpoints
    geom_point(
        data = points_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = state,
            colour = state
        ),
        shape = 21,
        size = 2.55,
        stroke = 0.60,
        alpha = 1,
        show.legend = FALSE
    ) +

    # tiny blue+magenta flanking marks for highlighted B/DR right-side sites
    geom_point(
        data = right_hi_bdr_blue,
        aes(
            x = x2,
            y = y
        ),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.90,
        colour = state_fill["B/B"],
        alpha = 0.98
    ) +

    geom_point(
        data = right_hi_bdr_magenta,
        aes(
            x = x2,
            y = y
        ),
        inherit.aes = FALSE,
        shape = 16,
        size = 0.90,
        colour = state_fill["DR/DR"],
        alpha = 0.98
    ) +

    # midpoint support diamonds
    geom_point(
        data = midpoints_hi,
        aes(
            x = position_Mb,
            y = y,
            fill = support_class,
            colour = support_class
        ),
        shape = 23,
        size = 2.20,
        stroke = 0.78,
        alpha = 1,
        show.legend = FALSE
    ) +

    scale_fill_manual(
        values = c(
            state_fill,
            support_colours
        ),
        drop = FALSE
    ) +

    scale_colour_manual(
        values = c(
            state_line,
            support_colours,
            trait_colours
        ),
        drop = FALSE
    )


# ============================================================
# SAVE
# ============================================================

plots <- list(
    variant_167A_midpoint = plot_167A,
    variant_167B_split = plot_167B
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
        width = 9.35,
        height = 5.35
    )

    draw_png(
        final_grobs[[nm]],
        png_file,
        width = 9.35,
        height = 5.35
    )
}

cairo_pdf(
    filename = comparison_pdf,
    width = 9.35,
    height = 5.35,
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
    "ANALYSIS 167 COMPLETE\n",
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
                "Science_variant_167A_midpoint.pdf",
                "Science_variant_167A_midpoint.png",
                "Science_variant_167B_split.pdf",
                "Science_variant_167B_split.png"
            )
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)
