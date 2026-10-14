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
    "science_BDR_rugs_LR_rows"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

pdf_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_rows_exact_background_highBA_highlight.pdf"
)

png_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_rows_exact_background_highBA_highlight.png"
)

audit_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_rows_exact_background_highBA_highlight.audit.tsv"
)

summary_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_rows_exact_background_highBA_highlight.summary.tsv"
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

deduplicate_by_position <- function(d) {
    if (nrow(d) == 0L) {
        return(d)
    }

    d <- d[
        order(
            d$position_Mb
        ),
        ,
        drop = FALSE
    ]

    d <- d[
        !duplicated(
            paste(
                d$trait,
                round(d$position_Mb * 1e6),
                sep = ":"
            )
        ),
        ,
        drop = FALSE
    ]

    d
}

make_exact_frame <- function(
    x,
    trait_name,
    left_state,
    right_state,
    ba_candidates
) {
    if (nrow(x) == 0L) {
        return(
            data.frame(
                trait = character(),
                position_Mb = numeric(),
                BA = numeric(),
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

    ba_col <- optional_column(
        x,
        ba_candidates
    )

    position_Mb <- suppressWarnings(
        as.numeric(
            x[[pos_col]]
        )
    )

    if (is.na(ba_col)) {
        BA <- rep(
            NA_real_,
            length(position_Mb)
        )
    } else {
        BA <- suppressWarnings(
            as.numeric(
                x[[ba_col]]
            )
        )
    }

    d <- data.frame(
        trait = trait_name,
        position_Mb = position_Mb,
        BA = BA,
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

    deduplicate_by_position(d)
}


# ============================================================
# STATE / COLOUR DEFINITIONS
# ============================================================

state_colours <- c(
    "B/B" = "#3A7CA5",
    "B/DR" = "#7A7A7A",
    "DR/DR" = "#B05A8A"
)

neutral_light <- "#DFE3E7"
neutral_mid <- "#AAB2BA"
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
    B = 5.00,
    D = 4.25,
    R = 3.50
)

ba$rug_base <- unname(
    rug_bases[
        as.character(ba$trait)
    ]
)

ba$tick_height <- ifelse(
    ba$BA > 0.99,
    0.52,
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
# LOWER PANEL DATA
#
# Background = all exact focal-pattern matches
# Highlight  = exact focal-pattern matches with BA > 0.99
#
# B: Left A/A -> Right P/A  ==  Left B/B -> Right B/DR
# D: Left P/P -> Right P/A  ==  Left DR/DR -> Right B/DR
# R: Left P/P -> Right P/A  ==  Left DR/DR -> Right B/DR
# ============================================================

exact_B_all <- make_exact_frame(
    band_correct,
    trait_name = "B",
    left_state = "B/B",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Band_BA",
        "B_BA",
        "band_ba"
    )
)

exact_D_all <- make_exact_frame(
    dr_correct,
    trait_name = "D",
    left_state = "DR/DR",
    right_state = "B/DR",
    ba_candidates = c(
        "Dennis_BA",
        "D_BA",
        "dennis_ba",
        "d_ba"
    )
)

exact_R_all <- make_exact_frame(
    dr_correct,
    trait_name = "R",
    left_state = "DR/DR",
    right_state = "B/DR",
    ba_candidates = c(
        "Ray_BA",
        "R_BA",
        "ray_ba",
        "r_ba"
    )
)

highlight_B <- make_exact_frame(
    high_B,
    trait_name = "B",
    left_state = "B/B",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Band_BA",
        "B_BA",
        "band_ba"
    )
)

highlight_D <- make_exact_frame(
    high_D,
    trait_name = "D",
    left_state = "DR/DR",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Dennis_BA",
        "D_BA",
        "dennis_ba",
        "d_ba"
    )
)

highlight_R <- make_exact_frame(
    high_R,
    trait_name = "R",
    left_state = "DR/DR",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Ray_BA",
        "R_BA",
        "ray_ba",
        "r_ba"
    )
)

exact_all <- rbind(
    exact_B_all,
    exact_D_all,
    exact_R_all
)

highlight_all <- rbind(
    highlight_B,
    highlight_D,
    highlight_R
)

if (nrow(exact_all) == 0L) {
    stop(
        "No exact focal-pattern sites were recovered."
    )
}

highlight_all$key <- paste(
    highlight_all$trait,
    round(highlight_all$position_Mb * 1e6),
    sep = ":"
)

exact_all$key <- paste(
    exact_all$trait,
    round(exact_all$position_Mb * 1e6),
    sep = ":"
)

exact_all$is_high_BA <- (
    exact_all$key %in% highlight_all$key
)


# ============================================================
# BUILD PAIRED L / R ROWS
# ============================================================

row_y <- data.frame(
    trait = c("B", "B", "D", "D", "R", "R"),
    side = c("Left", "Right", "Left", "Right", "Left", "Right"),
    y = c(2.30, 2.02, 1.30, 1.02, 0.30, 0.02),
    stringsAsFactors = FALSE
)

left_points <- merge(
    exact_all[
        ,
        c(
            "trait",
            "position_Mb",
            "BA",
            "left_state",
            "right_state",
            "key",
            "is_high_BA"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Left",
        ,
        drop = FALSE
    ],
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

left_points$state <- left_points$left_state

right_points <- merge(
    exact_all[
        ,
        c(
            "trait",
            "position_Mb",
            "BA",
            "left_state",
            "right_state",
            "key",
            "is_high_BA"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Right",
        ,
        drop = FALSE
    ],
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

right_points$state <- right_points$right_state

lr_points <- rbind(
    left_points,
    right_points
)

lr_points$side <- factor(
    lr_points$side,
    levels = c("Left", "Right")
)

lr_points$state <- factor(
    lr_points$state,
    levels = c("B/B", "B/DR", "DR/DR")
)

highlight_left <- merge(
    highlight_all[
        ,
        c(
            "trait",
            "position_Mb",
            "BA"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Left",
        ,
        drop = FALSE
    ],
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

highlight_left$state <- factor(
    c(
        rep("B/B", nrow(highlight_B)),
        rep("DR/DR", nrow(highlight_D)),
        rep("DR/DR", nrow(highlight_R))
    ),
    levels = c("B/B", "B/DR", "DR/DR")
)

highlight_right <- merge(
    highlight_all[
        ,
        c(
            "trait",
            "position_Mb",
            "BA"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Right",
        ,
        drop = FALSE
    ],
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

highlight_right$state <- factor(
    rep("B/DR", nrow(highlight_all)),
    levels = c("B/B", "B/DR", "DR/DR")
)

highlight_points <- rbind(
    highlight_left,
    highlight_right
)

highlight_points$side <- factor(
    highlight_points$side,
    levels = c("Left", "Right")
)


# connectors
exact_connectors <- merge(
    exact_all[
        ,
        c(
            "trait",
            "position_Mb",
            "BA",
            "key",
            "is_high_BA"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Left",
        c("trait", "y"),
        drop = FALSE
    ],
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

names(exact_connectors)[
    names(exact_connectors) == "y"
] <- "y_left"

tmp_right <- row_y[
    row_y$side == "Right",
    c("trait", "y"),
    drop = FALSE
]

names(tmp_right)[
    names(tmp_right) == "y"
] <- "y_right"

exact_connectors <- merge(
    exact_connectors,
    tmp_right,
    by = "trait",
    all.x = TRUE,
    sort = FALSE
)

highlight_connectors <- exact_connectors[
    exact_connectors$is_high_BA,
    ,
    drop = FALSE
]


# ============================================================
# GUIDE ROWS / LABELS
# ============================================================

trait_headers <- data.frame(
    trait = c("B", "D", "R"),
    x = science_x_limits[1] + 0.0015,
    y = c(2.46, 1.46, 0.46),
    label = c("B", "D", "R"),
    stringsAsFactors = FALSE
)

side_labels <- data.frame(
    trait = c("B", "B", "D", "D", "R", "R"),
    x = science_x_limits[1] + 0.0060,
    y = row_y$y,
    label = c("Left", "Right", "Left", "Right", "Left", "Right"),
    stringsAsFactors = FALSE
)

published_regions <- data.frame(
    region = c("B", "R1", "R2", "D1", "D2"),
    xmin = c(1112651, 1156746, 1174547, 1215210, 1218557) / 1e6,
    xmax = c(1114314, 1160220, 1187833, 1216124, 1221494) / 1e6,
    trait = c("B", "R", "R", "D", "D"),
    stringsAsFactors = FALSE
)

published_regions$y <- 5.88
published_regions$label_y <- 5.98


# ============================================================
# AUDIT / SUMMARY
# ============================================================

audit <- exact_all[
    ,
    c(
        "trait",
        "position_Mb",
        "BA",
        "left_state",
        "right_state",
        "is_high_BA"
    ),
    drop = FALSE
]

audit$pos <- as.integer(
    round(audit$position_Mb * 1e6)
)

audit <- audit[
    order(
        match(audit$trait, c("B", "D", "R")),
        audit$position_Mb
    ),
    c(
        "trait",
        "pos",
        "position_Mb",
        "BA",
        "left_state",
        "right_state",
        "is_high_BA"
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
    trait = c("B", "D", "R"),
    exact_background_sites = c(
        nrow(exact_B_all),
        nrow(exact_D_all),
        nrow(exact_R_all)
    ),
    exact_high_BA_sites = c(
        nrow(highlight_B),
        nrow(highlight_D),
        nrow(highlight_R)
    ),
    left_state = c(
        "B/B",
        "DR/DR",
        "DR/DR"
    ),
    right_state = c(
        "B/DR",
        "B/DR",
        "B/DR"
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
# MAIN PLOT
# ============================================================

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
        linewidth = 2.35,
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
        size = 2.35,
        fontface = "bold",
        vjust = 0,
        show.legend = FALSE
    ) +

    # optix
    geom_segment(
        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 5.48,
            yend = 5.48
        ),
        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 2.9,
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
        y = 5.56,
        label = "optix",
        colour = optix_colour,
        size = 2.45,
        fontface = "italic",
        vjust = 0
    ) +

    # upper rug baselines
    geom_segment(
        data = data.frame(
            trait = c("B", "D", "R"),
            y = c(5.00, 4.25, 3.50),
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
        linewidth = 0.27
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

    # high BA sites highlighted
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

    # B / D / R track labels
    geom_text(
        data = data.frame(
            trait = c("B", "D", "R"),
            x = science_x_limits[1] + 0.0015,
            y = c(5.30, 4.55, 3.80),
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

    # section separator
    geom_segment(
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = 2.72,
            yend = 2.72
        ),
        inherit.aes = FALSE,
        colour = "#D4D7DA",
        linewidth = 0.35
    ) +

    annotate(
        "text",
        x = science_x_limits[1] + 0.0015,
        y = 2.82,
        label = "L vs R exact focal matches",
        hjust = 0,
        vjust = 0,
        size = 2.35,
        fontface = "bold",
        colour = "#505050"
    ) +

    # lower row guides
    geom_segment(
        data = row_y,
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = y,
            yend = y
        ),
        inherit.aes = FALSE,
        colour = neutral_light,
        linewidth = 0.24
    ) +

    # row labels
    geom_text(
        data = trait_headers,
        aes(
            x = x,
            y = y,
            label = label,
            colour = trait
        ),
        hjust = 0,
        fontface = "bold",
        size = 3.0,
        show.legend = FALSE
    ) +

    geom_text(
        data = side_labels,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 0,
        size = 2.10,
        colour = "#5B636B"
    ) +

    # background exact connectors
    geom_segment(
        data = exact_connectors,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left
        ),
        colour = "#BFC5CB",
        linewidth = 0.28,
        alpha = 0.65
    ) +

    # background exact points (state-coloured)
    geom_point(
        data = lr_points,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.75,
        stroke = 0.23,
        colour = "#666666",
        alpha = 0.92
    ) +

    # highlighted exact high-BA connectors
    geom_segment(
        data = highlight_connectors,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left,
            colour = trait
        ),
        linewidth = 0.62,
        alpha = 0.95,
        show.legend = FALSE
    ) +

    # highlighted exact high-BA points
    geom_point(
        data = highlight_points,
        aes(
            x = position_Mb,
            y = y,
            fill = state,
            colour = trait
        ),
        shape = 21,
        size = 2.55,
        stroke = 0.75,
        alpha = 1,
        show.legend = FALSE
    ) +

    scale_fill_manual(
        values = state_colours,
        breaks = c("B/B", "B/DR", "DR/DR"),
        labels = c(
            "B/B  (B only)",
            "B/DR  (both)",
            "DR/DR  (DR only)"
        ),
        name = "Genotype state"
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
            mult = c(0, 0)
        )
    ) +

    coord_cartesian(
        xlim = science_x_limits,
        ylim = c(-0.18, 6.15),
        clip = "off"
    ) +

    labs(
        title = "B, Dennis and Ray association with L vs R exact focal matches",
        subtitle = paste(
            "Upper panel: BA rugs across the region.",
            "Lower panel: all exact focal-pattern matches in the background,",
            "with BA > 0.99 exact matches highlighted."
        ),
        caption = paste(
            "Exact focal patterns:",
            "B = Left A/A -> Right P/A;",
            "D and R = Left P/P -> Right P/A.",
            "Genotype key uses B/B = B only, B/DR = both, DR/DR = DR only."
        )
    ) +

    theme_classic(
        base_size = 8
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 9.6,
            colour = text_dark,
            margin = margin(b = 2)
        ),

        plot.subtitle = element_text(
            size = 7.1,
            colour = "#505050",
            margin = margin(b = 4)
        ),

        plot.caption = element_text(
            size = 6.5,
            colour = "#555555",
            hjust = 0,
            margin = margin(t = 4)
        ),

        axis.title.x = element_text(
            size = 7.5,
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


# ============================================================
# LEFT SIDECAR: CORRECTED GENOTYPE CARTOON KEY
# ============================================================

make_state_sidecar <- function() {
    grobTree(

        textGrob(
            "Genotype key",
            x = unit(0.08, "npc"),
            y = unit(0.96, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8,
                fontface = "bold",
                col = "#303030"
            )
        ),

        # B/B
        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.18,
                y = 0.78,
                width = 0.22,
                height = 0.16
            )
        ),

        textGrob(
            "B/B",
            x = unit(0.38, "npc"),
            y = unit(0.81, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = state_colours["B/B"]
            )
        ),

        textGrob(
            "B only",
            x = unit(0.38, "npc"),
            y = unit(0.74, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.3,
                col = "#4F4F4F"
            )
        ),

        # B/DR = both
        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.14,
                y = 0.50,
                width = 0.16,
                height = 0.12
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.24,
                y = 0.50,
                width = 0.16,
                height = 0.12
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.38, "npc"),
            y = unit(0.53, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = state_colours["B/DR"]
            )
        ),

        textGrob(
            "both",
            x = unit(0.38, "npc"),
            y = unit(0.46, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.3,
                col = "#4F4F4F"
            )
        ),

        # DR/DR
        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.18,
                y = 0.22,
                width = 0.22,
                height = 0.16
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.38, "npc"),
            y = unit(0.25, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = state_colours["DR/DR"]
            )
        ),

        textGrob(
            "DR only",
            x = unit(0.38, "npc"),
            y = unit(0.18, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.3,
                col = "#4F4F4F"
            )
        )
    )
}

add_sidecar <- function(
    plot,
    sidecar_width = unit(1.55, "in")
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
        grobs = make_state_sidecar(),
        t = panel_t,
        b = panel_b,
        l = panel_l,
        r = panel_l,
        clip = "off",
        name = "state-sidecar"
    )

    g
}

final_grob <- add_sidecar(p)


# ============================================================
# SAVE
# ============================================================

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

draw_pdf(
    final_grob,
    pdf_file,
    width = 8.95,
    height = 5.35
)

draw_png(
    final_grob,
    png_file,
    width = 8.95,
    height = 5.35
)


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 164 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nAudit:\n",
    audit_file,
    "\n\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
