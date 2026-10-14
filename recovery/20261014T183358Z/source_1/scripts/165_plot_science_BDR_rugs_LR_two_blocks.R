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
    "science_BDR_rugs_LR_two_blocks"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

pdf_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_two_blocks_exact_background_highBA_highlight.pdf"
)

png_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_two_blocks_exact_background_highBA_highlight.png"
)

audit_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_two_blocks_exact_background_highBA_highlight.audit.tsv"
)

summary_file <- file.path(
    output_dir,
    "Science_BDR_BA_rugs_LR_two_blocks_exact_background_highBA_highlight.summary.tsv"
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
                d$block,
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
    block_name,
    left_state,
    right_state,
    ba_candidates = character()
) {
    if (nrow(x) == 0L) {
        return(
            data.frame(
                block = character(),
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
        block = block_name,
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

positions_to_key <- function(x, prefix) {
    paste(
        prefix,
        round(x * 1e6),
        sep = ":"
    )
}


# ============================================================
# COLOURS
# ============================================================

state_colours <- c(
    "B/B" = "#3A7CA5",
    "B/DR" = "#8A8A8A",
    "DR/DR" = "#C46A97"
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
    B = 5.15,
    D = 4.35,
    R = 3.55
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
        0.39,
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
# Two L-vs-R blocks:
#
# 1) Band-compatible
#      Left  = B/B
#      Right = B/DR
#
# 2) Dennis/Ray-compatible
#      Left  = DR/DR
#      Right = B/DR
#
# Background = all exact focal-compatible sites
# Highlight  = trait-specific BA > 0.99 subsets
# ============================================================

band_bg <- make_exact_frame(
    band_correct,
    block_name = "Band-compatible",
    left_state = "B/B",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Band_BA",
        "B_BA",
        "band_ba"
    )
)

dr_bg <- make_exact_frame(
    dr_correct,
    block_name = "Dennis-Ray-compatible",
    left_state = "DR/DR",
    right_state = "B/DR",
    ba_candidates = c(
        "BA",
        "Dennis_BA",
        "Ray_BA",
        "D_BA",
        "R_BA",
        "dennis_ba",
        "ray_ba"
    )
)

highlight_B <- make_exact_frame(
    high_B,
    block_name = "Band-compatible",
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
    block_name = "Dennis-Ray-compatible",
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
    block_name = "Dennis-Ray-compatible",
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

if (nrow(band_bg) == 0L) {
    stop(
        "No exact Band-compatible focal sites found."
    )
}

if (nrow(dr_bg) == 0L) {
    stop(
        "No exact Dennis/Ray-compatible focal sites found."
    )
}

band_bg$highlight_B <- positions_to_key(
    band_bg$position_Mb,
    "Band-compatible"
) %in% positions_to_key(
    highlight_B$position_Mb,
    "Band-compatible"
)

band_bg$highlight_D <- FALSE
band_bg$highlight_R <- FALSE

dr_bg$highlight_B <- FALSE
dr_bg$highlight_D <- positions_to_key(
    dr_bg$position_Mb,
    "Dennis-Ray-compatible"
) %in% positions_to_key(
    highlight_D$position_Mb,
    "Dennis-Ray-compatible"
)

dr_bg$highlight_R <- positions_to_key(
    dr_bg$position_Mb,
    "Dennis-Ray-compatible"
) %in% positions_to_key(
    highlight_R$position_Mb,
    "Dennis-Ray-compatible"
)

background_sites <- rbind(
    band_bg,
    dr_bg
)

background_sites$pos <- as.integer(
    round(
        background_sites$position_Mb * 1e6
    )
)


# ------------------------------------------------------------
# Row geometry
# ------------------------------------------------------------

row_y <- data.frame(
    block = c(
        "Band-compatible",
        "Band-compatible",
        "Dennis-Ray-compatible",
        "Dennis-Ray-compatible"
    ),
    side = c(
        "Left",
        "Right",
        "Left",
        "Right"
    ),
    y = c(
        2.22,
        1.90,
        0.92,
        0.60
    ),
    stringsAsFactors = FALSE
)

block_headers <- data.frame(
    block = c(
        "Band-compatible",
        "Dennis-Ray-compatible"
    ),
    x = science_x_limits[1] + 0.0015,
    y = c(
        2.42,
        1.12
    ),
    label = c(
        "Band-compatible",
        "Dennis / Ray-compatible"
    ),
    stringsAsFactors = FALSE
)

side_labels <- data.frame(
    block = c(
        "Band-compatible",
        "Band-compatible",
        "Dennis-Ray-compatible",
        "Dennis-Ray-compatible"
    ),
    x = science_x_limits[1] + 0.0090,
    y = c(
        2.22,
        1.90,
        0.92,
        0.60
    ),
    label = c(
        "Left",
        "Right",
        "Left",
        "Right"
    ),
    stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Background connectors and points
# ------------------------------------------------------------

background_connectors <- merge(
    background_sites[
        ,
        c(
            "block",
            "position_Mb",
            "pos",
            "highlight_B",
            "highlight_D",
            "highlight_R"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Left",
        c("block", "y"),
        drop = FALSE
    ],
    by = "block",
    all.x = TRUE,
    sort = FALSE
)

names(background_connectors)[
    names(background_connectors) == "y"
] <- "y_left"

tmp_right <- row_y[
    row_y$side == "Right",
    c("block", "y"),
    drop = FALSE
]

names(tmp_right)[
    names(tmp_right) == "y"
] <- "y_right"

background_connectors <- merge(
    background_connectors,
    tmp_right,
    by = "block",
    all.x = TRUE,
    sort = FALSE
)

background_left <- merge(
    background_sites[
        ,
        c(
            "block",
            "position_Mb",
            "pos",
            "left_state",
            "highlight_B",
            "highlight_D",
            "highlight_R"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Left",
        ,
        drop = FALSE
    ],
    by = "block",
    all.x = TRUE,
    sort = FALSE
)

background_left$state <- factor(
    background_left$left_state,
    levels = c("B/B", "B/DR", "DR/DR")
)

background_right <- merge(
    background_sites[
        ,
        c(
            "block",
            "position_Mb",
            "pos",
            "right_state",
            "highlight_B",
            "highlight_D",
            "highlight_R"
        ),
        drop = FALSE
    ],
    row_y[
        row_y$side == "Right",
        ,
        drop = FALSE
    ],
    by = "block",
    all.x = TRUE,
    sort = FALSE
)

background_right$state <- factor(
    background_right$right_state,
    levels = c("B/B", "B/DR", "DR/DR")
)

background_points <- rbind(
    data.frame(
        block = background_left$block,
        side = background_left$side,
        y = background_left$y,
        position_Mb = background_left$position_Mb,
        pos = background_left$pos,
        state = background_left$state,
        highlight_B = background_left$highlight_B,
        highlight_D = background_left$highlight_D,
        highlight_R = background_left$highlight_R,
        stringsAsFactors = FALSE
    ),
    data.frame(
        block = background_right$block,
        side = background_right$side,
        y = background_right$y,
        position_Mb = background_right$position_Mb,
        pos = background_right$pos,
        state = background_right$state,
        highlight_B = background_right$highlight_B,
        highlight_D = background_right$highlight_D,
        highlight_R = background_right$highlight_R,
        stringsAsFactors = FALSE
    )
)

background_points$side <- factor(
    background_points$side,
    levels = c("Left", "Right")
)

background_points$state <- factor(
    background_points$state,
    levels = c("B/B", "B/DR", "DR/DR")
)


# ------------------------------------------------------------
# Highlight layers: B, D, and R as separate overlays
# ------------------------------------------------------------

make_highlight_layer <- function(
    x,
    block_name,
    trait_name,
    left_state,
    right_state
) {
    if (nrow(x) == 0L) {
        return(
            list(
                connectors = data.frame(),
                points = data.frame()
            )
        )
    }

    pos <- as.integer(
        round(
            x$position_Mb * 1e6
        )
    )

    connectors <- data.frame(
        block = block_name,
        position_Mb = x$position_Mb,
        pos = pos,
        trait = trait_name,
        stringsAsFactors = FALSE
    )

    connectors <- merge(
        connectors,
        row_y[
            row_y$side == "Left",
            c("block", "y"),
            drop = FALSE
        ],
        by = "block",
        all.x = TRUE,
        sort = FALSE
    )

    names(connectors)[
        names(connectors) == "y"
    ] <- "y_left"

    tmp <- row_y[
        row_y$side == "Right",
        c("block", "y"),
        drop = FALSE
    ]

    names(tmp)[
        names(tmp) == "y"
    ] <- "y_right"

    connectors <- merge(
        connectors,
        tmp,
        by = "block",
        all.x = TRUE,
        sort = FALSE
    )

    points <- rbind(
        data.frame(
            block = block_name,
            side = "Left",
            y = row_y$y[
                row_y$block == block_name &
                row_y$side == "Left"
            ],
            position_Mb = x$position_Mb,
            pos = pos,
            state = left_state,
            trait = trait_name,
            stringsAsFactors = FALSE
        ),
        data.frame(
            block = block_name,
            side = "Right",
            y = row_y$y[
                row_y$block == block_name &
                row_y$side == "Right"
            ],
            position_Mb = x$position_Mb,
            pos = pos,
            state = right_state,
            trait = trait_name,
            stringsAsFactors = FALSE
        )
    )

    points$state <- factor(
        points$state,
        levels = c("B/B", "B/DR", "DR/DR")
    )

    points$side <- factor(
        points$side,
        levels = c("Left", "Right")
    )

    list(
        connectors = connectors,
        points = points
    )
}

B_layers <- make_highlight_layer(
    highlight_B,
    block_name = "Band-compatible",
    trait_name = "B",
    left_state = "B/B",
    right_state = "B/DR"
)

D_layers <- make_highlight_layer(
    highlight_D,
    block_name = "Dennis-Ray-compatible",
    trait_name = "D",
    left_state = "DR/DR",
    right_state = "B/DR"
)

R_layers <- make_highlight_layer(
    highlight_R,
    block_name = "Dennis-Ray-compatible",
    trait_name = "R",
    left_state = "DR/DR",
    right_state = "B/DR"
)

highlight_connectors <- rbind(
    B_layers$connectors,
    D_layers$connectors,
    R_layers$connectors
)

highlight_points <- rbind(
    B_layers$points,
    D_layers$points,
    R_layers$points
)

highlight_connectors$trait <- factor(
    highlight_connectors$trait,
    levels = c("B", "D", "R")
)

highlight_points$trait <- factor(
    highlight_points$trait,
    levels = c("B", "D", "R")
)


# ------------------------------------------------------------
# Published regions
# ------------------------------------------------------------

published_regions <- data.frame(
    region = c("B", "R1", "R2", "D1", "D2"),
    xmin = c(1112651, 1156746, 1174547, 1215210, 1218557) / 1e6,
    xmax = c(1114314, 1160220, 1187833, 1216124, 1221494) / 1e6,
    trait = c("B", "R", "R", "D", "D"),
    stringsAsFactors = FALSE
)

published_regions$y <- 5.98
published_regions$label_y <- 6.08


# ============================================================
# AUDIT / SUMMARY
# ============================================================

audit <- background_sites[
    ,
    c(
        "block",
        "pos",
        "position_Mb",
        "left_state",
        "right_state",
        "highlight_B",
        "highlight_D",
        "highlight_R"
    ),
    drop = FALSE
]

audit <- audit[
    order(
        match(
            audit$block,
            c(
                "Band-compatible",
                "Dennis-Ray-compatible"
            )
        ),
        audit$position_Mb
    ),
    ,
    drop = FALSE
]

write.table(
    audit,
    file = audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

DR_overlap_n <- sum(
    positions_to_key(
        highlight_D$position_Mb,
        "Dennis-Ray-compatible"
    ) %in%
    positions_to_key(
        highlight_R$position_Mb,
        "Dennis-Ray-compatible"
    )
)

summary_table <- data.frame(
    block = c(
        "Band-compatible",
        "Dennis-Ray-compatible"
    ),
    exact_background_sites = c(
        nrow(band_bg),
        nrow(dr_bg)
    ),
    B_high_BA_sites = c(
        nrow(highlight_B),
        0
    ),
    D_high_BA_sites = c(
        0,
        nrow(highlight_D)
    ),
    R_high_BA_sites = c(
        0,
        nrow(highlight_R)
    ),
    left_state = c(
        "B/B",
        "DR/DR"
    ),
    right_state = c(
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

    # --------------------------------------------------------
    # Published intervals
    # --------------------------------------------------------

    geom_segment(
        data = published_regions,
        aes(
            x = xmin,
            xend = xmax,
            y = y,
            yend = y,
            colour = trait
        ),
        linewidth = 2.40,
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

    # --------------------------------------------------------
    # optix
    # --------------------------------------------------------

    geom_segment(
        aes(
            x = optix_start_science,
            xend = optix_end_science,
            y = 5.58,
            yend = 5.58
        ),
        inherit.aes = FALSE,
        colour = optix_colour,
        linewidth = 2.95,
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

    # --------------------------------------------------------
    # Upper rug baselines
    # --------------------------------------------------------

    geom_segment(
        data = data.frame(
            trait = c("B", "D", "R"),
            y = c(5.15, 4.35, 3.55),
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

    # BA > 0.99 sites highlighted
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
        linewidth = 0.76,
        alpha = 1,
        lineend = "butt",
        show.legend = FALSE
    ) +

    # trait labels
    geom_text(
        data = data.frame(
            trait = c("B", "D", "R"),
            x = science_x_limits[1] + 0.0015,
            y = c(5.33, 4.53, 3.73),
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
        size = 3.05,
        show.legend = FALSE
    ) +

    # section separator
    geom_segment(
        aes(
            x = science_x_limits[1],
            xend = science_x_limits[2],
            y = 2.78,
            yend = 2.78
        ),
        inherit.aes = FALSE,
        colour = "#D4D7DA",
        linewidth = 0.36
    ) +

    annotate(
        "text",
        x = science_x_limits[1] + 0.0015,
        y = 2.90,
        label = "Focal L vs R exact-compatible sites",
        hjust = 0,
        vjust = 0,
        size = 2.38,
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
        linewidth = 0.26
    ) +

    # block labels
    geom_text(
        data = block_headers,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 0,
        fontface = "bold",
        size = 2.90,
        colour = "#444444"
    ) +

    # side labels
    geom_text(
        data = side_labels,
        aes(
            x = x,
            y = y,
            label = label
        ),
        hjust = 0,
        size = 2.12,
        colour = "#5D646B"
    ) +

    # background connectors
    geom_segment(
        data = background_connectors,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left
        ),
        colour = "#C3C9CF",
        linewidth = 0.26,
        alpha = 0.68
    ) +

    # background points
    geom_point(
        data = background_points,
        aes(
            x = position_Mb,
            y = y,
            fill = state
        ),
        shape = 21,
        size = 1.52,
        stroke = 0.20,
        colour = "#6C6C6C",
        alpha = 0.92
    ) +

    # highlighted connectors
    geom_segment(
        data = highlight_connectors,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y_right,
            yend = y_left,
            colour = trait
        ),
        linewidth = 0.70,
        alpha = 0.96,
        show.legend = FALSE
    ) +

    # highlighted points
    geom_point(
        data = highlight_points,
        aes(
            x = position_Mb,
            y = y,
            fill = state,
            colour = trait
        ),
        shape = 21,
        size = 2.52,
        stroke = 0.78,
        alpha = 1,
        show.legend = FALSE
    ) +

    scale_fill_manual(
        values = state_colours,
        breaks = c("B/B", "B/DR", "DR/DR"),
        labels = c(
            "B/B (B only)",
            "B/DR (both)",
            "DR/DR (DR only)"
        )
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
        ylim = c(0.35, 6.28),
        clip = "off"
    ) +

    labs(
        title = "B, Dennis and Ray association with focal Left-vs-Right genotype matches",
        subtitle = paste(
            "Upper panel: compact BA rugs across the region.",
            "Lower panel: two explicit Left-vs-Right genotype blocks;",
            "all exact-compatible sites are shown, with BA > 0.99 trait-specific sites highlighted."
        ),
        caption = paste(
            "Band-compatible block: Left = B/B, Right = B/DR.",
            "Dennis-Ray-compatible block: Left = DR/DR, Right = B/DR.",
            "B highlights n =", nrow(highlight_B), ";",
            "D highlights n =", nrow(highlight_D), ";",
            "R highlights n =", nrow(highlight_R), ";",
            "D/R overlap =", DR_overlap_n, "."
        )
    ) +

    theme_classic(
        base_size = 8
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 9.7,
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
# LEFT SIDECAR: CORRECT GENOTYPE KEY
# ============================================================

make_state_sidecar <- function() {
    grobTree(

        textGrob(
            "Genotype key",
            x = unit(0.08, "npc"),
            y = unit(0.96, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8.0,
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
                height = 0.15
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

        # B/DR
        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.14,
                y = 0.50,
                width = 0.16,
                height = 0.11
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.24,
                y = 0.50,
                width = 0.16,
                height = 0.11
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
                height = 0.15
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
    width = 9.15,
    height = 5.45
)

draw_png(
    final_grob,
    png_file,
    width = 9.15,
    height = 5.45
)


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 165 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nDennis/Ray high-BA overlap sites:",
    DR_overlap_n,
    "\n",
    sep = ""
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
