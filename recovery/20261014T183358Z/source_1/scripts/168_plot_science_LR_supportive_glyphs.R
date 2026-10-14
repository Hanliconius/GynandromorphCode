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
    "science_LR_supportive_glyphs"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

summary_file <- file.path(
    output_dir,
    "Science_LR_supportive_glyphs.summary.tsv"
)

audit_file <- file.path(
    output_dir,
    "Science_LR_supportive_glyphs.audit.tsv"
)

comparison_pdf <- file.path(
    output_dir,
    "Science_LR_supportive_glyphs.comparison.pdf"
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

extract_position_table <- function(x) {
    if (nrow(x) == 0L) {
        return(
            data.frame(
                pos = integer(),
                position_Mb = numeric(),
                stringsAsFactors = FALSE
            )
        )
    }

    pos_col <- optional_column(
        x,
        c("pos")
    )

    posMb_col <- resolve_column(
        x,
        c(
            "position_Mb",
            "position_mb",
            "pos_Mb"
        )
    )

    position_Mb <- suppressWarnings(
        as.numeric(
            x[[posMb_col]]
        )
    )

    keep <- is.finite(position_Mb) &
        position_Mb >= science_x_limits[1] &
        position_Mb <= science_x_limits[2]

    position_Mb <- position_Mb[keep]

    if (is.na(pos_col)) {
        pos <- as.integer(
            round(position_Mb * 1e6)
        )
    } else {
        pos <- suppressWarnings(
            as.integer(
                x[[pos_col]][keep]
            )
        )

        missing_pos <- !is.finite(pos)

        pos[missing_pos] <- as.integer(
            round(
                position_Mb[missing_pos] * 1e6
            )
        )
    }

    d <- data.frame(
        pos = pos,
        position_Mb = position_Mb,
        stringsAsFactors = FALSE
    )

    d <- d[
        order(d$position_Mb),
        ,
        drop = FALSE
    ]

    d[
        !duplicated(d$pos),
        ,
        drop = FALSE
    ]
}

extract_support_table <- function(
    x,
    support_group
) {
    if (nrow(x) == 0L) {
        return(
            data.frame(
                support_group = character(),
                pos = integer(),
                position_Mb = numeric(),
                BA = numeric(),
                BA_is_one = logical(),
                stringsAsFactors = FALSE
            )
        )
    }

    pos_col <- optional_column(
        x,
        c("pos")
    )

    posMb_col <- resolve_column(
        x,
        c(
            "position_Mb",
            "position_mb",
            "pos_Mb"
        )
    )

    ba_col <- optional_column(
        x,
        c(
            "BA",
            "balanced_accuracy",
            "Band_BA",
            "Dennis_BA",
            "Ray_BA",
            "B_BA",
            "D_BA",
            "R_BA",
            "band_ba",
            "dennis_ba",
            "ray_ba"
        )
    )

    position_Mb <- suppressWarnings(
        as.numeric(
            x[[posMb_col]]
        )
    )

    keep <- is.finite(position_Mb) &
        position_Mb >= science_x_limits[1] &
        position_Mb <= science_x_limits[2]

    position_Mb <- position_Mb[keep]

    if (is.na(pos_col)) {
        pos <- as.integer(
            round(position_Mb * 1e6)
        )
    } else {
        pos <- suppressWarnings(
            as.integer(
                x[[pos_col]][keep]
            )
        )

        missing_pos <- !is.finite(pos)

        pos[missing_pos] <- as.integer(
            round(
                position_Mb[missing_pos] * 1e6
            )
        )
    }

    if (is.na(ba_col)) {
        BA <- rep(
            NA_real_,
            length(pos)
        )
    } else {
        BA <- suppressWarnings(
            as.numeric(
                x[[ba_col]][keep]
            )
        )
    }

    d <- data.frame(
        support_group = support_group,
        pos = pos,
        position_Mb = position_Mb,
        BA = BA,
        BA_is_one = is.finite(BA) & abs(BA - 1) < 1e-12,
        stringsAsFactors = FALSE
    )

    d <- d[
        order(d$position_Mb),
        ,
        drop = FALSE
    ]

    d[
        !duplicated(d$pos),
        ,
        drop = FALSE
    ]
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

blue_B <- "#3A7CA5"
magenta_DR <- "#C46A97"
grey_bg <- "#C9CED4"
grey_dark <- "#69727C"
text_dark <- "#242424"
separator_grey <- "#D4D7DA"
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
# LOWER PANEL DATA
# ============================================================

band_background <- extract_position_table(
    band_correct
)

dr_background <- extract_position_table(
    dr_correct
)

background_positions <- rbind(
    band_background,
    dr_background
)

background_positions <- background_positions[
    order(background_positions$position_Mb),
    ,
    drop = FALSE
]

background_positions <- background_positions[
    !duplicated(background_positions$pos),
    ,
    drop = FALSE
]

support_B <- extract_support_table(
    high_B,
    "B"
)

support_D <- extract_support_table(
    high_D,
    "D"
)

support_R <- extract_support_table(
    high_R,
    "R"
)

support_DR_union <- merge(
    support_D[
        ,
        c("pos", "position_Mb", "BA", "BA_is_one"),
        drop = FALSE
    ],
    support_R[
        ,
        c("pos", "position_Mb", "BA", "BA_is_one"),
        drop = FALSE
    ],
    by = "pos",
    all = TRUE,
    suffixes = c("_D", "_R"),
    sort = TRUE
)

support_DR_union$position_Mb <- ifelse(
    is.finite(support_DR_union$position_Mb_D),
    support_DR_union$position_Mb_D,
    support_DR_union$position_Mb_R
)

support_DR_union$position_Mb <- suppressWarnings(
    as.numeric(
        support_DR_union$position_Mb
    )
)

support_DR_union$subtype <- ifelse(
    is.finite(support_DR_union$BA_D) &
    is.finite(support_DR_union$BA_R),
    "D+R",
    ifelse(
        is.finite(support_DR_union$BA_D),
        "D",
        "R"
    )
)

support_DR_union$BA <- pmax(
    support_DR_union$BA_D,
    support_DR_union$BA_R,
    na.rm = TRUE
)

support_DR_union$BA_is_one <- (
    ifelse(
        is.na(support_DR_union$BA_is_one_D),
        FALSE,
        support_DR_union$BA_is_one_D
    ) |
    ifelse(
        is.na(support_DR_union$BA_is_one_R),
        FALSE,
        support_DR_union$BA_is_one_R
    )
)

support_DR_union <- support_DR_union[
    ,
    c(
        "pos",
        "position_Mb",
        "subtype",
        "BA",
        "BA_is_one"
    ),
    drop = FALSE
]

support_DR_union <- support_DR_union[
    order(support_DR_union$position_Mb),
    ,
    drop = FALSE
]

left_y <- 1.82
right_y <- 1.28

bg_half_height <- 0.06
support_half_height <- 0.10
both_offset <- 0.00033
ba1_cap_y_left <- left_y + 0.16
ba1_cap_y_right <- right_y + 0.16
ba1_cap_half_width <- 0.00025

background_ticks <- rbind(
    data.frame(
        row = "Left",
        position_Mb = background_positions$position_Mb,
        y0 = left_y - bg_half_height,
        y1 = left_y + bg_half_height,
        stringsAsFactors = FALSE
    ),
    data.frame(
        row = "Right",
        position_Mb = background_positions$position_Mb,
        y0 = right_y - bg_half_height,
        y1 = right_y + bg_half_height,
        stringsAsFactors = FALSE
    )
)

# 168A / 168B core supportive glyphs
B_left_ticks <- data.frame(
    position_Mb = support_B$position_Mb,
    y0 = left_y - support_half_height,
    y1 = left_y + support_half_height,
    BA_is_one = support_B$BA_is_one,
    stringsAsFactors = FALSE
)

B_right_ticks_blue <- data.frame(
    position_Mb = support_B$position_Mb - both_offset,
    y0 = right_y - support_half_height,
    y1 = right_y + support_half_height,
    BA_is_one = support_B$BA_is_one,
    stringsAsFactors = FALSE
)

B_right_ticks_magenta <- data.frame(
    position_Mb = support_B$position_Mb + both_offset,
    y0 = right_y - support_half_height,
    y1 = right_y + support_half_height,
    BA_is_one = support_B$BA_is_one,
    stringsAsFactors = FALSE
)

DR_left_ticks <- data.frame(
    position_Mb = support_DR_union$position_Mb,
    y0 = left_y - support_half_height,
    y1 = left_y + support_half_height,
    BA_is_one = support_DR_union$BA_is_one,
    subtype = support_DR_union$subtype,
    stringsAsFactors = FALSE
)

DR_right_ticks_blue <- data.frame(
    position_Mb = support_DR_union$position_Mb - both_offset,
    y0 = right_y - support_half_height,
    y1 = right_y + support_half_height,
    BA_is_one = support_DR_union$BA_is_one,
    subtype = support_DR_union$subtype,
    stringsAsFactors = FALSE
)

DR_right_ticks_magenta <- data.frame(
    position_Mb = support_DR_union$position_Mb + both_offset,
    y0 = right_y - support_half_height,
    y1 = right_y + support_half_height,
    BA_is_one = support_DR_union$BA_is_one,
    subtype = support_DR_union$subtype,
    stringsAsFactors = FALSE
)

B_ba1_caps_left <- data.frame(
    x0 = support_B$position_Mb - ba1_cap_half_width,
    x1 = support_B$position_Mb + ba1_cap_half_width,
    y = ba1_cap_y_left,
    stringsAsFactors = FALSE
)

B_ba1_caps_left <- B_ba1_caps_left[
    support_B$BA_is_one,
    ,
    drop = FALSE
]

B_ba1_caps_right <- data.frame(
    x0 = support_B$position_Mb - ba1_cap_half_width,
    x1 = support_B$position_Mb + ba1_cap_half_width,
    y = ba1_cap_y_right,
    stringsAsFactors = FALSE
)

B_ba1_caps_right <- B_ba1_caps_right[
    support_B$BA_is_one,
    ,
    drop = FALSE
]

DR_ba1_caps_left <- data.frame(
    x0 = support_DR_union$position_Mb - ba1_cap_half_width,
    x1 = support_DR_union$position_Mb + ba1_cap_half_width,
    y = ba1_cap_y_left,
    stringsAsFactors = FALSE
)

DR_ba1_caps_left <- DR_ba1_caps_left[
    support_DR_union$BA_is_one,
    ,
    drop = FALSE
]

DR_ba1_caps_right <- data.frame(
    x0 = support_DR_union$position_Mb - ba1_cap_half_width,
    x1 = support_DR_union$position_Mb + ba1_cap_half_width,
    y = ba1_cap_y_right,
    stringsAsFactors = FALSE
)

DR_ba1_caps_right <- DR_ba1_caps_right[
    support_DR_union$BA_is_one,
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

published_regions$y <- 5.90
published_regions$label_y <- 6.01


# ============================================================
# SUMMARY / AUDIT
# ============================================================

summary_table <- data.frame(
    statistic = c(
        "background_positions_unique",
        "B_supportive_sites",
        "D_supportive_sites",
        "R_supportive_sites",
        "D_R_supportive_union",
        "D_R_supportive_overlap",
        "B_BA_equals_1",
        "D_BA_equals_1",
        "R_BA_equals_1",
        "D_R_union_BA_equals_1"
    ),
    value = c(
        nrow(background_positions),
        nrow(support_B),
        nrow(support_D),
        nrow(support_R),
        nrow(support_DR_union),
        sum(support_D$pos %in% support_R$pos),
        sum(support_B$BA_is_one),
        sum(support_D$BA_is_one),
        sum(support_R$BA_is_one),
        sum(support_DR_union$BA_is_one)
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

audit <- rbind(
    data.frame(
        class = "background",
        support_group = "all",
        pos = background_positions$pos,
        position_Mb = background_positions$position_Mb,
        subtype = NA_character_,
        BA = NA_real_,
        BA_is_one = NA,
        left_symbol = "grey tick",
        right_symbol = "grey tick",
        stringsAsFactors = FALSE
    ),
    data.frame(
        class = "supportive",
        support_group = "B",
        pos = support_B$pos,
        position_Mb = support_B$position_Mb,
        subtype = "B",
        BA = support_B$BA,
        BA_is_one = support_B$BA_is_one,
        left_symbol = "single blue tick",
        right_symbol = "paired blue+magenta ticks",
        stringsAsFactors = FALSE
    ),
    data.frame(
        class = "supportive",
        support_group = "DR_union",
        pos = support_DR_union$pos,
        position_Mb = support_DR_union$position_Mb,
        subtype = support_DR_union$subtype,
        BA = support_DR_union$BA,
        BA_is_one = support_DR_union$BA_is_one,
        left_symbol = "single magenta tick",
        right_symbol = "paired blue+magenta ticks",
        stringsAsFactors = FALSE
    )
)

audit <- audit[
    order(
        audit$position_Mb,
        audit$class
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


# ============================================================
# BASE PLOT
# ============================================================

build_base_plot <- function(
    title_text,
    subtitle_text,
    caption_text
) {
    ggplot() +

        # Published intervals
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

        # Upper B/D/R rug baselines
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
            colour = "#AEB6BE",
            linewidth = 0.28
        ) +

        # All BA sites
        geom_segment(
            data = ba,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = rug_base,
                yend = rug_top
            ),
            inherit.aes = FALSE,
            colour = grey_dark,
            alpha = 0.32,
            linewidth = 0.23
        ) +

        # BA > 0.99
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

        # section separator
        geom_segment(
            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = 2.62,
                yend = 2.62
            ),
            inherit.aes = FALSE,
            colour = separator_grey,
            linewidth = 0.35
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = 2.72,
            label = "Focal Left vs Right supportive-SNP track",
            hjust = 0,
            vjust = 0,
            size = 2.36,
            fontface = "bold",
            colour = "#505050"
        ) +

        annotate(
            "text",
            x = science_x_limits[2] - 0.0015,
            y = 2.72,
            label = "grey = all focal-concordant sites; coloured glyphs = supportive SNPs",
            hjust = 1,
            vjust = 0,
            size = 2.02,
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
            colour = "#E0E4E8",
            linewidth = 0.28
        ) +

        geom_segment(
            aes(
                x = science_x_limits[1],
                xend = science_x_limits[2],
                y = right_y,
                yend = right_y
            ),
            inherit.aes = FALSE,
            colour = "#E0E4E8",
            linewidth = 0.28
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = left_y + 0.14,
            label = "L",
            hjust = 0,
            fontface = "bold",
            size = 3.2,
            colour = "#444444"
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
            y = right_y + 0.14,
            label = "R",
            hjust = 0,
            fontface = "bold",
            size = 3.2,
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

        scale_colour_manual(
            values = trait_colours,
            drop = FALSE
        ) +

        coord_cartesian(
            xlim = science_x_limits,
            ylim = c(1.00, 6.20),
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
                colour = text_dark,
                margin = margin(b = 2)
            ),
            plot.subtitle = element_text(
                size = 7.0,
                colour = "#505050",
                margin = margin(b = 4)
            ),
            plot.caption = element_text(
                size = 6.2,
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
}


# ============================================================
# SIDECAR LEGENDS
# ============================================================

make_sidecar <- function(include_ba1 = FALSE) {
    gt <- grobTree(

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
            "Lower-panel key",
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
            "Background",
            x = unit(0.06, "npc"),
            y = unit(0.90, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(c(0.10, 0.10), "npc"),
            x1 = unit(c(0.10, 0.10), "npc"),
            y0 = unit(c(0.82, 0.74), "npc"),
            y1 = unit(c(0.86, 0.78), "npc"),
            gp = gpar(
                col = grey_bg,
                lwd = 2.0
            )
        ),

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.84, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.76, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "grey ticks = all focal-concordant sites",
            x = unit(0.18, "npc"),
            y = unit(0.80, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "B-supportive SNP",
            x = unit(0.06, "npc"),
            y = unit(0.66, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.60, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.52, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        segmentsGrob(
            x0 = unit(0.10, "npc"),
            x1 = unit(0.10, "npc"),
            y0 = unit(0.57, "npc"),
            y1 = unit(0.63, "npc"),
            gp = gpar(
                col = blue_B,
                lwd = 3.2
            )
        ),

        segmentsGrob(
            x0 = unit(c(0.095, 0.105), "npc"),
            x1 = unit(c(0.095, 0.105), "npc"),
            y0 = unit(c(0.49, 0.49), "npc"),
            y1 = unit(c(0.55, 0.55), "npc"),
            gp = gpar(
                col = c(blue_B, magenta_DR),
                lwd = 3.0
            )
        ),

        textGrob(
            "Left = B/B",
            x = unit(0.18, "npc"),
            y = unit(0.60, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = blue_B
            )
        ),

        textGrob(
            "Right = B/DR",
            x = unit(0.18, "npc"),
            y = unit(0.52, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "D/R-supportive SNP",
            x = unit(0.06, "npc"),
            y = unit(0.38, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.32, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.24, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.3,
                fontface = "bold",
                col = "#4F4F4F"
            )
        ),

        segmentsGrob(
            x0 = unit(0.10, "npc"),
            x1 = unit(0.10, "npc"),
            y0 = unit(0.29, "npc"),
            y1 = unit(0.35, "npc"),
            gp = gpar(
                col = magenta_DR,
                lwd = 3.2
            )
        ),

        segmentsGrob(
            x0 = unit(c(0.095, 0.105), "npc"),
            x1 = unit(c(0.095, 0.105), "npc"),
            y0 = unit(c(0.21, 0.21), "npc"),
            y1 = unit(c(0.27, 0.27), "npc"),
            gp = gpar(
                col = c(blue_B, magenta_DR),
                lwd = 3.0
            )
        ),

        textGrob(
            "Left = DR/DR",
            x = unit(0.18, "npc"),
            y = unit(0.32, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = magenta_DR
            )
        ),

        textGrob(
            "Right = B/DR",
            x = unit(0.18, "npc"),
            y = unit(0.24, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#4F4F4F"
            )
        ),

        textGrob(
            "Cartoon guide",
            x = unit(0.06, "npc"),
            y = unit(0.14, "npc"),
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
                x = 0.10,
                y = 0.09,
                width = 0.10,
                height = 0.07
            )
        ),

        textGrob(
            "B/B",
            x = unit(0.17, "npc"),
            y = unit(0.09, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = blue_B
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.09,
                y = 0.05,
                width = 0.07,
                height = 0.05
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.13,
                y = 0.05,
                width = 0.07,
                height = 0.05
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.17, "npc"),
            y = unit(0.05, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#4F4F4F"
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.10,
                y = 0.015,
                width = 0.10,
                height = 0.07
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.17, "npc"),
            y = unit(0.015, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = magenta_DR
            )
        )
    )

    if (include_ba1) {
        gt <- grobTree(
            gt,
            textGrob(
                "Black cap = BA = 1.0",
                x = unit(0.06, "npc"),
                y = unit(0.705, "npc"),
                hjust = 0,
                gp = gpar(
                    fontsize = 5.8,
                    col = "#222222"
                )
            ),
            segmentsGrob(
                x0 = unit(0.095, "npc"),
                x1 = unit(0.105, "npc"),
                y0 = unit(0.705, "npc"),
                y1 = unit(0.705, "npc"),
                gp = gpar(
                    col = "#111111",
                    lwd = 1.8
                )
            )
        )
    }

    gt
}

add_sidecar <- function(
    plot,
    include_ba1 = FALSE,
    sidecar_width = unit(2.30, "in")
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
        grobs = make_sidecar(
            include_ba1 = include_ba1
        ),
        t = panel_t,
        b = panel_b,
        l = panel_l,
        r = panel_l,
        clip = "off",
        name = "lower-panel-sidecar"
    )

    g
}


# ============================================================
# 168A: CLEAN SUPPORTIVE GLYPHS
# ============================================================

plot_168A <- build_base_plot(
    title_text = "Compact B/D/R association panel with focal Left–Right supportive-SNP track",
    subtitle_text = paste(
        "168A: lower panel shows all focal-concordant sites as grey ticks;",
        "supportive SNPs are drawn as interpretable Left/Right glyphs."
    ),
    caption_text = paste(
        "Lower panel glyph grammar:",
        "B-supportive = Left single blue tick (B/B), Right paired blue+magenta ticks (B/DR).",
        "D/R-supportive = Left single magenta tick (DR/DR), Right paired blue+magenta ticks (B/DR)."
    )
) +

    geom_segment(
        data = background_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = grey_bg,
        linewidth = 0.26,
        alpha = 0.70
    ) +

    geom_segment(
        data = B_left_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.92,
        alpha = 1
    ) +

    geom_segment(
        data = B_right_ticks_blue,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = B_right_ticks_magenta,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = DR_left_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.92,
        alpha = 1
    ) +

    geom_segment(
        data = DR_right_ticks_blue,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = DR_right_ticks_magenta,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.90,
        alpha = 1
    )


# ============================================================
# 168B: SAME, PLUS BA = 1.0 CAPS
# ============================================================

plot_168B <- build_base_plot(
    title_text = "Compact B/D/R association panel with focal Left–Right supportive-SNP track",
    subtitle_text = paste(
        "168B: as in 168A, with black caps marking supportive SNPs with BA = 1.0."
    ),
    caption_text = paste(
        "Lower panel glyph grammar as in 168A.",
        "Black caps mark supportive SNPs with BA = 1.0."
    )
) +

    geom_segment(
        data = background_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = grey_bg,
        linewidth = 0.26,
        alpha = 0.70
    ) +

    geom_segment(
        data = B_left_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.92,
        alpha = 1
    ) +

    geom_segment(
        data = B_right_ticks_blue,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = B_right_ticks_magenta,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = DR_left_ticks,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.92,
        alpha = 1
    ) +

    geom_segment(
        data = DR_right_ticks_blue,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = blue_B,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = DR_right_ticks_magenta,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = magenta_DR,
        linewidth = 0.90,
        alpha = 1
    ) +

    geom_segment(
        data = B_ba1_caps_left,
        aes(
            x = x0,
            xend = x1,
            y = y,
            yend = y
        ),
        inherit.aes = FALSE,
        colour = "#111111",
        linewidth = 0.62
    ) +

    geom_segment(
        data = B_ba1_caps_right,
        aes(
            x = x0,
            xend = x1,
            y = y,
            yend = y
        ),
        inherit.aes = FALSE,
        colour = "#111111",
        linewidth = 0.62
    ) +

    geom_segment(
        data = DR_ba1_caps_left,
        aes(
            x = x0,
            xend = x1,
            y = y,
            yend = y
        ),
        inherit.aes = FALSE,
        colour = "#111111",
        linewidth = 0.62
    ) +

    geom_segment(
        data = DR_ba1_caps_right,
        aes(
            x = x0,
            xend = x1,
            y = y,
            yend = y
        ),
        inherit.aes = FALSE,
        colour = "#111111",
        linewidth = 0.62
    )


# ============================================================
# SAVE
# ============================================================

plot_list <- list(
    variant_168A_clean = add_sidecar(
        plot_168A,
        include_ba1 = FALSE
    ),
    variant_168B_BA1caps = add_sidecar(
        plot_168B,
        include_ba1 = TRUE
    )
)

for (nm in names(plot_list)) {
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
        plot_list[[nm]],
        pdf_file,
        width = 9.15,
        height = 5.25
    )

    draw_png(
        plot_list[[nm]],
        png_file,
        width = 9.15,
        height = 5.25
    )
}

cairo_pdf(
    filename = comparison_pdf,
    width = 9.15,
    height = 5.25,
    onefile = TRUE
)

for (nm in names(plot_list)) {
    grid.newpage()
    grid.draw(
        plot_list[[nm]]
    )
}

dev.off()


# ============================================================
# TERMINAL REPORT
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 168 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nComparison PDF:\n",
    comparison_pdf,
    "\n\nSummary:\n",
    summary_file,
    "\n\nAudit:\n",
    audit_file,
    "\n",
    sep = ""
)

cat(
    "\nIndividual files:\n",
    paste(
        file.path(
            output_dir,
            c(
                "Science_variant_168A_clean.pdf",
                "Science_variant_168A_clean.png",
                "Science_variant_168B_BA1caps.pdf",
                "Science_variant_168B_BA1caps.png"
            )
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)
