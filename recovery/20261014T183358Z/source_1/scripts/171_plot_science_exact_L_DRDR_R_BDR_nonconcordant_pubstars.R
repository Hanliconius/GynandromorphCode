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

mixture_file <- file.path(
    project,
    "analysis",
    "original_502",
    "reciprocal_L_R_mixture",
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing required source script:",
            source_file
        )
    )
}

if (!file.exists(mixture_file)) {
    stop(
        paste(
            "Missing reciprocal-mixture focal genotype file:",
            mixture_file
        )
    )
}


# ============================================================
# LOAD ESTABLISHED ANALYSIS-159 OBJECTS
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
            paste(missing_objects, collapse = ", ")
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
    "science_exact_L_DRDR_R_BDR_nonconcordant_pubstars"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

pdf_file <- file.path(
    output_dir,
    "Science_171_exact_L_DRDR_R_BDR_nonconcordant_pubstars.pdf"
)

png_file <- file.path(
    output_dir,
    "Science_171_exact_L_DRDR_R_BDR_nonconcordant_pubstars.png"
)

summary_file <- file.path(
    output_dir,
    "Science_171_exact_L_DRDR_R_BDR_nonconcordant_pubstars.summary.tsv"
)

audit_file <- file.path(
    output_dir,
    "Science_171_exact_L_DRDR_R_BDR_nonconcordant_pubstars.audit.tsv"
)


# ============================================================
# HELPERS
# ============================================================

normalize_name <- function(x) {
    gsub(
        "[^a-z0-9]+",
        "",
        tolower(as.character(x)),
        perl = TRUE
    )
}

resolve_column <- function(x, candidates) {
    norm_cols <- normalize_name(names(x))

    for (candidate in candidates) {
        target <- normalize_name(candidate)
        hit <- which(norm_cols == target)

        if (length(hit) == 1L) {
            return(names(x)[hit])
        }
    }

    stop(
        paste0(
            "Could not resolve any of: ",
            paste(candidates, collapse = ", "),
            "\nAvailable columns:\n",
            paste(names(x), collapse = "\n")
        )
    )
}

optional_column <- function(x, candidates) {
    norm_cols <- normalize_name(names(x))

    for (candidate in candidates) {
        target <- normalize_name(candidate)
        hit <- which(norm_cols == target)

        if (length(hit) == 1L) {
            return(names(x)[hit])
        }
    }

    NA_character_
}

draw_pdf <- function(grob, filename, width, height) {
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

draw_png <- function(grob, filename, width, height, dpi = 450) {
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

normalize_gt <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- gsub("\\s+", "", x)
    x <- gsub("\\|", "/", x, perl = TRUE)

    out <- rep(NA_character_, length(x))

    out[x %in% c("P/P", "1/1", "DR/DR", "PP", "DRDR")] <- "DR/DR"
    out[x %in% c("P/A", "A/P", "0/1", "1/0", "B/DR", "DR/B", "PA", "AP", "BDR", "DRB", "HET")] <- "B/DR"
    out[x %in% c("A/A", "0/0", "B/B", "AA", "BB")] <- "B/B"

    out
}

extract_position_table <- function(x) {
    pos_col <- optional_column(
        x,
        c("pos", "POS", "position")
    )

    posMb_col <- optional_column(
        x,
        c("position_Mb", "position_mb", "pos_Mb")
    )

    if (is.na(pos_col) && is.na(posMb_col)) {
        stop("Could not resolve position columns.")
    }

    if (!is.na(pos_col)) {
        pos <- suppressWarnings(as.integer(x[[pos_col]]))
    } else {
        pos <- rep(NA_integer_, nrow(x))
    }

    if (!is.na(posMb_col)) {
        position_Mb <- suppressWarnings(as.numeric(x[[posMb_col]]))
    } else {
        position_Mb <- rep(NA_real_, nrow(x))
    }

    if (all(!is.finite(position_Mb)) && any(is.finite(pos))) {
        position_Mb <- pos / 1e6
    }

    if (all(!is.finite(pos)) && any(is.finite(position_Mb))) {
        pos <- as.integer(round(position_Mb * 1e6))
    }

    d <- data.frame(
        pos = pos,
        position_Mb = position_Mb,
        stringsAsFactors = FALSE
    )

    d <- d[
        is.finite(d$pos) &
        is.finite(d$position_Mb) &
        d$position_Mb >= science_x_limits[1] &
        d$position_Mb <= science_x_limits[2],
        ,
        drop = FALSE
    ]

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

collapse_by_position <- function(x) {
    x <- x[
        order(x$position_Mb),
        ,
        drop = FALSE
    ]

    split_x <- split(
        x,
        x$pos,
        drop = TRUE
    )

    out <- lapply(split_x, function(d) {
        d <- d[1, , drop = FALSE]

        data.frame(
            pos = d$pos,
            position_Mb = d$position_Mb,
            left_state = d$left_state,
            right_state = d$right_state,
            supportive = any(
                x$pos == d$pos &
                x$callable &
                x$supportive
            ),
            callable = any(
                x$pos == d$pos &
                x$callable
            ),
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, out)

    out <- out[
        order(out$position_Mb),
        ,
        drop = FALSE
    ]

    rownames(out) <- NULL
    out
}

region_hit_label <- function(position_Mb, published_regions) {
    hit <- published_regions$region[
        position_Mb >= published_regions$xmin &
        position_Mb <= published_regions$xmax
    ]

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    paste(hit, collapse = ";")
}


# ============================================================
# COLOURS
# ============================================================

# supplement genoplot scheme
HOM_AA_colour <- "#2C9BCB"
HET_colour    <- "#F28E2B"
HOM_PP_colour <- "#B61F72"

# focal states used in the lower panel
DR_colour <- HOM_PP_colour   # DR/DR
BDR_colour <- HET_colour     # B/DR

non_supportive_grey <- "#C8CDD3"
neutral_dark <- "#66707A"
neutral_mid <- "#AEB6BE"
neutral_light <- "#D9DEE3"
separator_grey <- "#D4D7DA"
text_colour <- "#242424"
optix_colour <- "#3A7CA5"


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

published_regions$y <- 5.86
published_regions$label_y <- 5.97


# ============================================================
# UPPER B/D/R BA TRACKS
# ============================================================

ba_trait_column <- resolve_column(
    all_ba_science,
    c("trait")
)

ba_position_column <- resolve_column(
    all_ba_science,
    c("position_Mb", "position_mb")
)

ba_value_column <- resolve_column(
    all_ba_science,
    c("BA", "balanced_accuracy")
)

ba <- data.frame(
    trait = as.character(all_ba_science[[ba_trait_column]]),
    position_Mb = suppressWarnings(
        as.numeric(all_ba_science[[ba_position_column]])
    ),
    BA = suppressWarnings(
        as.numeric(all_ba_science[[ba_value_column]])
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

ba$pos <- as.integer(round(ba$position_Mb * 1e6))

ba$trait <- factor(
    ba$trait,
    levels = c("B", "D", "R"),
    ordered = TRUE
)

rug_bases <- c(
    B = 5.04,
    D = 4.27,
    R = 3.50
)

ba$rug_base <- unname(rug_bases[as.character(ba$trait)])

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

ba$rug_top <- ba$rug_base + ba$tick_height

high_ba <- ba[
    ba$BA > 0.99,
    ,
    drop = FALSE
]


# ============================================================
# SUPPORTIVE EXACT POSITIONS
#
# We preserve the established exact supportive set from dr_correct:
#   Left = DR/DR
#   Right = B/DR
# ============================================================

exact_focal <- extract_position_table(dr_correct)

if (nrow(exact_focal) == 0L) {
    stop("No exact supportive focal positions were recovered from dr_correct.")
}

exact_focal$supportive <- TRUE


# ============================================================
# ALL FOCAL L/R POSITIONS FROM MIXTURE FILE
# ============================================================

mixture <- read.table(
    mixture_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

mix_pos_col <- resolve_column(
    mixture,
    c("pos", "POS", "position")
)

mix_posMb_col <- optional_column(
    mixture,
    c("position_Mb", "position_mb", "pos_Mb")
)

mix_left_col <- resolve_column(
    mixture,
    c(
        "mixture_left_GT",
        "left_state",
        "left_GT",
        "left_gt"
    )
)

mix_right_col <- resolve_column(
    mixture,
    c(
        "mixture_right_GT",
        "right_state",
        "right_GT",
        "right_gt"
    )
)

posterior95_col <- resolve_column(
    mixture,
    c(
        "posterior95"
    )
)

mix_pos <- suppressWarnings(as.integer(mixture[[mix_pos_col]]))

if (is.na(mix_posMb_col)) {
    mix_posMb <- mix_pos / 1e6
} else {
    mix_posMb <- suppressWarnings(as.numeric(mixture[[mix_posMb_col]]))
    missing_mb <- !is.finite(mix_posMb) & is.finite(mix_pos)
    mix_posMb[missing_mb] <- mix_pos[missing_mb] / 1e6
}

mix_left <- normalize_gt(mixture[[mix_left_col]])
mix_right <- normalize_gt(mixture[[mix_right_col]])

mix_posterior95 <- toupper(
    trimws(
        as.character(
            mixture[[posterior95_col]]
        )
    )
) %in% c(
    "TRUE",
    "T",
    "1",
    "YES"
)

mixture_lr <- data.frame(
    pos = mix_pos,
    position_Mb = mix_posMb,
    left_state = mix_left,
    right_state = mix_right,
    posterior95 = mix_posterior95,
    stringsAsFactors = FALSE
)

mixture_lr <- mixture_lr[
    is.finite(mixture_lr$pos) &
    is.finite(mixture_lr$position_Mb) &
    mixture_lr$position_Mb >= science_x_limits[1] &
    mixture_lr$position_Mb <= science_x_limits[2],
    ,
    drop = FALSE
]

mixture_lr$callable <- (
    mixture_lr$posterior95 &
    !is.na(mixture_lr$left_state) &
    !is.na(mixture_lr$right_state)
)

mixture_lr$supportive <- (
    mixture_lr$callable &
    mixture_lr$left_state == "DR/DR" &
    mixture_lr$right_state == "B/DR"
)

mixture_unique <- collapse_by_position(mixture_lr)

callable_positions <- mixture_unique[
    mixture_unique$callable,
    ,
    drop = FALSE
]

non_supportive_positions <- callable_positions[
    !callable_positions$supportive,
    c("pos", "position_Mb"),
    drop = FALSE
]

# Optional validation only; do not let slight parsing differences redefine the supportive set
mixture_supportive_positions <- callable_positions[
    callable_positions$supportive,
    c("pos", "position_Mb"),
    drop = FALSE
]


# ============================================================
# HIGH-BA SUPPORT AT EXACT SUPPORTIVE POSITIONS
# ============================================================

exact_key <- exact_focal$pos

high_exact <- high_ba[
    high_ba$pos %in% exact_key,
    ,
    drop = FALSE
]

high_exact <- high_exact[
    order(high_exact$position_Mb, high_exact$trait),
    ,
    drop = FALSE
]

support_for_position <- function(position) {
    traits <- as.character(high_exact$trait[high_exact$pos == position])
    traits <- unique(traits)

    if (length(traits) == 0L) {
        return("none")
    }

    paste(traits, collapse = "+")
}

exact_focal$high_BA_support <- vapply(
    exact_focal$pos,
    support_for_position,
    character(1)
)

exact_focal$high_BA_any <- exact_focal$high_BA_support != "none"

exact_focal$published_region <- vapply(
    exact_focal$position_Mb,
    region_hit_label,
    character(1),
    published_regions = published_regions
)

exact_focal$in_published_region <- !is.na(exact_focal$published_region)


# ============================================================
# LOWER-PANEL GEOMETRY
# ============================================================

left_y <- 1.78
right_y <- 1.24

grey_half_height <- 0.040
supportive_half_height <- 0.060
high_half_height <- 0.125
star_y <- left_y + 0.22

non_supportive_left <- data.frame(
    position_Mb = non_supportive_positions$position_Mb,
    y0 = left_y - grey_half_height,
    y1 = left_y + grey_half_height,
    stringsAsFactors = FALSE
)

non_supportive_right <- data.frame(
    position_Mb = non_supportive_positions$position_Mb,
    y0 = right_y - grey_half_height,
    y1 = right_y + grey_half_height,
    stringsAsFactors = FALSE
)

supportive_all_left <- data.frame(
    position_Mb = exact_focal$position_Mb,
    y0 = left_y - supportive_half_height,
    y1 = left_y + supportive_half_height,
    stringsAsFactors = FALSE
)

supportive_all_right <- data.frame(
    position_Mb = exact_focal$position_Mb,
    y0 = right_y - supportive_half_height,
    y1 = right_y + supportive_half_height,
    stringsAsFactors = FALSE
)

high_supportive <- exact_focal[
    exact_focal$high_BA_any,
    ,
    drop = FALSE
]

high_left <- data.frame(
    position_Mb = high_supportive$position_Mb,
    y0 = left_y - high_half_height,
    y1 = left_y + high_half_height,
    stringsAsFactors = FALSE
)

high_right <- data.frame(
    position_Mb = high_supportive$position_Mb,
    y0 = right_y - high_half_height,
    y1 = right_y + high_half_height,
    stringsAsFactors = FALSE
)

pubstar_df <- exact_focal[
    exact_focal$in_published_region,
    c("pos", "position_Mb", "published_region", "high_BA_support"),
    drop = FALSE
]

pubstar_df$y <- star_y


# ============================================================
# AUDIT + SUMMARY
# ============================================================

audit <- rbind(
    data.frame(
        pos = non_supportive_positions$pos,
        position_Mb = non_supportive_positions$position_Mb,
        focal_class = "non_supportive",
        left_state = NA_character_,
        right_state = NA_character_,
        high_BA_support = "none",
        high_BA_any = FALSE,
        in_published_region = FALSE,
        published_region = NA_character_,
        stringsAsFactors = FALSE
    ),
    data.frame(
        pos = exact_focal$pos,
        position_Mb = exact_focal$position_Mb,
        focal_class = "supportive",
        left_state = "DR/DR",
        right_state = "B/DR",
        high_BA_support = exact_focal$high_BA_support,
        high_BA_any = exact_focal$high_BA_any,
        in_published_region = exact_focal$in_published_region,
        published_region = exact_focal$published_region,
        stringsAsFactors = FALSE
    )
)

audit <- audit[
    order(audit$position_Mb, audit$focal_class),
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

summary_table <- data.frame(
    statistic = c(
        "callable_focal_positions_from_mixture",
        "non_supportive_callable_positions",
        "supportive_positions_from_dr_correct",
        "supportive_positions_from_mixture",
        "supportive_overlap_dr_correct_vs_mixture",
        "supportive_positions_with_any_BA_gt_0.99",
        "supportive_positions_in_published_regions"
    ),
    value = c(
        nrow(callable_positions),
        nrow(non_supportive_positions),
        nrow(exact_focal),
        nrow(mixture_supportive_positions),
        sum(exact_focal$pos %in% mixture_supportive_positions$pos),
        sum(exact_focal$high_BA_any),
        sum(exact_focal$in_published_region)
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
# BASE PLOT
# ============================================================

build_base_plot <- function() {
    ggplot() +

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
            x = mean(c(optix_start_science, optix_end_science)),
            y = 5.62,
            label = "optix",
            colour = optix_colour,
            size = 2.42,
            fontface = "italic",
            vjust = 0
        ) +

        # BA baselines
        geom_segment(
            data = data.frame(
                trait = c("B", "D", "R"),
                y = unname(rug_bases),
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

        # all BA
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

        # BA > 0.99
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

        # trait labels
        geom_text(
            data = data.frame(
                trait = c("B", "D", "R"),
                x = science_x_limits[1] + 0.0015,
                y = unname(rug_bases) + 0.25,
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
                y = 2.53,
                yend = 2.53
            ),
            inherit.aes = FALSE,
            colour = separator_grey,
            linewidth = 0.35
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0015,
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
            x = science_x_limits[2] - 0.0010,
            y = 2.64,
            label = "grey = non-supportive; coloured = L DR/DR -> R B/DR; * = within published interval",
            hjust = 1,
            vjust = 0,
            size = 2.00,
            colour = "#666666"
        ) +

        # L/R guide rows
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
            x = science_x_limits[1] + 0.0013,
            y = left_y + 0.17,
            label = "L",
            hjust = 0,
            fontface = "bold",
            size = 3.45,
            colour = "#333333"
        ) +

        annotate(
            "text",
            x = science_x_limits[1] + 0.0013,
            y = right_y + 0.17,
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
            labels = sprintf("%.2f", science_x_breaks),
            expand = expansion(mult = c(0, 0))
        ) +

        coord_cartesian(
            xlim = science_x_limits,
            ylim = c(0.95, 6.15),
            clip = "off"
        ) +

        labs(
            title = "B, Dennis and Ray association with focal Left–Right genotype support",
            subtitle = paste(
                "171: non-supportive focal SNPs are shown in pale grey;",
                "supportive exact L = DR/DR, R = B/DR SNPs are coloured;",
                "high-BA supportive SNPs are taller;",
                "asterisks mark supportive SNPs inside previously annotated intervals."
            ),
            caption = paste(
                "Lower panel colour grammar follows the supplement genoplots:",
                "DR/DR = magenta, B/DR = orange."
            )
        ) +

        theme_classic(base_size = 8) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 9.5,
                colour = text_colour,
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
            axis.ticks.x = element_line(linewidth = 0.25),
            axis.line.x = element_line(linewidth = 0.30),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            plot.margin = margin(8, 18, 7, 10)
        )
}


# ============================================================
# SIDECAR KEY
# ============================================================

make_supportive_key <- function() {
    grobTree(
        textGrob(
            "Lower-panel key",
            x = unit(0.06, "npc"),
            y = unit(0.95, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8.2,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "Non-supportive focal SNP",
            x = unit(0.06, "npc"),
            y = unit(0.84, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.8,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.77, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.5,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.71, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.5,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(c(0.12, 0.12), "npc"),
            x1 = unit(c(0.12, 0.12), "npc"),
            y0 = unit(c(0.745, 0.685), "npc"),
            y1 = unit(c(0.795, 0.735), "npc"),
            gp = gpar(
                col = non_supportive_grey,
                lwd = 2.2
            )
        ),

        textGrob(
            "grey = callable focal SNP that does not match the exact supportive pattern",
            x = unit(0.18, "npc"),
            y = unit(0.74, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.7,
                col = "#555555"
            )
        ),

        textGrob(
            "Supportive focal SNP",
            x = unit(0.06, "npc"),
            y = unit(0.57, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.8,
                fontface = "bold",
                col = "#444444"
            )
        ),

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.49, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.5,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.45, "npc"),
            y1 = unit(0.53, "npc"),
            gp = gpar(
                col = DR_colour,
                lwd = 3.0
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.26,
                y = 0.49,
                width = 0.20,
                height = 0.14
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.39, "npc"),
            y = unit(0.51, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.7,
                fontface = "bold",
                col = DR_colour
            )
        ),

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.36, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 6.5,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.32, "npc"),
            y1 = unit(0.40, "npc"),
            gp = gpar(
                col = BDR_colour,
                lwd = 3.0
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.23,
                y = 0.36,
                width = 0.14,
                height = 0.10
            )
        ),

        textGrob(
            "+",
            x = unit(0.29, "npc"),
            y = unit(0.36, "npc"),
            gp = gpar(
                fontsize = 7.2,
                fontface = "bold",
                col = "#555555"
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.35,
                y = 0.36,
                width = 0.14,
                height = 0.10
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.45, "npc"),
            y = unit(0.38, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.7,
                fontface = "bold",
                col = BDR_colour
            )
        ),

        textGrob(
            "exact supportive pattern = Left DR/DR, Right B/DR",
            x = unit(0.06, "npc"),
            y = unit(0.24, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        textGrob(
            "*",
            x = unit(0.09, "npc"),
            y = unit(0.15, "npc"),
            gp = gpar(
                fontsize = 11,
                fontface = "bold",
                col = "#333333"
            )
        ),

        textGrob(
            "asterisk = supportive SNP within a previously annotated B/D/R interval",
            x = unit(0.16, "npc"),
            y = unit(0.15, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        textGrob(
            "taller coloured ticks = supportive SNPs with BA > 0.99",
            x = unit(0.06, "npc"),
            y = unit(0.08, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        )
    )
}

add_sidecar <- function(
    plot,
    sidecar_width = unit(2.25, "in")
) {
    g <- ggplotGrob(plot)

    panel_rows <- which(grepl("^panel", g$layout$name))

    if (length(panel_rows) == 0L) {
        stop("Could not identify ggplot panel.")
    }

    panel_l <- min(g$layout$l[panel_rows])
    panel_t <- min(g$layout$t[panel_rows])
    panel_b <- max(g$layout$b[panel_rows])

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
        name = "supportive-focal-sidecar"
    )

    g
}


# ============================================================
# BUILD FIGURE
# ============================================================

plot_main <- build_base_plot() +

    # non-supportive callable focal SNPs in grey
    geom_segment(
        data = non_supportive_left,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = non_supportive_grey,
        linewidth = 0.30,
        alpha = 0.45
    ) +

    geom_segment(
        data = non_supportive_right,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = non_supportive_grey,
        linewidth = 0.30,
        alpha = 0.45
    ) +

    # all supportive exact positions, faintly
    geom_segment(
        data = supportive_all_left,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = DR_colour,
        linewidth = 0.38,
        alpha = 0.35
    ) +

    geom_segment(
        data = supportive_all_right,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = BDR_colour,
        linewidth = 0.38,
        alpha = 0.35
    ) +

    # high-BA supportive subset, strongly
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
        data = high_right,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = y0,
            yend = y1
        ),
        inherit.aes = FALSE,
        colour = BDR_colour,
        linewidth = 0.95,
        alpha = 1
    ) +

    # asterisks for supportive sites inside published intervals
    geom_text(
        data = pubstar_df,
        aes(
            x = position_Mb,
            y = y,
            label = "*"
        ),
        inherit.aes = FALSE,
        colour = "#333333",
        size = 3.0,
        fontface = "bold",
        vjust = 0.5
    )

grob_main <- add_sidecar(plot_main)


# ============================================================
# SAVE
# ============================================================

draw_pdf(
    grob_main,
    pdf_file,
    width = 9.15,
    height = 5.20
)

draw_png(
    grob_main,
    png_file,
    width = 9.15,
    height = 5.20
)


# ============================================================
# REPORT
# ============================================================

cat(
    "\n========================================\n",
    "ANALYSIS 171 COMPLETE\n",
    "========================================\n\n",
    sep = ""
)

cat(
    "Lower-panel focal logic:\n",
    "  supportive     = Left DR/DR, Right B/DR\n",
    "  non-supportive = all other callable focal L/R SNPs\n\n",
    sep = ""
)

print(summary_table, row.names = FALSE)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nSummary:\n",
    summary_file,
    "\n\nAudit:\n",
    audit_file,
    "\n",
    sep = ""
)
